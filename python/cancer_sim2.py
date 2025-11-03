"""
Stage I Cancer Model — Multi-Therapy, Surgery@t0, RT Repop Delay, Clonogens, Downloads
--------------------------------------------------------------------------------------
What this script does
- Models Stage I tumor growth (Gompertz) with multiple therapy schedules:
  * Pulse chemo: q3w x 6 (cisplatin-like)
  * Metronomic chemo: daily x 6 weeks (dose-matched)
  * Weekly dose-dense chemo: q1w x 8 (dose-matched)
  * External-beam radiotherapy: 50 Gy in 25 fractions (5 weeks, weekdays)
- Adds SURGERY at t=0 (debulking by a configurable fraction).
- Adds RT REPOPULATION DELAY: during RT, proliferation stays suppressed until Tk after RT start,
  then ramps up to baseline over a configurable number of days.
- Tracks total tumor cells AND clonogens (fixed fraction), with a 1-clonogen threshold line.
- Saves each therapy graph as a separate PNG in ./outputs and writes a compact JSON with series + config.

Notes on modeling choices (aligned with literature)
- Growth law: Gompertz often captures decelerating macroscopic growth in solid tumors.
- Chemo: One-compartment PK → Emax PD, with Norton–Simon (kill ∝ untreated growth rate) to compare schedules.
- RT: Linear–Quadratic per fraction for 1–3 Gy/fx; repopulation “kick-off” modeled via Tk + ramp on rho.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Callable, List, Tuple, Optional, Dict
import numpy as np
import matplotlib.pyplot as plt
import json
import os

# =========================
# ===== Data Classes ======
# =========================

@dataclass
class GrowthParams:
    law: str = "gompertz"
    K: float = 1e6
    rho: float = 0.16  # intrinsic/day

@dataclass
class PKParams:
    kel: float
    F: float
    V: float
    doses: List[Tuple[float, float]]  # (day, dose)

@dataclass
class PDParams:
    Emax: float
    EC50: float

@dataclass
class ChemoParams:
    pk: PKParams
    pd: PDParams
    intensity: float = 1.0  # optional scalar

@dataclass
class RTParams:
    alpha: float
    beta: float
    d: float
    schedule: Callable[[float], bool]     # weekday fx in [start, end)
    start: float                          # RT course start (day)
    end: float                            # RT course end (day)
    Tk: float = 28.0                      # repopulation kick-off (days after RT start)
    ramp: float = 7.0                     # days to ramp back to baseline rho
    rho_scale_preTk: float = 0.0          # rho multiplier before Tk (0.0 = full arrest)

@dataclass
class StageParams:
    growth: GrowthParams
    alpha_beta: Tuple[float, float]

@dataclass
class InitialState:
    p0: float


# =========================
# ===== Helper Functions ==
# =========================

def weekday_schedule(start: float, end: float) -> Callable[[float], bool]:
    """True on weekdays in [start, end). Day 0 is treated as weekday index 0."""
    def sched(t: float) -> bool:
        if not (start <= t < end):
            return False
        dow = int(t) % 7
        return dow not in (5, 6)  # 0..4 weekdays; 5,6 weekend
    return sched

def make_chemo_doses(start_day: float, num_cycles: int, cycle_len: float, dose_per_cycle: float):
    return [(start_day + i * cycle_len, dose_per_cycle) for i in range(num_cycles)]

def growth_rate(p: float, gp: GrowthParams) -> float:
    if p <= 0:
        return 0.0
    if gp.law.lower() == "gompertz":
        return gp.rho * p * np.log(max(gp.K / p, 1e-12))
    # logistic fallback
    return gp.rho * p * (1.0 - p / gp.K)

def rho_with_repopulation_delay(t: float, rho_base: float, rt: Optional[RTParams]) -> float:
    """
    During RT: for t in [start, end), suppress growth until Tk after RT start,
    then ramp linearly to rho_base over rt.ramp days.
    Outside RT window: use rho_base.
    """
    if rt is None:
        return rho_base
    if t < rt.start or t >= rt.end:
        # outside the course -> baseline rho
        return rho_base
    # within RT course
    t_since_start = t - rt.start
    if t_since_start < rt.Tk:
        return rt.rho_scale_preTk * rho_base
    # ramp from Tk to Tk + ramp
    ramp_pos = (t_since_start - rt.Tk) / max(rt.ramp, 1e-9)
    ramp_pos = np.clip(ramp_pos, 0.0, 1.0)
    return rho_base * (rt.rho_scale_preTk + (1.0 - rt.rho_scale_preTk) * ramp_pos)

def pk_concentration(t: float, pk: PKParams) -> float:
    C = 0.0
    for tj, Dj in pk.doses:
        if t >= tj:
            C += pk.F * Dj / pk.V * np.exp(-pk.kel * (t - tj))
    return C

def pd_emax(C: float, pd: PDParams) -> float:
    return pd.Emax * C / (pd.EC50 + C)

def ns_kill(growth_untr: float, effect: float, intensity: float) -> float:
    """Norton–Simon kill: proportional to untreated growth rate."""
    return intensity * effect * growth_untr

def apply_rt_LQ(x: float, alpha: float, beta: float, d: float) -> float:
    """One fraction via LQ: S = exp(-αd - βd²)."""
    S = np.exp(-(alpha * d + beta * d * d))
    return x * S


# =========================
# ===== Simulation ========
# =========================

def run_sim(days: np.ndarray,
            init: InitialState,
            stage: StageParams,
            surgery_fraction: float = 0.0,  # fraction removed at t=0 (e.g., 0.95)
            chemo: Optional[ChemoParams] = None,
            rt: Optional[RTParams] = None) -> np.ndarray:
    """
    Returns an array P(t) of tumor cell counts.
    - Surgery at t=0 removes 'surgery_fraction' of initial mass.
    - During RT, growth rate rho is modified by repopulation delay (Tk + ramp).
    """
    # Apply surgery at t=0
    p = init.p0 * (1.0 - np.clip(surgery_fraction, 0.0, 1.0))
    P = [p]

    # Copy growth params (we'll vary rho dynamically)
    gp_base = stage.growth

    for i in range(1, len(days)):
        t_prev = days[i - 1]
        dt = days[i] - days[i - 1]

        # Dynamic rho with repopulation delay (only during RT window)
        rho_eff = rho_with_repopulation_delay(t_prev, gp_base.rho, rt)
        gp_dyn = GrowthParams(law=gp_base.law, K=gp_base.K, rho=rho_eff)

        # Untreated growth at this instant
        g_p = growth_rate(p, gp_dyn)

        # Chemo kill (Norton–Simon: proportional to the untreated growth)
        kill_p = 0.0
        if chemo is not None:
            C = pk_concentration(t_prev, chemo.pk)
            eff = pd_emax(C, chemo.pd)
            kill_p = ns_kill(g_p, eff, chemo.intensity)

        # Continuous update
        p = max(p + (g_p - kill_p) * dt, 0.0)

        # Discrete RT fraction at end of step if scheduled for t_prev
        if rt is not None and rt.schedule(t_prev):
            alpha, beta = stage.alpha_beta
            p = apply_rt_LQ(p, alpha, beta, rt.d)

        P.append(p)

    return np.array(P)


# =========================
# ===== Visualization =====
# =========================

def plot_and_save(days: np.ndarray,
                  P: np.ndarray,
                  title: str,
                  outdir: str,
                  filename: str,
                  f_clon: float = 0.01,
                  chemo: Optional[ChemoParams] = None,
                  rt: Optional[RTParams] = None,
                  surgery_fraction: float = 0.0):
    """Minimal, clean plot with clonogens, 1-clonogen line, and treatment markers."""
    os.makedirs(outdir, exist_ok=True)
    N_clon = f_clon * P

    plt.figure(figsize=(9, 5))
    plt.plot(days, P, linewidth=2.0, label="Tumor cells")
    plt.plot(days, N_clon, linestyle="--", linewidth=2.0, label=f"Clonogens ({int(f_clon*100)}%)")
    plt.axhline(1.0, linestyle=":", linewidth=1.5, label="1 clonogen threshold")

    # Annotate <1 clonogen crossing
    below = N_clon < 1.0
    if np.any(below):
        t_cross = days[np.argmax(below)]
        plt.axvline(t_cross, linestyle="--", linewidth=1.0)
        plt.text(t_cross + 1, 2, f"<1 clonogen @ day {t_cross:.0f}", fontsize=8)

    # Treatment markers
    if chemo is not None:
        for t, _ in chemo.pk.doses:
            plt.axvline(t, linestyle=":", linewidth=0.8, alpha=0.7)

    if rt is not None:
        # Shade 5-day treatment weeks from rt.start to rt.end
        s = rt.start
        while s < rt.end:
            plt.axvspan(s, min(s + 5.0, rt.end), alpha=0.08)
            s += 7.0

    # Surgery marker at t=0 if any
    if surgery_fraction > 0:
        plt.axvline(0.0, linestyle="-.", linewidth=1.0)
        plt.text(0.5, 2, f"Surgery @ t=0 (removed {int(surgery_fraction*100)}%)", fontsize=8)

    plt.yscale("log")
    plt.xlabel("Days since diagnosis")
    plt.ylabel("Cell count (log scale)")
    plt.title(title)
    plt.grid(True, which="both", linestyle=":")
    plt.legend()
    outpath = os.path.join(outdir, filename)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    return outpath


# =========================
# ====== Main Config ======
# =========================

if __name__ == "__main__":
    # Output directory
    OUTDIR = "outputs"
    os.makedirs(OUTDIR, exist_ok=True)

    # Time grid — daily sampling for clarity
    DT = 1.0
    TMAX = 180.0
    DAYS = np.arange(0.0, TMAX + DT, DT)

    # Stage I parameters and initial state
    STAGE1 = StageParams(growth=GrowthParams("gompertz", K=1.0e6, rho=0.16),
                         alpha_beta=(0.30, 0.03))
    INIT = InitialState(p0=1.2e4)

    # SURGERY at t=0 — set how much to remove (e.g., 0.95 = 95% debulking)
    SURGERY_FRAC = 0.95

    # PK/PD defaults (cisplatin-like)
    Emax = 0.9
    EC50 = 5.0
    kel = np.log(2) / 1.0  # half-life ~1 day
    Vapp = 20.0
    Fbio = 1.0

    # ===== Chemo regimens =====
    # 1) Pulse q3w x 6 (cisplatin 75 mg/m^2, BSA 1.8 => ~135 mg/cycle)
    dose_cycle = 75.0 * 1.8
    doses_pulse = make_chemo_doses(0.0, 6, 21.0, dose_cycle)
    CHEMO_PULSE = ChemoParams(
        pk=PKParams(kel=kel, F=Fbio, V=Vapp, doses=doses_pulse),
        pd=PDParams(Emax=Emax, EC50=EC50),
        intensity=1.0
    )

    # 2) Metronomic daily x 6 wk (dose-matched to total pulse)
    total_pulse = dose_cycle * 6
    daily_metro = total_pulse / 42.0
    doses_metro = [(float(t), daily_metro) for t in np.arange(0.0, 42.0, 1.0)]
    CHEMO_METRO = ChemoParams(
        pk=PKParams(kel=kel, F=Fbio, V=Vapp, doses=doses_metro),
        pd=PDParams(Emax=Emax, EC50=EC50),
        intensity=1.0
    )

    # 3) Weekly dose-dense q1w x 8 (dose-matched to total pulse)
    dose_weekly = total_pulse / 8.0
    doses_weekly = [(7.0 * i, dose_weekly) for i in range(8)]
    CHEMO_WEEKLY = ChemoParams(
        pk=PKParams(kel=kel, F=Fbio, V=Vapp, doses=doses_weekly),
        pd=PDParams(Emax=Emax, EC50=EC50),
        intensity=1.0
    )

    # ===== Radiotherapy params =====
    RT_START, RT_END = 0.0, 35.0  # 5 weeks
    RT = RTParams(
        alpha=0.30, beta=0.03, d=2.0,
        schedule=weekday_schedule(RT_START, RT_END),
        start=RT_START, end=RT_END,
        Tk=28.0,       # kick-off after 4 weeks from RT start
        ramp=7.0,      # 1-week ramp back to baseline
        rho_scale_preTk=0.0  # full arrest before Tk
    )

    # =========================
    # ==== Run simulations ====
    # =========================

    P_pulse = run_sim(DAYS, INIT, STAGE1, SURGERY_FRAC, chemo=CHEMO_PULSE, rt=None)
    P_metro = run_sim(DAYS, INIT, STAGE1, SURGERY_FRAC, chemo=CHEMO_METRO, rt=None)
    P_weekly = run_sim(DAYS, INIT, STAGE1, SURGERY_FRAC, chemo=CHEMO_WEEKLY, rt=None)
    P_rt = run_sim(DAYS, INIT, STAGE1, SURGERY_FRAC, chemo=None, rt=RT)

    # =========================
    # ===== Save figures  =====
    # =========================

    fp_pulse = plot_and_save(
        DAYS, P_pulse,
        title="Pulse Chemo (q3w ×6) — Tumor & Clonogens (Surgery@t0, RT repop delay off)",
        outdir=OUTDIR, filename="pulse_chemo.png",
        chemo=CHEMO_PULSE, rt=None, surgery_fraction=SURGERY_FRAC
    )
    fp_metro = plot_and_save(
        DAYS, P_metro,
        title="Metronomic Chemo (daily ×6wk, dose-matched) — Tumor & Clonogens (Surgery@t0)",
        outdir=OUTDIR, filename="metronomic_chemo.png",
        chemo=CHEMO_METRO, rt=None, surgery_fraction=SURGERY_FRAC
    )
    fp_weekly = plot_and_save(
        DAYS, P_weekly,
        title="Weekly Dose-Dense Chemo (q1w ×8) — Tumor & Clonogens (Surgery@t0)",
        outdir=OUTDIR, filename="weekly_chemo.png",
        chemo=CHEMO_WEEKLY, rt=None, surgery_fraction=SURGERY_FRAC
    )
    fp_rt = plot_and_save(
        DAYS, P_rt,
        title="Radiotherapy (50 Gy / 25 fx) — Tumor & Clonogens (Surgery@t0, Tk+Ramp)",
        outdir=OUTDIR, filename="radiotherapy.png",
        chemo=None, rt=RT, surgery_fraction=SURGERY_FRAC
    )

    # =========================
    # ===== Export JSON  ======
    # =========================

    results = {
        "days": DAYS.tolist(),
        "tumor": {
            "chemo_pulse": P_pulse.tolist(),
            "chemo_metronomic": P_metro.tolist(),
            "chemo_weekly": P_weekly.tolist(),
            "radiotherapy": P_rt.tolist()
        },
        "clonogens_fraction": 0.01,
        "surgery_fraction_t0": SURGERY_FRAC,
        "config": {
            "growth": asdict(STAGE1.growth),
            "rt": {
                "alpha": RT.alpha, "beta": RT.beta, "d": RT.d,
                "start": RT.start, "end": RT.end,
                "Tk": RT.Tk, "ramp": RT.ramp, "rho_scale_preTk": RT.rho_scale_preTk
            },
            "chemo": {
                "pulse": {"pk": asdict(CHEMO_PULSE.pk), "pd": asdict(CHEMO_PULSE.pd)},
                "metronomic": {"pk": asdict(CHEMO_METRO.pk), "pd": asdict(CHEMO_METRO.pd)},
                "weekly": {"pk": asdict(CHEMO_WEEKLY.pk), "pd": asdict(CHEMO_WEEKLY.pd)}
            }
        },
        "files": {
            "pulse_chemo_png": fp_pulse,
            "metronomic_chemo_png": fp_metro,
            "weekly_chemo_png": fp_weekly,
            "radiotherapy_png": fp_rt
        }
    }

    with open(os.path.join(OUTDIR, "stage1_surgery_repop_outputs.json"), "w") as f:
        json.dump(results, f, indent=2)

    print("Saved to", OUTDIR)
    for k, v in results["files"].items():
        print("-", k, "->", v)
