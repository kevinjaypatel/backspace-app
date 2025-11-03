"""
Cancer progression simulator with therapy and parameter variability
==================================================================

What this code does (big picture)
---------------------------------
- Simulates primary tumor and total metastatic burden over time.
- Uses a biologically motivated growth law (Gompertz by default).
- Applies chemotherapy using the Norton–Simon hypothesis (kill proportional
  to the *unperturbed* tumor growth rate) [Ref. 2–5].
- Applies radiotherapy with the Linear–Quadratic (LQ) model per fraction,
  a clinical standard for cell survival after ionizing radiation [Ref. 6–8].
- Allows Monte Carlo variability across key parameters to visualize uncertainty.
- Summarizes uncertainty via 5–95% envelopes and medians.

Key modeling choices (and why)
------------------------------
1) Growth law = Gompertz (optionally Logistic).
   Gompertz often describes decelerating tumor growth better than Logistic in
   many settings and is widely used in oncology modeling [Ref. 1,9–11].
2) Chemo kill = Norton–Simon hypothesis.
   Drug effect is proportional to the instantaneous *untreated* growth rate;
   this usually better captures response than constant-kill models [Ref. 2–5].
3) Radiotherapy = Linear–Quadratic model.
   Per fraction: survival S = exp(-α d - β d^2). This is the default radiobiology
   workhorse for conventional fractionation and often adequate for daily 1–3 Gy
   dosing; be cautious for very large fractions (SBRT), where extensions exist
   [Ref. 6–8].
4) Metastatic emission = λ0 * p**beta.
   The exponent beta=2/3 is often used as a surface-area scaling ansatz, but the
   exponent is context-specific and remains tunable here [Ref. 12–13].

References (inline numbers used throughout)
-------------------------------------------
[1] “The Model Muddle: In Search of Tumor Growth Laws”, Cancer Research (2013).  # Ref. 1
[2] West et al., "Chemotherapeutic Dose Scheduling Based on Tumor Growth..." Cancer Research (2017).  # Ref. 2
[3] López et al., “The dose-dense principle in chemotherapy” (arXiv, 2017).  # Ref. 3
[4] López & Ventura, "Modelling the Norton-Simon hypothesis", Commun Nonlinear Sci (2019).  # Ref. 4
[5] Review notes on Norton–Simon hypothesis (multiple sources).  # Ref. 5
[6] Brenner, “Point: The linear–quadratic model is appropriate...”, Int J Radiat Oncol (2008).  # Ref. 6
[7] Miyakawa et al., “Applicability of the LQ model...”, J Radiat Res (2013).  # Ref. 7
[8] van Leeuwen et al., “The alpha and beta of tumours: a review...”, Radiat Oncol (2018).  # Ref. 8
[9] Tjørve & Tjørve, “The use of Gompertz models in growth analyses”, PLoS One (2017).  # Ref. 9
[10] Vaghi et al., “Population modeling of tumor growth curves...”, CPT PSP (2020).  # Ref. 10
[11] Heesterman et al., “Mathematical Models for Tumor Growth...”, Cancers (2018).  # Ref. 11
[12] Franssen et al., “A Mathematical Framework for Modelling the Metastatic Process”, J Theor Biol (2019).  # Ref. 12
[13] Survey statements that metastatic emission often scales ~ p^(2/3).  # Ref. 13

NOTE: See bottom of file for parameter guidance and tips.
"""

import json
from dataclasses import dataclass
from typing import Callable, Dict, Tuple, Literal, Optional

import numpy as np
import matplotlib.pyplot as plt

# =========================
# ====== Parameters =======
# =========================

@dataclass
class GrowthParams:
    K: float                   # carrying capacity (cells)
    rho: float                 # intrinsic growth rate (1/day) for Gompertz/Logistic
    law: Literal["gompertz", "logistic"] = "gompertz"  # choose growth law
    # [Ref. 1,9–11] Gompertz widely used for tumor growth; logistic optional.


@dataclass
class ChemoParams:
    k_NS: float                # proportionality for Norton–Simon kill (unitless scaling)
    schedule: Callable[[float], float]  # returns intensity in [0,1] at a given day
    # [Ref. 2–5] kill = k_NS * g_untreated(p). If schedule=0, no chemo that day.


@dataclass
class RTParams:
    alpha: float               # alpha (Gy^-1)
    beta: float                # beta  (Gy^-2)
    dose_per_fraction: float   # d (Gy per fraction)
    schedule: Callable[[float], bool]  # True if a fraction is delivered that day
    # [Ref. 6–8] Survival S = exp(-alpha*d - beta*d^2) applied multiplicatively.


@dataclass
class MetastasisParams:
    lambda0: float             # baseline emission rate
    beta_exp: float            # exponent in p**beta (often ~2/3 as a surface proxy) [Ref. 12–13]
    # m growth uses same form/parameters as primary but can differ if desired:
    gamma_m: float             # metastatic growth rate parameter (1/day)


@dataclass
class InitialState:
    p0: float                  # primary initial cells
    m0: float                  # total metastases initial cells


@dataclass
class SimConfig:
    dt: float
    T_max: float
    n_sims: int
    seed: Optional[int] = 123


# =========================
# == Therapy Schedules  ==
# =========================

def rectangular_schedule(start: float, end: float, intensity: float = 1.0) -> Callable[[float], float]:
    """
    Rectangular (on/off) schedule that returns `intensity` during [start, end) and 0 otherwise.
    Use for chemotherapy intensity in [0,1].
    """
    def sched(t: float) -> float:
        return intensity if (start <= t < end) else 0.0
    return sched


def weekday_rt_schedule(start: float, end: float) -> Callable[[float], bool]:
    """
    Radiotherapy schedule delivering one fraction per *calendar day* in [start, end)
    but skipping weekends (Sat/Sun) for conventional fractionation.
    You can turn weekends on by editing this function if needed.
    """
    def sched(t: float) -> bool:
        # Treat t=0 as day 0 (Monday can be assumed arbitrarily); approximate weekday via modulo 7.
        # 0..6 cycle; set 5,6 to weekend off.
        in_window = (start <= t < end)
        day_of_week = int(t) % 7
        is_weekday = day_of_week not in (5, 6)
        return bool(in_window and is_weekday)
    return sched


# =========================
# ====== Core Model  ======
# =========================

def growth_rate(p: float, gp: GrowthParams) -> float:
    """
    Tumor growth law g(p).
    - Gompertz: dp/dt = rho * p * ln(K/p)         [Ref. 1,9–11]
    - Logistic: dp/dt = rho * p * (1 - p/K)
    """
    if p <= 0:
        return 0.0
    if gp.law == "gompertz":
        return gp.rho * p * np.log(max(gp.K / p, 1e-12))
    elif gp.law == "logistic":
        return gp.rho * p * (1.0 - p / gp.K)
    else:
        raise ValueError("Unknown growth law")


def apply_chemo_NS(p: float, g_untreated: float, chemo: ChemoParams) -> float:
    """
    Norton–Simon hypothesis for chemotherapy:
      kill term = k_NS * g_untreated(p) * u(t), where u in [0,1] is schedule intensity.
    This means higher kill when the (hypothetical) untreated growth would be faster [Ref. 2–5].
    Returns the *rate* to subtract from dp/dt (i.e., units of cells/day).
    """
    return chemo.k_NS * g_untreated


def apply_rt_LQ(p: float, rt: RTParams) -> float:
    """
    Apply one radiotherapy fraction via Linear–Quadratic survival:
      S = exp(-alpha*d - beta*d^2)  [Ref. 6–8]
    This is an *instantaneous* multiplicative drop in cell number when a fraction is delivered.
    """
    d = rt.dose_per_fraction
    S = np.exp(-(rt.alpha * d + rt.beta * d * d))
    return p * S


def step(
    p: float, m: float,
    t: float, dt: float,
    gp: GrowthParams,
    chemo: Optional[ChemoParams],
    rt: Optional[RTParams],
    mp: MetastasisParams
) -> Tuple[float, float]:
    """
    One Euler–Maruyama style step (deterministic here; add noise if desired).
    1) Compute untreated growth rates for p and m.
    2) Apply chemotherapy as a reduction in the *rate* (Norton–Simon).
    3) Update p, m continuously over dt.
    4) If RT scheduled, apply a discrete LQ survival hit at the *end* of the day.
    5) Metastatic emission λ0 * p**beta_exp, with independent Gompertz/Logistic growth for m.

    NOTE: If you prefer exact treatment-at-beginning-of-day, swap the order.
    """
    # Untreated growth rates
    g_p = growth_rate(p, gp)
    g_m = mp.gamma_m * m if gp.law == "logistic" else (mp.gamma_m * m * np.log(max(gp.K / max(m, 1e-12), 1e-12)))

    # Chemo intensity
    u = chemo.schedule(t) if chemo else 0.0

    # Norton–Simon kill terms (reduce rates)
    kill_p = u * apply_chemo_NS(p, g_p, chemo) if chemo else 0.0
    kill_m = u * apply_chemo_NS(m, g_m, chemo) if chemo else 0.0  # assume same mechanism on m

    # Net continuous dynamics (Euler)
    dp = (g_p - kill_p) * dt
    dm = (mp.lambda0 * (p ** mp.beta_exp) + (g_m - kill_m)) * dt

    p = max(p + dp, 0.0)
    m = max(m + dm, 0.0)

    # Apply discrete RT at end of step if scheduled
    if rt and rt.schedule(t):
        p = apply_rt_LQ(p, rt)
        m = apply_rt_LQ(m, rt)   # irradiating mets may or may not be clinically intended; included for completeness

    return p, m


# =========================
# ====== Simulator =========
# =========================

def run_sim(
    days: np.ndarray,
    init: InitialState,
    gp: GrowthParams,
    chemo: Optional[ChemoParams],
    rt: Optional[RTParams],
    mp: MetastasisParams
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Run a single trajectory.
    Returns arrays of p(t) and m(t) for t in `days`.
    """
    p, m = init.p0, init.m0
    P, M = [p], [m]
    for i in range(1, len(days)):
        t_prev = days[i-1]
        dt = days[i] - days[i-1]
        p, m = step(p, m, t_prev, dt, gp, chemo, rt, mp)
        P.append(p)
        M.append(m)
    return np.array(P), np.array(M)


def monte_carlo(
    days: np.ndarray,
    init: InitialState,
    gp_base: GrowthParams,
    chemo_base: Optional[ChemoParams],
    rt_base: Optional[RTParams],
    mp_base: MetastasisParams,
    ranges: Dict[str, Tuple[float, float]],
    n_sims: int,
    seed: Optional[int] = 123
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Monte Carlo over parameter ranges (uniform sampling).
    `ranges` may include keys: "rho", "gamma_m", "k_NS", "alpha", "beta", "dose_per_fraction",
                               "lambda0", "beta_exp", "K".
    Returns arrays of shape (n_sims, len(days)) for primary and metastases.
    """
    rng = np.random.default_rng(seed)
    all_P, all_M = [], []
    for _ in range(n_sims):
        # Copy base and perturb
        gp = GrowthParams(
            K = rng.uniform(*ranges.get("K", (gp_base.K, gp_base.K))),
            rho = rng.uniform(*ranges.get("rho", (gp_base.rho, gp_base.rho))),
            law = gp_base.law
        )
        mp = MetastasisParams(
            lambda0 = rng.uniform(*ranges.get("lambda0", (mp_base.lambda0, mp_base.lambda0))),
            beta_exp = rng.uniform(*ranges.get("beta_exp", (mp_base.beta_exp, mp_base.beta_exp))),
            gamma_m = rng.uniform(*ranges.get("gamma_m", (mp_base.gamma_m, mp_base.gamma_m))),
        )
        chemo = None
        if chemo_base:
            chemo = ChemoParams(
                k_NS = rng.uniform(*ranges.get("k_NS", (chemo_base.k_NS, chemo_base.k_NS))),
                schedule = chemo_base.schedule
            )
        rt = None
        if rt_base:
            rt = RTParams(
                alpha = rng.uniform(*ranges.get("alpha", (rt_base.alpha, rt_base.alpha))),
                beta = rng.uniform(*ranges.get("beta", (rt_base.beta, rt_base.beta))),
                dose_per_fraction = rng.uniform(*ranges.get("dose_per_fraction", (rt_base.dose_per_fraction, rt_base.dose_per_fraction))),
                schedule = rt_base.schedule
            )

        P, M = run_sim(days, init, gp, chemo, rt, mp)
        all_P.append(P); all_M.append(M)
    return np.vstack(all_P), np.vstack(all_M)


# =========================
# ====== Defaults =========
# =========================

# Time grid
dt = 0.5
T_max = 200
days = np.arange(0.0, T_max + dt, dt)

# Growth (Gompertz by default) [Ref. 1,9–11]
gp_base = GrowthParams(K=1e6, rho=0.15, law="gompertz")

# Initial sizes
init = InitialState(p0=1e4, m0=0.0)

# Chemo schedule: day 50–100 continuous intensity 1.0 (edit as needed)
chemo_sched = rectangular_schedule(50, 100, intensity=1.0)
# Norton–Simon proportionality [Ref. 2–5]
chemo_base = ChemoParams(k_NS=1.0, schedule=chemo_sched)

# RT schedule: 2 Gy/fr weekdays from day 120–150 [Ref. 6–8]
rt_sched = weekday_rt_schedule(120, 150)
rt_base = RTParams(alpha=0.3, beta=0.03, dose_per_fraction=2.0, schedule=rt_sched)
# NOTE: alpha/beta here implies α/β ≈ 10 Gy, typical for early-responding tissues/many tumors [Ref. 8].

# Metastatic emission & growth [Ref. 12–13]
mp_base = MetastasisParams(lambda0=1e-6, beta_exp=2/3, gamma_m=0.12)

# Monte Carlo ranges (± ~20–40% around base where sensible)
ranges = {
    "rho": (0.10, 0.22),
    "gamma_m": (0.08, 0.18),
    "k_NS": (0.6, 1.2),
    "alpha": (0.25, 0.4),
    "beta": (0.02, 0.05),
    "dose_per_fraction": (1.8, 2.2),
    "lambda0": (1e-7, 5e-6),
    "beta_exp": (0.55, 0.8),
    "K": (8e5, 1.2e6),
}

# =========================
# ====== Run sims =========
# =========================

N_simulations = 200
all_primary, all_metastases = monte_carlo(
    days, init, gp_base, chemo_base, rt_base, mp_base,
    ranges, n_sims=N_simulations, seed=123
)

# Summaries
p_low, p_med, p_high = np.percentile(all_primary, [5, 50, 95], axis=0)
m_low, m_med, m_high = np.percentile(all_metastases, [5, 50, 95], axis=0)

# =========================
# ====== Plotting =========
# =========================

plt.figure(figsize=(11, 6))
plt.fill_between(days, p_low, p_high, alpha=0.2, label="Primary 5–95%")
plt.plot(days, p_med, lw=2, label="Primary median")
plt.fill_between(days, m_low, m_high, alpha=0.2, label="Metastases 5–95%")
plt.plot(days, m_med, lw=2, label="Metastases median")

# Shade therapy windows
plt.axvspan(50, 100, alpha=0.15, label="Chemo window")
plt.axvspan(120, 150, alpha=0.15, label="RT window (wkdays)")

plt.yscale("log")  # tumor burdens span orders of magnitude
plt.xlabel("Days")
plt.ylabel("Cells (log scale)")
plt.title("Cancer Progression with Norton–Simon Chemo and LQ Radiotherapy (Gompertz growth)")
plt.grid(True, which="both", ls=":")
plt.legend()
plt.tight_layout()
plt.show()

# =========================
# === JSON export (opt) ===
# =========================
results = {
    "days": days.tolist(),
    "primary_median": p_med.tolist(),
    "primary_5": p_low.tolist(),
    "primary_95": p_high.tolist(),
    "met_median": m_med.tolist(),
    "met_5": m_low.tolist(),
    "met_95": m_high.tolist(),
    "config": {
        "growth": gp_base.__dict__,
        "chemo": {"k_NS": chemo_base.k_NS} if chemo_base else None,
        "rt": {"alpha": rt_base.alpha, "beta": rt_base.beta, "dose_per_fraction": rt_base.dose_per_fraction} if rt_base else None,
        "metastasis": mp_base.__dict__
    }
}
print(json.dumps(results)[:250] + " ...")  # preview
