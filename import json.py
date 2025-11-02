import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import json

# -----------------------------
# Simulation settings
# -----------------------------
K = 1e6        # carrying capacity
dt = 0.5       # timestep in days
T_max = 730    # simulate 2 years
days = np.arange(0, T_max + dt, dt)
start_date = datetime(2025, 11, 1)
calendar_dates = [start_date + timedelta(days=float(d)) for d in days]

# -----------------------------
# Tumor profiles by stage + grade
# -----------------------------
profiles = {
    "1L": {"p0":1e4, "m0":0, "gamma_p":0.08, "gamma_m":0.05, "k_c":0.05, "k_r":0.05},
    "1H": {"p0":1e4, "m0":0, "gamma_p":0.15, "gamma_m":0.08, "k_c":0.05, "k_r":0.08},
    "2L": {"p0":5e4, "m0":1e3, "gamma_p":0.08, "gamma_m":0.06, "k_c":0.05, "k_r":0.05},
    "2H": {"p0":5e4, "m0":1e3, "gamma_p":0.18, "gamma_m":0.12, "k_c":0.05, "k_r":0.08},
    "3L": {"p0":2e5, "m0":5e3, "gamma_p":0.12, "gamma_m":0.08, "k_c":0.05, "k_r":0.05},
    "3H": {"p0":2e5, "m0":5e3, "gamma_p":0.25, "gamma_m":0.18, "k_c":0.03, "k_r":0.05},
    "4L": {"p0":5e5, "m0":1e4, "gamma_p":0.15, "gamma_m":0.10, "k_c":0.03, "k_r":0.05},
    "4H": {"p0":5e5, "m0":1e4, "gamma_p":0.30, "gamma_m":0.20, "k_c":0.02, "k_r":0.04}
}

# -----------------------------
# Lifestyle modifiers
# -----------------------------
lifestyles = {"H":1.05, "A":1.0, "P":0.95}

# -----------------------------
# User inputs
# -----------------------------
profile_choice = input("Choose profile (Stage+Grade, e.g., 1L, 2H): ").strip().upper()
profile = profiles.get(profile_choice, profiles["1L"])

lifestyle_choice = input("Lifestyle (H=Healthy, A=Average, P=Poor): ").strip().upper()
modifier = lifestyles.get(lifestyle_choice, 1.0)

def get_int_input(prompt, min_val, max_val):
    while True:
        try:
            val = int(input(prompt))
            return min(max(val, min_val), max_val)
        except ValueError:
            print(f"Enter an integer between {min_val} and {max_val}")

n_chemo = get_int_input("Number of chemo cycles (1–8): ", 1, 8)
chemo_len = get_int_input("Chemo length per cycle (14–28 days): ", 14, 28)
radio_len = get_int_input("Radiotherapy duration (35–56 days): ", 35, 56)

# -----------------------------
# Generate therapy intervals
# -----------------------------
chemo_cycles = []
start_day = 20
for _ in range(n_chemo):
    end_day = start_day + chemo_len
    chemo_cycles.append((start_day, end_day))
    start_day = end_day + 7  # 1-week rest

radio_start = chemo_cycles[-1][1] + 14
radio_cycles = [(radio_start, radio_start + radio_len)]

def chemo_func(day):
    return any(start <= day < end for start, end in chemo_cycles)

def radio_func(day):
    return any(start <= day < end for start, end in radio_cycles)

# -----------------------------
# Run stochastic simulations
# -----------------------------
n_simulations = 50
primary_matrix = np.zeros((n_simulations, len(days)))
meta_matrix = np.zeros((n_simulations, len(days)))

for sim in range(n_simulations):
    rng = np.random.default_rng(seed=sim)
    p, m = profile["p0"], profile["m0"]
    gamma_p, gamma_m = profile["gamma_p"], profile["gamma_m"]
    k_c, k_r = profile["k_c"]*modifier, profile["k_r"]*modifier
    
    primary_sim = []
    meta_sim = []

    for day in days:
        gp = gamma_p * rng.uniform(0.95, 1.05)
        gm = gamma_m * rng.uniform(0.95, 1.05)
        kc = k_c * rng.uniform(0.95, 1.05)
        kr = k_r * rng.uniform(0.95, 1.05)

        dp = gp*p*(1 - p/K) - kc*chemo_func(day)*p - kr*radio_func(day)*p
        dm = 1e-6*p**0.7 + gm*m*(1 - m/K) - kc*chemo_func(day)*m

        p += dp*dt
        m += dm*dt

        # Cure threshold
        if p < 10: p = 0
        if m < 10: m = 0

        primary_sim.append(p)
        meta_sim.append(m)

    primary_matrix[sim, :] = primary_sim
    meta_matrix[sim, :] = meta_sim

# -----------------------------
# Median and 95% interval
# -----------------------------
median_primary = np.median(primary_matrix, axis=0)
lower_primary = np.percentile(primary_matrix, 2.5, axis=0)
upper_primary = np.percentile(primary_matrix, 97.5, axis=0)

median_meta = np.median(meta_matrix, axis=0)
lower_meta = np.percentile(meta_matrix, 2.5, axis=0)
upper_meta = np.percentile(meta_matrix, 97.5, axis=0)

# -----------------------------
# Save JSON
# -----------------------------
data_out = [
    {
        "day": int(days[i]),
        "primary_median": float(median_primary[i]),
        "primary_lower": float(lower_primary[i]),
        "primary_upper": float(upper_primary[i]),
        "meta_median": float(median_meta[i]),
        "meta_lower": float(lower_meta[i]),
        "meta_upper": float(upper_meta[i])
    }
    for i in range(len(days))
]

with open("cancer_progression.json", "w") as f:
    json.dump(data_out, f, indent=4)

print("Simulation complete. Output saved to cancer_progression.json")

# -----------------------------
# Plot results
# -----------------------------
plt.figure(figsize=(12,6))
plt.plot(calendar_dates, median_primary, color='red', label="Primary Tumor (median)")
plt.fill_between(calendar_dates, lower_primary, upper_primary, color='red', alpha=0.2, label="Primary Tumor ±95%")
plt.plot(calendar_dates, median_meta, color='blue', label="Metastases (median)")
plt.fill_between(calendar_dates, lower_meta, upper_meta, color='blue', alpha=0.2, label="Metastases ±95%")

# Therapy spans
for idx, (start, end) in enumerate(chemo_cycles):
    lbl = "Chemo" if idx == 0 else None
    plt.axvspan(calendar_dates[int(start/dt)], calendar_dates[min(int(end/dt), len(calendar_dates)-1)], color='orange', alpha=0.2, label=lbl)
for idx, (start, end) in enumerate(radio_cycles):
    lbl = "Radiotherapy" if idx == 0 else None
    plt.axvspan(calendar_dates[int(start/dt)], calendar_dates[min(int(end/dt), len(calendar_dates)-1)], color='green', alpha=0.2, label=lbl)

plt.xlabel("Date")
plt.ylabel("Tumor Cells")
plt.title(f"Cancer Progression ({profile_choice}, {lifestyle_choice} lifestyle)")
plt.grid(True)
plt.xticks(rotation=30)
plt.legend()
plt.tight_layout()
plt.show()
