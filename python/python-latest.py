import numpy as np
import matplotlib.pyplot as plt
import os
import json

# -----------------------------
# Model Parameters
# -----------------------------
N0_tumor = 1e8         # initial tumor cell number
clonogenic_fraction = 1e-5  # fraction of clonogenic cells within tumor
r_mean, r_std = 0.1, 0.02   # growth rate/day (mean ± SD)
K_mean, K_std = 1e10, 1e9   # carrying capacity ± variance
kill_mean, kill_std = 0.98, 0.03  # treatment kill fraction ± SD
treat_days = [20, 40, 60, 80, 100, 120]
dt, t_final = 0.1, 160
n_runs = 500
threshold = 1.0  # control threshold for clonogens

t = np.arange(0, t_final + dt, dt)

# -----------------------------
# Helper Functions
# -----------------------------
def simulate_tumor():
    """Simulate one stochastic tumor + clonogen trajectory."""
    r = np.random.normal(r_mean, r_std)
    K = np.random.normal(K_mean, K_std)
    kill = np.clip(np.random.normal(kill_mean, kill_std), 0, 1)
    
    N_tumor = np.zeros_like(t)
    N_clono = np.zeros_like(t)
    N_tumor[0] = N0_tumor
    N_clono[0] = N0_tumor * clonogenic_fraction

    for i in range(1, len(t)):
        # Logistic growth
        N_tumor[i] = N_tumor[i-1] + r * N_tumor[i-1] * (1 - N_tumor[i-1]/K) * dt
        N_clono[i] = N_clono[i-1] + r * N_clono[i-1] * (1 - N_clono[i-1]/K) * dt
        
        # Apply treatment on scheduled days
        if any(abs(t[i] - td) < dt/2 for td in treat_days):
            N_tumor[i] *= (1 - kill)
            N_clono[i] *= (1 - kill)
        
        # prevent negative or zero for log scale
        N_tumor[i] = max(N_tumor[i], 1)
        N_clono[i] = max(N_clono[i], 1)
    return N_tumor, N_clono

# -----------------------------
# Monte Carlo Simulation
# -----------------------------
tumor_mat = np.zeros((n_runs, len(t)))
clono_mat = np.zeros((n_runs, len(t)))

for run in range(n_runs):
    tumor_mat[run, :], clono_mat[run, :] = simulate_tumor()

# Compute statistics
tumor_median = np.median(tumor_mat, axis=0)
tumor_p25, tumor_p75 = np.percentile(tumor_mat, [25, 75], axis=0)
clono_median = np.median(clono_mat, axis=0)
clono_p25, clono_p75 = np.percentile(clono_mat, [25, 75], axis=0)

# -----------------------------
# Plotting
# -----------------------------
fig, ax = plt.subplots(figsize=(10, 6))

# Tumor total
ax.plot(t, tumor_median, color='tab:blue', label='Tumor (median)')
ax.fill_between(t, tumor_p25, tumor_p75, color='tab:blue', alpha=0.2, label='Tumor IQR')

# Clonogenic cells
ax.plot(t, clono_median, color='tab:green', label='Clonogenic cells (median)')
ax.fill_between(t, clono_p25, clono_p75, color='tab:green', alpha=0.2, label='Clonogenic IQR')

# Threshold line
ax.axhline(y=threshold, color='red', linestyle='--', label='Control Threshold (1 clonogen)')

# Treatment markers
for td in treat_days:
    ax.axvline(x=td, color='grey', linestyle=':', alpha=0.5)

# Scientific formatting
ax.set_yscale('log')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Cell number (log scale)')
ax.set_title('Tumor and Clonogenic Cell Dynamics with Stochastic Treatments')
ax.legend()
ax.grid(True, which='both', ls=':')
plt.tight_layout()

# Save plot to outputs folder
output_dir = os.path.join(os.path.dirname(__file__), '..', 'outputs')
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, 'python_latest_simulation.png')
plt.savefig(output_path, dpi=150)
plt.close()
print(f"Plot saved to: {output_path}")

# Save JSON data for frontend
json_data = {
    "days": t.tolist(),
    "tumor": {
        "median": tumor_median.tolist(),
        "p25": tumor_p25.tolist(),
        "p75": tumor_p75.tolist()
    },
    "clonogens": {
        "median": clono_median.tolist(),
        "p25": clono_p25.tolist(),
        "p75": clono_p75.tolist()
    },
    "threshold": threshold,
    "treatment_days": treat_days,
    "parameters": {
        "N0_tumor": N0_tumor,
        "clonogenic_fraction": clonogenic_fraction,
        "r_mean": r_mean,
        "r_std": r_std,
        "K_mean": K_mean,
        "K_std": K_std,
        "kill_mean": kill_mean,
        "kill_std": kill_std,
        "n_runs": n_runs
    }
}

json_path = os.path.join(output_dir, 'python_latest_simulation.json')
with open(json_path, 'w') as f:
    json.dump(json_data, f, indent=2)
print(f"JSON data saved to: {json_path}")
