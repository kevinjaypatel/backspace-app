import numpy as np
import matplotlib.pyplot as plt

# --- Parameters & distributions ---
n_runs = 1000
dt = 0.1
t_final = 140
t = np.arange(0, t_final + dt, dt)

# Distributions
# Initial clonogenic cell number ~ log-normal
logN0_mean = 10      # log10 scale ~ 10 => ~1e10 clonogens
logN0_std  = 0.5
# Growth rate r ~ normal
r_mean = 0.1
r_std = 0.02
# Carrying capacity K ~ lognormal
logK_mean = 9       # ~1e9
logK_std  = 0.3
# Treatment schedule
treatment_days = [20, 40, 60, 80, 100]
# Kill‐fraction ~ Beta distribution (to keep within 0–1)
kill_shape_a = 9
kill_shape_b = 1

# Threshold for “controlled” in clonogen count
threshold = 1.0

# Storage
survivor_counts = np.zeros((n_runs, len(t)))
tcp_empirical   = np.zeros((len(t),), dtype=float)

for run in range(n_runs):
    # sample parameters for this run
    N0 = 10**(np.random.normal(logN0_mean, logN0_std))
    r  = np.random.normal(r_mean, r_std)
    K  = 10**(np.random.normal(logK_mean, logK_std))
    kill_frac = np.random.beta(kill_shape_a, kill_shape_b)
    
    N = np.zeros_like(t)
    N[0] = N0
    for i in range(1, len(t)):
        # deterministic logistic increment as approximation
        N[i] = N[i-1] + r * N[i-1] * (1 - N[i-1]/K) * dt
        
        # apply treatment if scheduled
        if any(abs(t[i] - td) < dt/2 for td in treatment_days):
            N[i] *= (1 - kill_frac)
        
        # cannot go below zero
        if N[i] < 0:
            N[i] = 0
    
    survivor_counts[run, :] = N

# Compute empirical TCP: fraction of runs where count < threshold
for i, ti in enumerate(t):
    tcp_empirical[i] = np.mean(survivor_counts[:, i] < threshold)

# Plotting
plt.figure(figsize=(12, 6))
# plot a subset of runs
for run in range(min(50, n_runs)):
    plt.plot(t, survivor_counts[run, :], color='gray', alpha=0.3)
plt.axhline(y=threshold, color='red', linestyle=':', label='Control threshold')
plt.title('Stochastic tumour‐growth + treatment model')
plt.xlabel('Time (days)')
plt.ylabel('Clonogenic cell count (log scale)')
plt.yscale('log')
plt.grid(True, which='both', ls=':')
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 4))
plt.plot(t, tcp_empirical, color='blue')
plt.title('Empirical Tumour Control Probability (TCP) vs Time')
plt.xlabel('Time (days)')
plt.ylabel('TCP (fraction of runs controlled)')
plt.grid(True)
plt.tight_layout()
plt.show()
