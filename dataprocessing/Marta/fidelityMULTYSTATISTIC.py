import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit
from math import erf, sqrt

# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v28_10_11_25.db"
run_first, run_last = 480, 499        # range of runs to aggregate
bins      = 200                      # histogram bins (try 120–200)
use_calibri = True                     # set False if Calibri isn't available
out_dir   = r"C:\Users\Public\Fidelity_AGG"  # folder to save the figure/CSV (optional)
save_png  = True
# =======================================================

os.makedirs(out_dir, exist_ok=True)
if use_calibri:
    mpl.rcParams['font.family'] = 'Calibri'

# ------------------------ helpers ------------------------
def _ravel(a):
    """Flatten arrays and lists-of-arrays."""
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()

def extract_amplitude_uV(ds):
    """
    Extract a demodulated amplitude trace from a QCoDeS dataset.
    Priority: 'R' (magnitude) -> sqrt(X^2+Y^2) -> |X|.
    Returns a 1D array in microvolts (µV).
    """
    pdata = ds.get_parameter_data()
    for dep_key, block in pdata.items():
        keys = list(block.keys())
        # direct amplitude?
        amp_keys = [k for k in keys if any(s in k.lower() for s in
                    [" sample r","sample r"," amplitude","magnitude","_r",":r","_ampl","amplitude"])]
        if amp_keys:
            R = _ravel(block[amp_keys[0]]) * 1e6  # V -> µV
            return R[np.isfinite(R)]
        # X/Y -> magnitude
        x_keys = [k for k in keys if any(s in k.lower() for s in [" sample x","sample x","_x",":x"])]
        y_keys = [k for k in keys if any(s in k.lower() for s in [" sample y","sample y","_y",":y"])]
        if x_keys and y_keys:
            X = _ravel(block[x_keys[0]]); Y = _ravel(block[y_keys[0]])
            R = np.sqrt(X**2 + Y**2) * 1e6
            return R[np.isfinite(R)]
        # fallback: X only
        if x_keys:
            X = _ravel(block[x_keys[0]]) * 1e6
            return np.abs(X[np.isfinite(X)])
    raise RuntimeError("No R/X/Y channel found in dataset.")

def make_hist(y, bins=120):
    """Make a non-normalized histogram. Returns (bin_centers, counts, bin_width)."""
    counts, edges = np.histogram(y, bins=bins)
    centers = 0.5*(edges[:-1] + edges[1:])
    bw = np.diff(edges).mean()
    return centers, counts.astype(float), bw

# ---------- Gaussian model (two peaks) and analytics ----------
def g1(x, mu, s, A):
    """Single Gaussian: height A (not area), center mu, std s."""
    return A * np.exp(-0.5 * ((x - mu)/s)**2)

def g2(x, mu0, mu1, s0, A0, s1, A1):
    """Sum of two Gaussians (no background)."""
    return g1(x, mu0, s0, A0) + g1(x, mu1, s1, A1)

def bayes_threshold(mu0, s0, mu1, s1):
    """
    Bayes-optimal threshold for equal priors, solving N(mu0,s0)=N(mu1,s1).
    Returns the root closest to the midpoint.
    """
    a = 1/(s1**2) - 1/(s0**2)
    b = -2*(mu1/(s1**2) - mu0/(s0**2))
    c = (mu1**2)/(s1**2) - (mu0**2)/(s0**2) - 2*np.log(s1/s0)
    if abs(a) < 1e-15:  # nearly equal variances
        return (mu0*s1**2 - mu1*s0**2) / (s1**2 - s0**2)
    roots = np.real(np.roots([a, b, c]))
    xm = 0.5*(mu0 + mu1)
    return float(roots[np.argmin(np.abs(roots - xm))])

def Phi(z):  # standard normal CDF
    return 0.5 * (1.0 + erf(z/np.sqrt(2.0)))

def fidelity_from_params(mu0, s0, mu1, s1, thr):
    """
    Fidelity with equal priors and threshold 'thr':
    F = 0.5 * [ P0(x<thr) + P1(x>thr) ].
    """
    if s0 <= 0 or s1 <= 0:
        return np.nan, np.nan, np.nan
    eps_lo = 1 - Phi((thr - mu0)/s0)   # low→high
    eps_hi =      Phi((thr - mu1)/s1)   # high→low
    F = 1 - 0.5*(eps_lo + eps_hi)
    return F, eps_lo, eps_hi

def snr_from_params(mu0, s0, mu1, s1):
    """SNR commonly used for two-lobe readout: 2*(Δμ)^2 / (σ0^2 + σ1^2)."""
    return 2.0 * (mu1 - mu0)**2 / (s0**2 + s1**2)

def two_means_threshold(v):
    """Fast 1D ISODATA two-means threshold to get robust initial guesses."""
    v = np.asarray(v); v = v[np.isfinite(v)]
    if v.size == 0: return np.nan
    t = np.median(v)
    for _ in range(100):
        lo = v[v <= t]; hi = v[v > t]
        if len(lo) == 0 or len(hi) == 0: break
        t_new = 0.5*(lo.mean() + hi.mean())
        if abs(t_new - t) < 1e-9: t = t_new; break
        t = t_new
    return float(t)

# ==================== AGGREGATE ALL RUNS ====================
initialise_or_create_database_at(db_path)

all_amp = []
for rid in range(run_first, run_last+1):
    try:
        ds = load_by_id(rid)
        R_uV = extract_amplitude_uV(ds)  # µV
        if R_uV.size:
            all_amp.append(R_uV)
            print(f"[run {rid}] collected {R_uV.size} samples")
        else:
            print(f"[run {rid}] empty trace (skipped)")
    except Exception as e:
        print(f"[run {rid}] skipped: {e}")

if not all_amp:
    raise RuntimeError("No valid data found across the requested runs.")

R_all = np.concatenate(all_amp)
R_all = R_all[np.isfinite(R_all)]
print(f"\nAggregated samples: {R_all.size}  |  min={np.min(R_all):.2f} µV, max={np.max(R_all):.2f} µV")

# ==================== HISTOGRAM + TWO-GAUSSIAN FIT ====================
ticks, counts, bw = make_hist(R_all, bins=bins)

# Robust initial guesses from the aggregated data
thr_guess = two_means_threshold(R_all)
lo = R_all[R_all <= thr_guess]; hi = R_all[R_all > thr_guess]
if lo.size < 10 or hi.size < 10:  # fallback
    p30, p70 = np.percentile(R_all, [30, 70])
    lo = R_all[R_all <= p30]; hi = R_all[R_all >= p70]

mu0, s0 = (np.mean(lo), np.std(lo, ddof=1)) if lo.size else (np.min(R_all), np.std(R_all, ddof=1))
mu1, s1 = (np.mean(hi), np.std(hi, ddof=1)) if hi.size else (np.max(R_all), np.std(R_all, ddof=1))
if mu0 > mu1:  # enforce low < high
    mu0, mu1, s0, s1 = mu1, mu0, s1, s0

A0 = counts[np.argmin(np.abs(ticks - mu0))] if len(ticks) else counts.max()
A1 = counts[np.argmin(np.abs(ticks - mu1))] if len(ticks) else counts.max()
p0 = [mu0, mu1, max(s0, 1e-6), max(A0, 1.0), max(s1, 1e-6), max(A1, 1.0)]

# Fit the two-Gaussian model to histogram counts
popt, pcov = curve_fit(g2, ticks, counts, p0=p0, maxfev=30000)
mu0_f, mu1_f, s0_f, A0_f, s1_f, A1_f = popt
if mu0_f > mu1_f:
    mu0_f, mu1_f = mu1_f, mu0_f
    s0_f,  s1_f  = s1_f, s0_f
    A0_f,  A1_f  = A1_f, A0_f

# Thresholds and metrics
thr_mid   = 0.5*(mu0_f + mu1_f)
thr_bayes = bayes_threshold(mu0_f, s0_f, mu1_f, s1_f)
F_mid,  eps_lo_m, eps_hi_m  = fidelity_from_params(mu0_f, s0_f, mu1_f, s1_f, thr_mid)
F_bay,  eps_lo_b, eps_hi_b  = fidelity_from_params(mu0_f, s0_f, mu1_f, s1_f, thr_bayes)
SNR_val = snr_from_params(mu0_f, s0_f, mu1_f, s1_f)

# ==================== PLOT (AGGREGATED) ====================
xplot = np.linspace(ticks.min(), ticks.max(), 2000)
fit_tot = g2(xplot, mu0_f, mu1_f, s0_f, A0_f, s1_f, A1_f)
fit_0   = g1(xplot, mu0_f, s0_f, A0_f)
fit_1   = g1(xplot, mu1_f, s1_f, A1_f)

plt.figure(figsize=(7.2, 4.6))
plt.bar(ticks, counts, width=bw, color="0.75", edgecolor="none", alpha=0.65, label=f"All runs {run_first}-{run_last}")
plt.plot(xplot, fit_tot, "-", lw=2.0, label="Fit (G1+G2)")
plt.plot(xplot, fit_0,  "--", lw=1.2, label="G1")
plt.plot(xplot, fit_1,  "--", lw=1.2, label="G2")
plt.axvline(thr_mid,   ls="--",  lw=1.4, color="tab:red",   label=f"Mid = {thr_mid:.1f} µV")
plt.axvline(thr_bayes, ls="-.",  lw=1.4, color="tab:green", label=f"Bayes = {thr_bayes:.1f} µV")
plt.xlabel("amplitude (µV)", fontsize=13)
plt.ylabel("counts", fontsize=13)
plt.title(f"Aggregated histogram (runs {run_first}–{run_last})", fontsize=13)
plt.legend(fontsize=10)
plt.grid(alpha=0.3)
plt.tight_layout()

if save_png:
    fig_path = os.path.join(out_dir, f"agg_hist_2G_runs_{run_first}_{run_last}.png")
    plt.savefig(fig_path, dpi=220)
    print("Saved figure ->", fig_path)

plt.show()

# ==================== PRINT SUMMARY ====================
print("\n=== Aggregated two-Gaussian fit ===")
print(f"mu_low  = {mu0_f:.2f} µV   sd_low  = {s0_f:.2f} µV")
print(f"mu_high = {mu1_f:.2f} µV   sd_high = {s1_f:.2f} µV")
print(f"Mid threshold  = {thr_mid:.2f} µV   |  Fidelity_mid  = {F_mid*100:.2f}%  "
      f"(err low→high={eps_lo_m*100:.2f}%, high→low={eps_hi_m*100:.2f}%)")
print(f"Bayes threshold= {thr_bayes:.2f} µV |  Fidelity_bayes= {F_bay*100:.2f}% "
      f"(err low→high={eps_lo_b*100:.2f}%, high→low={eps_hi_b*100:.2f}%)")
print(f"SNR = {SNR_val:.2f}")
