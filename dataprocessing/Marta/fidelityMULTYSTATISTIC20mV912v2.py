import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit
from scipy.stats import norm

# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v47_08_12_25.db"

bins_max = 200                # max bins; actual bins will be adaptive
use_calibri = True
units_in_volts = True         # True if dataset is in V; False if already in µV
out_dir = r"C:\Users\Public\Fidelity_AGG_X_signed_only"
save_png = True
show_group_plots = False      # show histogram+fit for each group?

fidelity_decimals = 18       # printing precision

# ------------------ integration groups (241..860, blocks of 10) ------------------
integration_groups = []

start_run  = 241
end_run    = 860
block_size = 10

t0_ns = 100.0      # 241–250
dt_ns = 150.0      # +150 ns every 10 runs

k = 0
run_first = start_run
while run_first <= end_run:
    run_last = min(run_first + block_size - 1, end_run)
    t_ns = t0_ns + k * dt_ns

    if t_ns >= 1000:
        label = f"{t_ns/1000:.2f} us"
    else:
        label = f"{t_ns:.0f} ns"

    integration_groups.append({
        "label": label,
        "t_ns": float(t_ns),
        "run_first": int(run_first),
        "run_last": int(run_last),
    })

    k += 1
    run_first += block_size

# sanity check expected last block
if integration_groups[-1]["run_first"] == 851 and integration_groups[-1]["run_last"] == 860:
    assert abs(integration_groups[-1]["t_ns"] - 9250.0) < 1e-12, "Last group is not 9.25 us (9250 ns)!"

# Optional manual initial guess for the fit (None = auto)
p_guess = None
# =======================================================

os.makedirs(out_dir, exist_ok=True)
if use_calibri:
    mpl.rcParams['font.family'] = 'Calibri'

# ------------------------ helpers ------------------------
def _ravel(a):
    """Flatten arrays and lists-of-arrays to 1D."""
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()

def extract_X_signed_uV(ds, units_in_volts=True, debug=False):
    """
    Strictly extract a SIGNED 'X' (Real) demod channel.
    - NO fallback to magnitude R.
    - If multiple 'X' candidates exist, return the first.
    Returns a 1D array (µV). Raises if not found.
    """
    scale = 1e6 if units_in_volts else 1.0
    pdata = ds.get_parameter_data()

    for dep_key, block in pdata.items():
        keys = list(block.keys())
        if debug:
            print(f"[dep='{dep_key}'] keys: {keys}")

        for k in keys:
            kl = k.lower().strip()
            if (
                kl == "x" or
                kl.endswith(" sample x") or
                "sample x" in kl or
                kl.endswith(":x") or
                kl.endswith("_x") or
                " real" in kl or
                "(x" in kl or
                kl.startswith("x ")
            ):
                X = _ravel(block[k]) * scale
                X = X[np.isfinite(X)]
                if X.size == 0:
                    continue
                return X

    raise RuntimeError("Signed X channel not found. Print keys and adjust the patterns above.")

def make_hist(y, bins=120):
    """Non-normalized histogram -> (centers, counts, bin_width)."""
    counts, edges = np.histogram(y, bins=bins)
    centers = 0.5*(edges[:-1] + edges[1:])
    bw = np.diff(edges).mean()
    return centers, counts.astype(float), bw

def two_means_threshold(v):
    """1D ISODATA-style threshold to split two lobes (for initial guesses)."""
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan
    t = np.median(v)
    for _ in range(100):
        lo = v[v <= t]
        hi = v[v > t]
        if len(lo) == 0 or len(hi) == 0:
            break
        t_new = 0.5*(lo.mean() + hi.mean())
        if abs(t_new - t) < 1e-9:
            t = t_new
            break
        t = t_new
    return float(t)

# Two-Gaussian model (sum of two unnormalized Gaussians)
def f(x, *p):
    # p = (x0, x1, w1, h1, w2, h2)
    x0, x1, w1, h1, w2, h2 = p
    return (h1*np.exp(-(x-x0)**2/(2*w1**2)) +
            h2*np.exp(-(x-x1)**2/(2*w2**2)))

def g1(x, x0, w, h):
    return h*np.exp(-(x-x0)**2/(2*w**2))

def optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2):
    """
    Threshold t where the two Gaussians are equal:
        h1 * exp(-(t-x0)^2/(2 w1^2)) = h2 * exp(-(t-x1)^2/(2 w2^2))
    (equal priors boundary)
    """
    if h1 <= 0 or h2 <= 0:
        return 0.5 * (x0 + x1)

    A = 1.0 / (w2**2) - 1.0 / (w1**2)
    B = -2.0 * x1 / (w2**2) + 2.0 * x0 / (w1**2)
    C = (x1**2) / (w2**2) - (x0**2) / (w1**2) - 2.0 * np.log(h2 / h1)

    if abs(A) < 1e-18:
        if abs(B) < 1e-18:
            return 0.5 * (x0 + x1)
        return float(-C / B)

    roots = np.roots([A, B, C])
    roots = roots[np.isreal(roots)].real
    if roots.size == 0:
        return 0.5 * (x0 + x1)

    mid = 0.5 * (x0 + x1)
    between = [r for r in roots if min(x0, x1) <= r <= max(x0, x1)]
    if between:
        return float(min(between, key=lambda r: abs(r - mid)))
    return float(min(roots, key=lambda r: abs(r - mid)))

def format_percent_floor(p, decimals=4):
    """Format probability p (0–1) as percent with floor (no rounding up)."""
    p = max(0.0, min(1.0, float(p)))
    factor = 10**decimals
    v = np.floor(p * 100.0 * factor) / factor
    return f"{v:.{decimals}f}"

def format_scientific_prob(p, decimals=6):
    """Format probability in scientific notation."""
    p = max(0.0, float(p))
    if p == 0.0:
        return f"0.{'0'*decimals}e+00"
    return f"{p:.{decimals}e}"

# ==================== ANALYSIS FOR ONE GROUP ====================
def analyze_group(run_first, run_last, group_label, t_ns, show_plot=False):
    print("\n" + "="*70)
    print(f"Integration time: {t_ns:.1f} ns  |  group '{group_label}'")
    print(f"Runs {run_first}–{run_last}")
    print("="*70)

    all_X = []
    for rid in range(run_first, run_last+1):
        try:
            ds = load_by_id(rid)
            debug = (rid in (run_first, run_first+1))
            X = extract_X_signed_uV(ds, units_in_volts=units_in_volts, debug=debug)
            if X.size:
                all_X.append(X)
                # print(f"[run {rid}] collected {X.size} samples (SIGNED X)")
            else:
                print(f"[run {rid}] empty trace (skipped)")
        except Exception as e:
            print(f"[run {rid}] skipped: {e}")

    if not all_X:
        raise RuntimeError(f"No X-signed data found for runs {run_first}-{run_last}.")

    X_all = np.concatenate(all_X)
    X_all = X_all[np.isfinite(X_all)]

    print(f"Aggregated samples: {X_all.size}  |  range = [{np.min(X_all):.2f}, {np.max(X_all):.2f}] µV")

    if np.min(X_all) >= 0:
        print("WARNING: no negative values found -> likely NOT signed X. Check extract patterns.")

    # adaptive bins (stabilizes fit at long integration times)
    bins_eff = min(bins_max, max(50, int(np.sqrt(X_all.size))))
    ticks, counts, bw = make_hist(X_all, bins=bins_eff)

    # Initial guess
    if p_guess is None:
        thr = two_means_threshold(X_all)
        lo = X_all[X_all <= thr]
        hi = X_all[X_all > thr]
        if lo.size < 10 or hi.size < 10:
            p30, p70 = np.percentile(X_all, [30, 70])
            lo = X_all[X_all <= p30]
            hi = X_all[X_all >= p70]

        x0 = np.mean(lo) if lo.size else np.min(X_all)
        x1 = np.mean(hi) if hi.size else np.max(X_all)
        if x0 > x1:
            x0, x1 = x1, x0

        w1 = np.std(lo, ddof=1) if lo.size > 1 else np.std(X_all, ddof=1)
        w2 = np.std(hi, ddof=1) if hi.size > 1 else np.std(X_all, ddof=1)

        def peak_height(xc):
            i = np.argmin(np.abs(ticks - xc))
            return max(counts[i], 1.0)

        h1 = peak_height(x0)
        h2 = peak_height(x1)
        p0 = np.array([x0, x1, max(w1, 1e-6), h1, max(w2, 1e-6), h2], dtype=float)
    else:
        p0 = np.array(p_guess, dtype=float)

    lb = [-np.inf, -np.inf, 1e-12, 0.0, 1e-12, 0.0]
    ub = [ np.inf,  np.inf,  np.inf, np.inf,  np.inf, np.inf]

    # Poisson weighting for histogram counts (stabilizes the fit)
    sigma = np.sqrt(np.maximum(counts, 1.0))

    popt, pcov = curve_fit(
        f, ticks, counts,
        p0=p0, bounds=(lb, ub),
        sigma=sigma, absolute_sigma=True,
        maxfev=60000
    )

    x0, x1, w1, h1, w2, h2 = popt

    # enforce x0 < x1
    if x0 > x1:
        x0, x1 = x1, x0
        w1, w2 = w2, w1
        h1, h2 = h2, h1

    threshold = optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2)

    # ---- Fidelity (analytic, stable) ----
    z0 = (threshold - x0) / w1
    z1 = (threshold - x1) / w2
    F0 = norm.cdf(z0)
    F1 = 1.0 - norm.cdf(z1)
    Fidelity = 0.5 * (F0 + F1)

    # SNR (your definition)
    SNR = 2.0 * (x1 - x0)**2 / (w1**2 + w2**2)

    print(f"x0={x0:.2f} µV, w1={w1:.2f} µV | x1={x1:.2f} µV, w2={w2:.2f} µV")
    print(f"thr={threshold:.2f} µV | F0={format_percent_floor(F0, fidelity_decimals)}% "
          f"| F1={format_percent_floor(F1, fidelity_decimals)}%")
    print(f"Fidelity={format_percent_floor(Fidelity, fidelity_decimals)}% | SNR={SNR:.2f}")

    # --------------------------- HISTOGRAM PLOT ---------------------------
    xplot = np.linspace(ticks.min(), ticks.max(), 4000)
    fit_total = f(xplot, x0, x1, w1, h1, w2, h2)
    g_low = g1(xplot, x0, w1, h1)
    g_high = g1(xplot, x1, w2, h2)

    fig = plt.figure(figsize=(7.4, 4.8))
    plt.bar(ticks, counts, width=bw, color="0.75", edgecolor="none", alpha=0.65,
            label=f"Runs {run_first}-{run_last} (N={X_all.size})")
    plt.yscale('log')
    plt.ylim(bottom=1)
    plt.plot(xplot, fit_total, "-", lw=2.0, label="Fit (G1+G2)")
    plt.plot(xplot, g_low, "--", lw=1.2, label="G1")
    plt.plot(xplot, g_high, "--", lw=1.2, label="G2")
    plt.axvline(threshold, ls="--", lw=1.4, color="tab:red",
                label=f"thr = {threshold:.1f} µV")
    plt.xlabel("Amplitude X (µV)", fontsize=13)
    plt.ylabel("Counts", fontsize=13)
    plt.title(f"Histogram & two-Gaussian fit ({group_label}, {t_ns:.0f} ns)", fontsize=13)
    plt.legend(fontsize=8)
    plt.tight_layout()

    if save_png:
        fig_path = os.path.join(
            out_dir,
            f"agg_hist_2G_Xsigned_{group_label.replace(' ','_')}_{int(t_ns)}ns_runs_{run_first}_{run_last}.png"
        )
        fig.savefig(fig_path, dpi=220)
        print("Saved figure ->", fig_path)

    if show_plot:
        plt.show()
    else:
        plt.close(fig)

    return Fidelity, SNR

# ==================== MAIN ====================
initialise_or_create_database_at(db_path)

all_times_ns = []
all_fidelities = []
all_snrs = []

for grp in integration_groups:
    Fid, snr = analyze_group(
        run_first=grp["run_first"],
        run_last=grp["run_last"],
        group_label=grp.get("label", f"{grp['t_ns']} ns"),
        t_ns=grp["t_ns"],
        show_plot=show_group_plots,
    )
    all_times_ns.append(grp["t_ns"])
    all_fidelities.append(Fid)
    all_snrs.append(snr)

all_times_ns = np.array(all_times_ns, dtype=float)
all_fidelities = np.array(all_fidelities, dtype=float)
all_snrs = np.array(all_snrs, dtype=float)

# Sort by integration time
order = np.argsort(all_times_ns)
all_times_ns = all_times_ns[order]
all_fidelities = all_fidelities[order]
all_snrs = all_snrs[order]

# ==================== PLOT: (1 - Fidelity) vs integration time ====================
errors = 1.0 - all_fidelities

# Avoid zeros for log plot (numerical floor only for plotting)
errors_plot = np.clip(errors, 1e-20, None)

plt.figure(figsize=(6.0, 4.2))
plt.plot(all_times_ns / 1000.0, errors_plot, "o-", lw=1.8, ms=5)  # µs on x-axis
plt.xlabel("Integration time (µs)", fontsize=13)
plt.ylabel("1 - Fidelity", fontsize=13)
plt.yscale("log")
plt.title("Readout error vs integration time", fontsize=13)
plt.grid(True, which="both", ls="--", alpha=0.4)
plt.tight_layout()

if save_png:
    fig_path = os.path.join(out_dir, "error_vs_integration_time.png")
    plt.savefig(fig_path, dpi=220)
    print("Saved figure ->", fig_path)

plt.show()

# ==================== TEXT SUMMARY ====================
print("\n=== SUMMARY BY INTEGRATION TIME ===")
for t, F, snr in zip(all_times_ns, all_fidelities, all_snrs):
    err = 1.0 - F
    err_str = format_scientific_prob(err, decimals=fidelity_decimals)
    t_us = t / 1000.0
    print(
        f"t = {t_us:7.3f} us | "
        f"Fidelity = {format_percent_floor(F, decimals=fidelity_decimals)}% "
        f"| 1 - F = {err_str} | SNR = {snr:.2f}"
    )
