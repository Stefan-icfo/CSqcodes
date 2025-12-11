import os  
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit

# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v38_28_11_25.db"
bins      = 160
use_calibri = True
units_in_volts = True                 # True if dataset is in V; False if already in µV
out_dir   = r"C:\Users\Public\Fidelity_AGG_X_signed_only"
save_png  = True

# Number of decimal places for Fidelity in percent
# Example: 6 -> 99.999993%
fidelity_decimals = 15

# List of integration-time groups.
# Each group corresponds to a different integration time with its run range.
integration_groups = [
    {
        "label": "29.9 ns",
        "t_ns": 29.9,
        "run_first": 740,
        "run_last": 746,
    },
    {
        "label": "49.7 ns",
        "t_ns": 49.7,
        "run_first": 747,
        "run_last": 753,
    },
    {
        "label": "60.7 ns",
        "t_ns": 60.7,
        "run_first": 754,
        "run_last": 760,
    },
    {
        "label": "102.6 ns",
        "t_ns": 102.6,
        "run_first": 761,
        "run_last": 767,
    },
    {
        "label": "130 ns",
        "t_ns": 130.0,
        "run_first": 768,
        "run_last": 774,
    },
    {
        "label": "150 ns",
        "t_ns": 150.0,
        "run_first": 775,
        "run_last": 781,
    },
]


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

        # candidate X key patterns (very permissive)
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
    """
    1D ISODATA-style threshold to split two lobes (for initial guesses).
    """
    v = np.asarray(v); v = v[np.isfinite(v)]
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

# Two-Gaussian model
def f(x, *p):
    # p = (x0, x1, w1, h1, w2, h2)
    x0, x1, w1, h1, w2, h2 = p
    return (h1*np.exp(-(x-x0)**2/(2*w1**2)) +
            h2*np.exp(-(x-x1)**2/(2*w2**2)))

def g1(x, x0, w, h):
    return h*np.exp(-(x-x0)**2/(2*w**2))

# ---------- Optimal threshold for unequal widths ----------
def optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2):
    """
    Return the discrimination threshold where the two Gaussians are equal:
        h1 * exp(-(t-x0)^2/(2 w1^2)) = h2 * exp(-(t-x1)^2/(2 w2^2))
    This is the Bayes-optimal boundary for equal priors.
    """
    if h1 <= 0 or h2 <= 0:
        return 0.5 * (x0 + x1)

    A = 1.0 / (w2**2) - 1.0 / (w1**2)
    B = -2.0 * x1 / (w2**2) + 2.0 * x0 / (w1**2)
    C = (x1**2) / (w2**2) - (x0**2) / (w1**2) - 2.0 * np.log(h2 / h1)

    if abs(A) < 1e-18:
        # Widths effectively equal -> linear equation
        if abs(B) < 1e-18:
            return 0.5 * (x0 + x1)
        t = -C / B
        return float(t)

    roots = np.roots([A, B, C])
    # Keep only real roots
    roots = roots[np.isreal(roots)].real
    if roots.size == 0:
        return 0.5 * (x0 + x1)

    mid = 0.5 * (x0 + x1)
    between = [r for r in roots if min(x0, x1) <= r <= max(x0, x1)]
    if between:
        return float(min(between, key=lambda r: abs(r - mid)))
    return float(min(roots, key=lambda r: abs(r - mid)))

# ---------- Percent formatter without rounding up to 100 ----------
def format_percent_floor(p, decimals=4):
    """
    Format a probability p (0–1) as a percentage with `decimals` decimals,
    using floor instead of round so that ~100% never prints as >100.0000%.
    Example: p = 0.999999 -> '99.9999' (for decimals=4).
    """
    p = max(0.0, min(1.0, float(p)))  # clamp to [0,1]
    factor = 10**decimals
    v = np.floor(p * 100.0 * factor) / factor
    return f"{v:.{decimals}f}"

# ---------- Scientific formatter for small probabilities ----------
def format_scientific_prob(p, decimals=6):
    """
    Format a small probability p (0–1) in scientific notation with a given
    number of decimals in the mantissa. Example:
        decimals = 6 -> '1.234567e-06'
    This gives you a precision consistent with the Fidelity formatting.
    """
    p = max(0.0, float(p))
    if p == 0.0:
        return f"0.{'0'*decimals}e+00"
    return f"{p:.{decimals}e}"

# ==================== ANALYSIS FOR ONE GROUP ====================
def analyze_group(run_first, run_last, group_label, t_ns):
    """
    Analyze runs [run_first, run_last] for a single integration time.
    Returns (Fidelity_float, SNR).
    """
    print("\n" + "="*70)
    print(f"Integration time: {t_ns} ns  |  group '{group_label}'")
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
                print(f"[run {rid}] collected {X.size} samples (SIGNED X)")
            else:
                print(f"[run {rid}] empty trace (skipped)")
        except Exception as e:
            print(f"[run {rid}] skipped: {e}")

    if not all_X:
        raise RuntimeError(f"No X-signed data found for runs {run_first}-{run_last}.")

    X_all = np.concatenate(all_X)
    X_all = X_all[np.isfinite(X_all)]
    print(f"\nAggregated samples: {X_all.size}  |  range = [{np.min(X_all):.2f}, {np.max(X_all):.2f}] µV")

    if np.min(X_all) >= 0:
        print("WARNING: no negative values found. You are likely not reading a signed X channel.")
        print("Check the printed keys and adjust the matching rules inside 'extract_X_signed_uV'.")

    # -------- HISTOGRAM + TWO-GAUSSIAN FIT --------
    ticks, counts, bw = make_hist(X_all, bins=bins)

    # Initial guess (auto) unless p_guess provided
    if p_guess is None:
        thr = two_means_threshold(X_all)
        lo = X_all[X_all <= thr]; hi = X_all[X_all > thr]
        if lo.size < 10 or hi.size < 10:
            p30, p70 = np.percentile(X_all, [30, 70])
            lo = X_all[X_all <= p30]; hi = X_all[X_all >= p70]
        x0 = np.mean(lo) if lo.size else np.min(X_all)
        x1 = np.mean(hi) if hi.size else np.max(X_all)
        if x0 > x1:
            x0, x1 = x1, x0
        w1 = np.std(lo, ddof=1) if lo.size > 1 else np.std(X_all, ddof=1)
        w2 = np.std(hi, ddof=1) if hi.size > 1 else np.std(X_all, ddof=1)

        def peak_height(xc):
            if len(ticks) == 0:
                return 1.0
            i = np.argmin(np.abs(ticks - xc))
            return max(counts[i], 1.0)

        h1 = peak_height(x0); h2 = peak_height(x1)
        p0 = np.array([x0, x1, max(w1, 1e-6), h1, max(w2, 1e-6), h2], dtype=float)
    else:
        p0 = np.array(p_guess, dtype=float)

    lb = [-np.inf, -np.inf, 1e-12, 0.0, 1e-12, 0.0]
    ub = [ np.inf,  np.inf,  np.inf, np.inf,  np.inf, np.inf]
    popt, pcov = curve_fit(f, ticks, counts, p0=p0, bounds=(lb, ub), maxfev=40000)
    x0, x1, w1, h1, w2, h2 = popt
    if x0 > x1:
        x0, x1 = x1, x0
        w1, w2 = w2, w1
        h1, h2 = h2, h1

    threshold = optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2)

    arr = np.linspace(ticks.min(), ticks.max(), 200_001)
    F0 = np.sum(g1(arr[arr < threshold], x0, w1, h1)) / np.sum(g1(arr, x0, w1, h1))
    F1 = np.sum(g1(arr[arr > threshold], x1, w2, h2)) / np.sum(g1(arr, x1, w2, h2))
    Fidelity = 0.5*(F0 + F1)
    SNR = 2.0 * (x1 - x0)**2 / (w1**2 + w2**2)

    print("\n=== Aggregated two-Gaussian fit (SIGNED X ONLY) ===")
    print(f"x0 = {x0:.2f} µV,  w1 = {w1:.2f} µV,  h1 = {h1:.1f}")
    print(f"x1 = {x1:.2f} µV,  w2 = {w2:.2f} µV,  h2 = {h2:.1f}")
    print(f"Optimal threshold = {threshold:.2f} µV")

    # Use floor-style formatting for Fidelity (percent) with user-defined decimals
    print(f"F0(left)  = {format_percent_floor(F0, decimals=fidelity_decimals)}%")
    print(f"F1(right) = {format_percent_floor(F1, decimals=fidelity_decimals)}%")
    print(f"Fidelity  = {format_percent_floor(Fidelity, decimals=fidelity_decimals)}%   |   SNR = {SNR:.2f}")

    # --------------------------- HISTOGRAM PLOT ---------------------------
    xplot = np.linspace(ticks.min(), ticks.max(), 20000)
    fit_total = f(xplot, x0, x1, w1, h1, w2, h2)
    g_low     = g1(xplot, x0, w1, h1)
    g_high    = g1(xplot, x1, w2, h2)

    plt.figure(figsize=(7.4, 4.8))
    plt.bar(ticks, counts, width=bw, color="0.75", edgecolor="none", alpha=0.65,
            label=f"Aggregated runs {run_first}-{run_last}")
    plt.plot(xplot, fit_total, "-", lw=2.0, label="Fit (G1+G2)")
    plt.plot(xplot, g_low,  "--", lw=1.2, label="G1")
    plt.plot(xplot, g_high, "--", lw=1.2, label="G2")
    plt.axvline(threshold, ls="--", lw=1.4, color="tab:red",
                label=f"opt. thr. = {threshold:.1f} µV")
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
        plt.savefig(fig_path, dpi=220)
        print("Saved figure ->", fig_path)

    plt.show()

    return Fidelity, SNR

# ==================== LOOP OVER ALL INTEGRATION TIMES ====================
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
    )
    all_times_ns.append(grp["t_ns"])
    all_fidelities.append(Fid)
    all_snrs.append(snr)

all_times_ns = np.array(all_times_ns, dtype=float)
all_fidelities = np.array(all_fidelities, dtype=float)
all_snrs = np.array(all_snrs, dtype=float)

# Sort by increasing integration time
order = np.argsort(all_times_ns)
all_times_ns = all_times_ns[order]
all_fidelities = all_fidelities[order]
all_snrs = all_snrs[order]

# ==================== PLOT 1 - (1 - Fidelity) vs integration time ====================
errors = 1.0 - all_fidelities   # 1 - Fidelity (error probability)

plt.figure(figsize=(6.0, 4.2))
plt.plot(all_times_ns, errors, "o-", lw=1.8, ms=5)
plt.xlabel("Integration time (ns)", fontsize=13)
plt.ylabel("1 - Fidelity", fontsize=13)
# If you want log scale for the error, uncomment the next line:
# plt.yscale("log")
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
    err_str = format_scientific_prob(err, decimals=fidelity_decimals)  # same precision as Fidelity
    print(
        f"t = {t:7.1f} ns | "
        f"Fidelity = {format_percent_floor(F, decimals=fidelity_decimals)}% "
        f"| 1 - F = {err_str} | SNR = {snr:.2f}"
    )

