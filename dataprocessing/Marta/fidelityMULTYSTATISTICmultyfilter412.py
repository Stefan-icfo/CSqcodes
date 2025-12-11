import os  
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit

# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v43_04_12_25.db"
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
    {"label": "30 ns",     "t_ns": 30.0,    "run_first": 13,  "run_last": 20},
    {"label": "50 ns",     "t_ns": 50.0,    "run_first": 21,  "run_last": 27},
    {"label": "60 ns",     "t_ns": 60.0,    "run_first": 27,  "run_last": 34},
    {"label": "102 ns",    "t_ns": 102.0,   "run_first": 35,  "run_last": 41},
    {"label": "250 ns",    "t_ns": 250.0,   "run_first": 42,  "run_last": 48},
    {"label": "400 ns",    "t_ns": 400.0,   "run_first": 49,  "run_last": 55},
    {"label": "550 ns",    "t_ns": 550.0,   "run_first": 56,  "run_last": 62},
    {"label": "700 ns",    "t_ns": 700.0,   "run_first": 72,  "run_last": 78},
    {"label": "850 ns",    "t_ns": 850.0,   "run_first": 79,  "run_last": 85},

    {"label": "1.0 us",    "t_ns": 1000.0,  "run_first": 96,  "run_last": 102},
    {"label": "1.3 us",    "t_ns": 1300.0,  "run_first": 103, "run_last": 109},
    {"label": "1.45 us",   "t_ns": 1450.0,  "run_first": 110, "run_last": 116},
    {"label": "1.6 us",    "t_ns": 1600.0,  "run_first": 117, "run_last": 123},
    {"label": "1.75 us",   "t_ns": 1750.0,  "run_first": 129, "run_last": 135},
    {"label": "1.9 us",    "t_ns": 1900.0,  "run_first": 136, "run_last": 142},
    {"label": "2.05 us",   "t_ns": 2050.0,  "run_first": 143, "run_last": 149},
    {"label": "2.2 us",    "t_ns": 2200.0,  "run_first": 151, "run_last": 157},
    {"label": "2.35 us",   "t_ns": 2350.0,  "run_first": 158, "run_last": 164},
    {"label": "2.5 us",    "t_ns": 2500.0,  "run_first": 165, "run_last": 171},
    {"label": "2.65 us",   "t_ns": 2650.0,  "run_first": 172, "run_last": 178},
    {"label": "2.8 us",    "t_ns": 2800.0,  "run_first": 188, "run_last": 195},
    {"label": "2.5 us",    "t_ns": 2500.0,  "run_first": 196, "run_last": 203},

    # -------------------------
    # PARTE "SISTEMATICA"
    # ogni 15 run aumento di 150 ns, coerente col tuo esempio 264–278
    # -------------------------
    {"label": "3.1 us",    "t_ns": 3100.0,  "run_first": 249, "run_last": 263},
    {"label": "3.25 us",   "t_ns": 3250.0,  "run_first": 264, "run_last": 278},
    {"label": "3.4 us",    "t_ns": 3400.0,  "run_first": 279, "run_last": 293},
    {"label": "3.55 us",   "t_ns": 3550.0,  "run_first": 294, "run_last": 308},
    {"label": "3.7 us",    "t_ns": 3700.0,  "run_first": 309, "run_last": 323},
    {"label": "3.85 us",   "t_ns": 3850.0,  "run_first": 324, "run_last": 338},
    {"label": "4.0 us",    "t_ns": 4000.0,  "run_first": 339, "run_last": 353},
    {"label": "4.15 us",   "t_ns": 4150.0,  "run_first": 354, "run_last": 368},
    {"label": "4.3 us",    "t_ns": 4300.0,  "run_first": 369, "run_last": 383},
    {"label": "4.45 us",   "t_ns": 4450.0,  "run_first": 384, "run_last": 398},
    {"label": "4.6 us",    "t_ns": 4600.0,  "run_first": 399, "run_last": 413},
    {"label": "4.75 us",   "t_ns": 4750.0,  "run_first": 414, "run_last": 428},
    {"label": "4.9 us",    "t_ns": 4900.0,  "run_first": 429, "run_last": 443},
    {"label": "5.05 us",   "t_ns": 5050.0,  "run_first": 444, "run_last": 458},
    {"label": "5.2 us",    "t_ns": 5200.0,  "run_first": 459, "run_last": 473},
    {"label": "5.35 us",   "t_ns": 5350.0,  "run_first": 474, "run_last": 488},
    {"label": "5.5 us",    "t_ns": 5500.0,  "run_first": 489, "run_last": 503},
    {"label": "5.65 us",   "t_ns": 5650.0,  "run_first": 504, "run_last": 518},
    {"label": "5.8 us",    "t_ns": 5800.0,  "run_first": 519, "run_last": 533},
    {"label": "5.95 us",   "t_ns": 5950.0,  "run_first": 534, "run_last": 548},
    {"label": "6.1 us",    "t_ns": 6100.0,  "run_first": 549, "run_last": 563},
    {"label": "6.25 us",   "t_ns": 6250.0,  "run_first": 564, "run_last": 578},
    {"label": "6.4 us",    "t_ns": 6400.0,  "run_first": 579, "run_last": 593},
    {"label": "6.55 us",   "t_ns": 6550.0,  "run_first": 594, "run_last": 608},
    {"label": "6.7 us",    "t_ns": 6700.0,  "run_first": 609, "run_last": 623},
    {"label": "6.85 us",   "t_ns": 6850.0,  "run_first": 624, "run_last": 634},
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

    arr = np.linspace(ticks.min(), ticks.max(), 2000_001) #here you can increse the tick to have a better resolution
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
    xplot = np.linspace(ticks.min(), ticks.max(), 2000)
    fit_total = f(xplot, x0, x1, w1, h1, w2, h2)
    g_low     = g1(xplot, x0, w1, h1)
    g_high    = g1(xplot, x1, w2, h2)

    plt.figure(figsize=(7.4, 4.8))
    plt.bar(ticks, counts, width=bw, color="0.75", edgecolor="none", alpha=0.65,
            label=f"Aggregated runs {run_first}-{run_last}")
    plt.yscale('log')
    plt.ylim(bottom=1)
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

