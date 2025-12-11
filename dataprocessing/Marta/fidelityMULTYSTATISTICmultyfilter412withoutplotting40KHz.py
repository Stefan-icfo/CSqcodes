import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit

# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v44_05_12_25.db"
bins = 200
use_calibri = True
units_in_volts = True                 # True if dataset is in V; False if already in µV
out_dir = r"C:\Users\Public\Fidelity_AGG_X_signed_only"
save_png = True

# Show histogram+fit plot for each integration time?
# If False: no windows appear; figures are saved (if save_png=True) and then closed.
show_group_plots = False

# Number of decimal places for Fidelity in percent
# Example: 6 -> 99.999993%
fidelity_decimals = 18

# List of integration-time groups.
# ------------------ integration groups (UPDATED) ------------------

integration_groups = []

# 102 ns block
integration_groups.append({"label": "102 ns", "t_ns": 102.0, "run_first": 39, "run_last": 50})

# Systematic blocks: from run 51 to run 677
# each block is 15 runs, and t_ns increases by 150 ns each block
start_run = 51
end_run = 677
block_size = 15

t0_ns = 250.0       # for runs 51–65
dt_ns = 150.0       # +150 ns per block

k = 0
run_first = start_run
while run_first <= end_run:
    run_last = min(run_first + block_size - 1, end_run)
    t_ns = t0_ns + k * dt_ns

    # label in ns or us (nice formatting)
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

# Extra short-time blocks at the end
integration_groups = [
    # 102 ns
    {"label": "102 ns", "t_ns": 102.0, "run_first": 39, "run_last": 50},

    # Sistematico: da 51 a 677, blocchi da 15 run, +150 ns ogni blocco
    {"label": "250 ns",  "t_ns": 250.0,  "run_first": 51,  "run_last": 65},
    {"label": "400 ns",  "t_ns": 400.0,  "run_first": 66,  "run_last": 80},
    {"label": "550 ns",  "t_ns": 550.0,  "run_first": 81,  "run_last": 95},
    {"label": "700 ns",  "t_ns": 700.0,  "run_first": 96,  "run_last": 110},
    {"label": "850 ns",  "t_ns": 850.0,  "run_first": 111, "run_last": 125},

    {"label": "1.00 us", "t_ns": 1000.0, "run_first": 126, "run_last": 140},
    {"label": "1.15 us", "t_ns": 1150.0, "run_first": 141, "run_last": 155},
    {"label": "1.30 us", "t_ns": 1300.0, "run_first": 156, "run_last": 170},
    {"label": "1.45 us", "t_ns": 1450.0, "run_first": 171, "run_last": 185},
    {"label": "1.60 us", "t_ns": 1600.0, "run_first": 186, "run_last": 200},
    {"label": "1.75 us", "t_ns": 1750.0, "run_first": 201, "run_last": 215},
    {"label": "1.90 us", "t_ns": 1900.0, "run_first": 216, "run_last": 230},
    {"label": "2.05 us", "t_ns": 2050.0, "run_first": 231, "run_last": 245},
    {"label": "2.20 us", "t_ns": 2200.0, "run_first": 246, "run_last": 260},
    {"label": "2.35 us", "t_ns": 2350.0, "run_first": 261, "run_last": 275},
    {"label": "2.50 us", "t_ns": 2500.0, "run_first": 276, "run_last": 290},
    {"label": "2.65 us", "t_ns": 2650.0, "run_first": 291, "run_last": 305},
    {"label": "2.80 us", "t_ns": 2800.0, "run_first": 306, "run_last": 320},
    {"label": "2.95 us", "t_ns": 2950.0, "run_first": 321, "run_last": 335},
    {"label": "3.10 us", "t_ns": 3100.0, "run_first": 336, "run_last": 350},
    {"label": "3.25 us", "t_ns": 3250.0, "run_first": 351, "run_last": 365},
    {"label": "3.40 us", "t_ns": 3400.0, "run_first": 366, "run_last": 380},
    {"label": "3.55 us", "t_ns": 3550.0, "run_first": 381, "run_last": 395},
    {"label": "3.70 us", "t_ns": 3700.0, "run_first": 396, "run_last": 410},
    {"label": "3.85 us", "t_ns": 3850.0, "run_first": 411, "run_last": 425},
    {"label": "4.00 us", "t_ns": 4000.0, "run_first": 426, "run_last": 440},
    {"label": "4.15 us", "t_ns": 4150.0, "run_first": 441, "run_last": 455},
    {"label": "4.30 us", "t_ns": 4300.0, "run_first": 456, "run_last": 470},
    {"label": "4.45 us", "t_ns": 4450.0, "run_first": 471, "run_last": 485},
    {"label": "4.60 us", "t_ns": 4600.0, "run_first": 486, "run_last": 500},
    {"label": "4.75 us", "t_ns": 4750.0, "run_first": 501, "run_last": 515},
    {"label": "4.90 us", "t_ns": 4900.0, "run_first": 516, "run_last": 530},
    {"label": "5.05 us", "t_ns": 5050.0, "run_first": 531, "run_last": 545},
    {"label": "5.20 us", "t_ns": 5200.0, "run_first": 546, "run_last": 560},
    {"label": "5.35 us", "t_ns": 5350.0, "run_first": 561, "run_last": 575},
    {"label": "5.50 us", "t_ns": 5500.0, "run_first": 576, "run_last": 590},
    {"label": "5.65 us", "t_ns": 5650.0, "run_first": 591, "run_last": 605},
    {"label": "5.80 us", "t_ns": 5800.0, "run_first": 606, "run_last": 620},
    {"label": "5.95 us", "t_ns": 5950.0, "run_first": 621, "run_last": 635},
    {"label": "6.10 us", "t_ns": 6100.0, "run_first": 636, "run_last": 650},
    {"label": "6.25 us", "t_ns": 6250.0, "run_first": 651, "run_last": 665},
    {"label": "6.40 us", "t_ns": 6400.0, "run_first": 666, "run_last": 677},  # ultimo bloccetto tronco

    # Blocchi finali
    {"label": "29.9 ns", "t_ns": 29.9, "run_first": 849, "run_last": 855},
    {"label": "49 ns",   "t_ns": 49.0, "run_first": 856, "run_last": 862},
    {"label": "60 ns",   "t_ns": 60.0, "run_first": 863, "run_last": 869},
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

# Two-Gaussian model
def f(x, *p):
    # p = (x0, x1, w1, h1, w2, h2)
    x0, x1, w1, h1, w2, h2 = p
    return (h1*np.exp(-(x-x0)**2/(2*w1**2)) +
            h2*np.exp(-(x-x1)**2/(2*w2**2)))

def g1(x, x0, w, h):
    return h*np.exp(-(x-x0)**2/(2*w**2))

def optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2):
    """
    Return threshold t where the two Gaussians are equal:
        h1 * exp(-(t-x0)^2/(2 w1^2)) = h2 * exp(-(t-x1)^2/(2 w2^2))
    (Bayes boundary for equal priors)
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
    """Format probability p (0–1) as percent with floor (no rounding up to >100%)."""
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
            if len(ticks) == 0:
                return 1.0
            i = np.argmin(np.abs(ticks - xc))
            return max(counts[i], 1.0)

        h1 = peak_height(x0)
        h2 = peak_height(x1)
        p0 = np.array([x0, x1, max(w1, 1e-6), h1, max(w2, 1e-6), h2], dtype=float)
    else:
        p0 = np.array(p_guess, dtype=float)

    lb = [-np.inf, -np.inf, 1e-12, 0.0, 1e-12, 0.0]
    ub = [ np.inf,  np.inf,  np.inf, np.inf,  np.inf, np.inf]

    popt, pcov = curve_fit(f, ticks, counts, p0=p0, bounds=(lb, ub), maxfev=40000)
    x0, x1, w1, h1, w2, h2 = popt

    # enforce x0 < x1
    if x0 > x1:
        x0, x1 = x1, x0
        w1, w2 = w2, w1
        h1, h2 = h2, h1

    threshold = optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2)

    # Fidelity and SNR
    arr = np.linspace(ticks.min(), ticks.max(), 2_000_001)  # high-res integration grid
    F0 = np.sum(g1(arr[arr < threshold], x0, w1, h1)) / np.sum(g1(arr, x0, w1, h1))
    F1 = np.sum(g1(arr[arr > threshold], x1, w2, h2)) / np.sum(g1(arr, x1, w2, h2))
    Fidelity = 0.5*(F0 + F1)
    SNR = 2.0 * (x1 - x0)**2 / (w1**2 + w2**2)

    print("\n=== Aggregated two-Gaussian fit (SIGNED X ONLY) ===")
    print(f"x0 = {x0:.2f} µV,  w1 = {w1:.2f} µV,  h1 = {h1:.1f}")
    print(f"x1 = {x1:.2f} µV,  w2 = {w2:.2f} µV,  h2 = {h2:.1f}")
    print(f"Optimal threshold = {threshold:.2f} µV")
    print(f"F0(left)  = {format_percent_floor(F0, decimals=fidelity_decimals)}%")
    print(f"F1(right) = {format_percent_floor(F1, decimals=fidelity_decimals)}%")
    print(f"Fidelity  = {format_percent_floor(Fidelity, decimals=fidelity_decimals)}%   |   SNR = {SNR:.2f}")

    # --------------------------- HISTOGRAM PLOT ---------------------------
    xplot = np.linspace(ticks.min(), ticks.max(), 2000)
    fit_total = f(xplot, x0, x1, w1, h1, w2, h2)
    g_low = g1(xplot, x0, w1, h1)
    g_high = g1(xplot, x1, w2, h2)

    fig = plt.figure(figsize=(7.4, 4.8))
    plt.bar(
        ticks, counts, width=bw, color="0.75", edgecolor="none", alpha=0.65,
        label=f"Aggregated runs {run_first}-{run_last}"
    )
    plt.yscale('log')
    plt.ylim(bottom=1)
    plt.plot(xplot, fit_total, "-", lw=2.0, label="Fit (G1+G2)")
    plt.plot(xplot, g_low, "--", lw=1.2, label="G1")
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

# Sort by increasing integration time
order = np.argsort(all_times_ns)
all_times_ns = all_times_ns[order]
all_fidelities = all_fidelities[order]
all_snrs = all_snrs[order]

# ==================== PLOT: (1 - Fidelity) vs integration time ====================
errors = 1.0 - all_fidelities

plt.figure(figsize=(6.0, 4.2))
plt.plot(all_times_ns, errors, "o-", lw=1.8, ms=5)
plt.xlabel("Integration time (ns)", fontsize=13)
plt.ylabel("1 - Fidelity", fontsize=13)
# plt.yscale("log")  # uncomment if you want log scale
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
    print(
        f"t = {t:7.1f} ns | "
        f"Fidelity = {format_percent_floor(F, decimals=fidelity_decimals)}% "
        f"| 1 - F = {err_str} | SNR = {snr:.2f}"
    )


