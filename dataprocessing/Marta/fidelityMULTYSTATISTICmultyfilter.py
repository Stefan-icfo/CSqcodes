import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict
from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit

# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v34_25_11_25.db"
run_first, run_last = 32, 101    # runs to aggregate
bins      = 160
use_calibri = True
units_in_volts = True                 # True if dataset is in V; False if already in µV
out_dir   = r"C:\Users\Public\Fidelity_AGG_filter_scan"
save_png  = True
# Optional manual initial guess (None = auto)
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


def parse_filter_khz_from_expname(exp_name):
    """
    Parse filter frequency (in kHz) from the experiment/run name,
    e.g. '_fidelity50mVL50kdemodtimetrace' -> 50 kHz.
    """
    if exp_name is None:
        return None
    s = str(exp_name).lower()
    m = re.search(r'(\d+)\s*k', s)   # finds '50k', '100k', etc.
    if not m:
        return None
    return float(m.group(1))


def extract_trace_by_filter_from_ds(ds, units_in_volts=True, debug=False):
    """
    For a single dataset:
      - read filter from ds.exp_name (string with 50k / 100k / ...),
      - take the FIRST dependent parameter in get_parameter_data()
        (typically 'y'),
      - return dict {filter_khz: 1D array in µV}.
    """
    # ---- filter from experiment name ----
    try:
        exp_name = ds.exp_name
    except AttributeError:
        exp_name = getattr(ds, "name", None)

    fk = parse_filter_khz_from_expname(exp_name)
    if debug:
        print(f"\n[DEBUG] run_id={ds.run_id}, exp_name={exp_name!r}, filter={fk} kHz")

    if fk is None:
        raise RuntimeError(f"Could not parse filter from exp_name={exp_name!r}")

    # ---- get data from first dependent parameter ----
    pdata = ds.get_parameter_data()

    if debug:
        print("[DEBUG] get_parameter_data keys:", [repr(k) for k in pdata.keys()])

    # take first dep_key
    dep_keys = list(pdata.keys())
    if not dep_keys:
        raise RuntimeError("No dependent parameters in dataset")

    dep_key = dep_keys[0]
    block = pdata[dep_key]

    if debug:
        print(f"[DEBUG] using dep_key={dep_key!r}, block keys={list(block.keys())}")

    # the actual data array is block[dep_key] if present, otherwise take first value
    if dep_key in block:
        arr = np.asarray(block[dep_key])
    else:
        first_vkey = list(block.keys())[0]
        arr = np.asarray(block[first_vkey])

    if arr.size == 0:
        raise RuntimeError("Empty data array")

    scale = 1e6 if units_in_volts else 1.0
    X = _ravel(arr) * scale
    X = X[np.isfinite(X)]
    if X.size == 0:
        raise RuntimeError("Only NaNs in data array")

    return {fk: X}


def make_hist(y, bins=120):
    """Non-normalized histogram -> (centers, counts, bin_width)."""
    counts, edges = np.histogram(y, bins=bins)
    centers = 0.5*(edges[:-1] + edges[1:])
    bw = np.diff(edges).mean()
    return centers, counts.astype(float), bw


def two_means_threshold(v):
    """1D ISODATA threshold to split two lobes (for initial guesses)."""
    v = np.asarray(v); v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan
    t = np.median(v)
    for _ in range(100):
        lo = v[v <= t]; hi = v[v > t]
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


def fit_fidelity_from_samples(X_all, bins=160, p_guess=None, label=""):
    """
    Do histogram + 2-Gaussian fit and return dict with results.
    """
    ticks, counts, bw = make_hist(X_all, bins=bins)

    # Initial guess
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
        p0 = np.array([x0, x1, max(w1, 1e-6), h1,
                       max(w2, 1e-6), h2], dtype=float)
    else:
        p0 = np.array(p_guess, dtype=float)

    # Fit
    lb = [-np.inf, -np.inf, 1e-12, 0.0, 1e-12, 0.0]
    ub = [ np.inf,  np.inf,  np.inf, np.inf,  np.inf, np.inf]
    popt, pcov = curve_fit(f, ticks, counts, p0=p0,
                           bounds=(lb, ub), maxfev=40000)
    x0, x1, w1, h1, w2, h2 = popt
    if x0 > x1:
        x0, x1 = x1, x0
        w1, w2 = w2, w1
        h1, h2 = h2, h1

    threshold = 0.5*(x0 + x1)
    arr = np.linspace(ticks.min(), ticks.max(), 200_001)
    F0 = np.sum(g1(arr[arr < threshold], x0, w1, h1)) / np.sum(g1(arr, x0, w1, h1))
    F1 = np.sum(g1(arr[arr > threshold], x1, w2, h2)) / np.sum(g1(arr, x1, w2, h2))
    Fidelity = 0.5*(F0 + F1)
    SNR = 2.0 * (x1 - x0)**2 / (w1**2 + w2**2)

    print(f"\n=== Two-Gaussian fit ({label}) ===")
    print(f"x0 = {x0:.2f} µV,  w1 = {w1:.2f} µV,  h1 = {h1:.1f}")
    print(f"x1 = {x1:.2f} µV,  w2 = {w2:.2f} µV,  h2 = {h2:.1f}")
    print(f"Midpoint threshold = {threshold:.2f} µV")
    print(f"F0(left) = {F0*100:.2f}%   |   F1(right) = {F1*100:.2f}%")
    print(f"Fidelity  = {Fidelity*100:.2f}%   |   SNR = {SNR:.2f}")

    # Plot for this filter
    xplot = np.linspace(ticks.min(), ticks.max(), 2000)
    fit_total = f(xplot, x0, x1, w1, h1, w2, h2)
    g_low     = g1(xplot, x0, w1, h1)
    g_high    = g1(xplot, x1, w2, h2)

    plt.figure(figsize=(7.4, 4.8))
    plt.bar(ticks, counts, width=bw, color="0.75", edgecolor="none", alpha=0.65,
            label=f"{label}")
    plt.plot(xplot, fit_total, "-", lw=2.0, label="Fit (G1+G2)")
    plt.plot(xplot, g_low,  "--", lw=1.2, label="G1")
    plt.plot(xplot, g_high, "--", lw=1.2, label="G2")
    plt.axvline(threshold, ls="--", lw=1.4, color="tab:red",
                label=f"mid = {threshold:.1f} µV")
    plt.xlabel("amplitude (µV)", fontsize=13)
    plt.ylabel("counts", fontsize=13)
    plt.title(f"Histogram ({label}) and two-Gaussian fit", fontsize=13)
    plt.legend(fontsize=8)
    plt.tight_layout()

    return {
        "x0": x0, "x1": x1, "w1": w1, "w2": w2,
        "h1": h1, "h2": h2,
        "threshold": threshold,
        "F0": F0*100.0,
        "F1": F1*100.0,
        "Fidelity": Fidelity*100.0,
        "SNR": SNR
    }


# ==================== AGGREGATE BY FILTER ====================
initialise_or_create_database_at(db_path)

all_X_by_filter = defaultdict(list)

for rid in range(run_first, run_last+1):
    try:
        ds = load_by_id(rid)
        debug = (rid in (run_first, run_first+1))  # debug per i primi due run
        filt_dict = extract_trace_by_filter_from_ds(ds, units_in_volts=units_in_volts, debug=debug)
        for fk, arr in filt_dict.items():
            all_X_by_filter[fk].append(arr)
        print(f"[run {rid}] collected samples for filter {list(filt_dict.keys())} kHz")
    except Exception as e:
        print(f"[run {rid}] skipped: {e}")

if not all_X_by_filter:
    raise RuntimeError("No data found across the requested runs.")


# concatena per filtro e info base
for fk in list(all_X_by_filter.keys()):
    X_all = np.concatenate(all_X_by_filter[fk])
    X_all = X_all[np.isfinite(X_all)]
    all_X_by_filter[fk] = X_all
    print(f"Filter {fk:.0f} kHz -> {X_all.size} samples, "
          f"range = [{np.min(X_all):.2f}, {np.max(X_all):.2f}] µV")


# ==================== FIT PER FILTER & SUMMARY TABLE ====================
results = []

for fk in sorted(all_X_by_filter.keys()):
    X_all = all_X_by_filter[fk]
    label = f"runs {run_first}-{run_last}, filter {fk:.0f} kHz"
    res = fit_fidelity_from_samples(X_all, bins=bins, p_guess=p_guess, label=label)
    results.append({"filter_khz": fk, **res})

    if save_png:
        fig_path = os.path.join(
            out_dir,
            f"hist_2G_{int(fk)}kHz_runs_{run_first}_{run_last}.png"
        )
        plt.savefig(fig_path, dpi=220)
        print("Saved figure ->", fig_path)

    plt.show()


# -------- table: filter vs fidelity --------
print("\n=== Fidelity vs filter ===")
print("{:>8}  {:>9}  {:>9}  {:>11}  {:>9}".format(
    "filter", "F0 (%)", "F1 (%)", "Fidelity (%)", "SNR"))
for r in sorted(results, key=lambda x: x["filter_khz"]):
    print("{:8.0f}  {:9.2f}  {:9.2f}  {:11.2f}  {:9.2f}".format(
        r["filter_khz"], r["F0"], r["F1"], r["Fidelity"], r["SNR"]))

# -------- plot: Fidelity vs filter --------
filters = [r["filter_khz"] for r in sorted(results, key=lambda x: x["filter_khz"])]
fids    = [r["Fidelity"]   for r in sorted(results, key=lambda x: x["filter_khz"])]

plt.figure(figsize=(5.5, 3.8))
plt.plot(filters, fids, "o-")
plt.xlabel("filter (kHz)", fontsize=13)
plt.ylabel("Fidelity (%)", fontsize=13)
plt.title(f"Fidelity vs filter (runs {run_first}-{run_last})", fontsize=13)
plt.grid(True)
plt.tight_layout()

if save_png:
    fig_path = os.path.join(out_dir,
        f"fidelity_vs_filter_runs_{run_first}_{run_last}.png")
    plt.savefig(fig_path, dpi=220)
    print("Saved figure ->", fig_path)

plt.show()
