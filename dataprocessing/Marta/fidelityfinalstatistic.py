import os, re, csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict
from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit

# ==================== USER SETTINGS ====================
db_path        = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v28_10_11_25.db"
run_first      = 400
run_last       = 670
units_in_volts = True       # True: convert V->µV with *1e6; False if already µV
bins           = 160        # histogram bins for the aggregated fit
use_calibri    = True
out_dir        = r"C:\Users\Public\Fidelity_vs_Filter"
save_png       = True
save_csv       = True
# =======================================================

os.makedirs(out_dir, exist_ok=True)
if use_calibri:
    mpl.rcParams['font.family'] = 'Calibri'

# ------------------------ helpers ------------------------
def _ravel(a):
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()

def extract_X_signed_uV(ds, units_in_volts=True, debug=False):
    """
    STRICTLY return a SIGNED 'X/Real' demod trace (µV).
    No fallback to magnitude R (so negatives are kept).
    If multiple candidates, the first match is returned.
    """
    scale = 1e6 if units_in_volts else 1.0
    pdata = ds.get_parameter_data()
    for dep_key, block in pdata.items():
        keys = list(block.keys())
        if debug:
            print(f"[dep='{dep_key}'] keys:", keys)
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
                if X.size:
                    return X
    raise RuntimeError("Signed X channel not found—adjust key patterns above (print keys with debug=True).")

def get_filter_k_from_ds(ds):
    """
    Parse the filter value 'NNk' from dataset/experiment name.
    Matches things like 'fidelity70k', 'Fidelity 120K', etc.
    Returns an integer (e.g., 70 for '70k') or None if not found.
    """
    names = []
    # Try common attributes
    for attr in ("name", "exp_name"):
        v = getattr(ds, attr, None)
        if isinstance(v, str):
            names.append(v)
    # Try Experiment object name
    try:
        en = getattr(getattr(ds, "experiment", None), "name", None)
        if isinstance(en, str):
            names.append(en)
    except Exception:
        pass
    # Combine and search
    joined = " | ".join(names)
    m = re.search(r"fidelity\s*([0-9]+)\s*[kK]", joined)
    if not m:
        m = re.search(r"([0-9]+)\s*[kK]", joined)  # fallback: any “NNk”
    if m:
        return int(m.group(1))
    return None

def make_hist(y, bins=120):
    counts, edges = np.histogram(y, bins=bins)
    centers = 0.5*(edges[:-1] + edges[1:])
    bw = np.diff(edges).mean()
    return centers, counts.astype(float), bw

def two_means_threshold(v):
    """ISODATA two-means threshold for robust initial guesses."""
    v = np.asarray(v); v = v[np.isfinite(v)]
    if v.size == 0: return np.nan
    t = np.median(v)
    for _ in range(100):
        lo = v[v <= t]; hi = v[v > t]
        if len(lo)==0 or len(hi)==0: break
        t_new = 0.5*(lo.mean() + hi.mean())
        if abs(t_new - t) < 1e-9: t = t_new; break
        t = t_new
    return float(t)

# ------ colleague-style Gaussian model + numeric fidelity ------
def f2g(x, x0, x1, w1, h1, w2, h2):
    return h1*np.exp(-(x-x0)**2/(2*w1**2)) + h2*np.exp(-(x-x1)**2/(2*w2**2))

def g1(x, x0, w, h):
    return h*np.exp(-(x-x0)**2/(2*w**2))

def fit_two_gaussians_and_metrics(samples, bins=160, label=""):
    """
    Aggregate histogram -> fit two Gaussians -> metrics:
    - threshold: midpoint (as in your colleague's snippet)
    - Fidelity: 0.5*(F0 + F1) via numerical sums
    - SNR: 2*(Δμ)^2/(σ1^2+σ2^2)
    Returns dict with fit params and metrics.
    """
    ticks, counts, bw = make_hist(samples, bins=bins)

    # robust init
    thr = two_means_threshold(samples)
    lo = samples[samples <= thr]; hi = samples[samples > thr]
    if lo.size < 10 or hi.size < 10:
        p30, p70 = np.percentile(samples, [30, 70])
        lo = samples[samples <= p30]; hi = samples[samples >= p70]
    x0 = np.mean(lo) if lo.size else np.min(samples)
    x1 = np.mean(hi) if hi.size else np.max(samples)
    if x0 > x1: x0, x1 = x1, x0
    w1 = np.std(lo, ddof=1) if lo.size > 1 else np.std(samples, ddof=1)
    w2 = np.std(hi, ddof=1) if hi.size > 1 else np.std(samples, ddof=1)
    def peak_height(xc):
        if len(ticks) == 0: return 1.0
        i = np.argmin(np.abs(ticks - xc))
        return max(counts[i], 1.0)
    h1 = peak_height(x0); h2 = peak_height(x1)
    p0 = np.array([x0, x1, max(w1, 1e-6), h1, max(w2, 1e-6), h2], dtype=float)

    # fit with positivity bounds on widths/heights
    lb = [-np.inf, -np.inf, 1e-12, 0.0, 1e-12, 0.0]
    ub = [ np.inf,  np.inf,  np.inf, np.inf,  np.inf, np.inf]
    popt, pcov = curve_fit(f2g, ticks, counts, p0=p0, bounds=(lb, ub), maxfev=40000)
    x0, x1, w1, h1, w2, h2 = popt
    if x0 > x1:   # enforce low<high
        x0, x1 = x1, x0;  w1, w2 = w2, w1;  h1, h2 = h2, h1

    # midpoint threshold & numerical fractions (colleague style)
    threshold = 0.5*(x0 + x1)
    grid = np.linspace(ticks.min(), ticks.max(), 200_001)
    F0 = np.sum(g1(grid[grid < threshold], x0, w1, h1)) / np.sum(g1(grid, x0, w1, h1))
    F1 = np.sum(g1(grid[grid > threshold], x1, w2, h2)) / np.sum(g1(grid, x1, w2, h2))
    fidelity = 0.5*(F0 + F1)
    snr = 2.0 * (x1 - x0)**2 / (w1**2 + w2**2)

    return dict(
        ticks=ticks, counts=counts, bw=bw,
        x0=x0, x1=x1, w1=w1, h1=h1, w2=w2, h2=h2,
        threshold=threshold, F0=F0, F1=F1, fidelity=fidelity, snr=snr
    )

# ==================== SCAN RUNS & GROUP BY FILTER ====================
initialise_or_create_database_at(db_path)

samples_by_filter = defaultdict(list)   # key: filter_k (int), value: list of arrays
missing_filter = []

for rid in range(run_first, run_last+1):
    try:
        ds = load_by_id(rid)
        # grab filter value from name(s)
        filt_k = get_filter_k_from_ds(ds)
        if filt_k is None:
            missing_filter.append(rid)
            continue
        # signed X samples
        debug = (rid in (run_first, run_first+1))
        X = extract_X_signed_uV(ds, units_in_volts=units_in_volts, debug=debug)
        if X.size:
            samples_by_filter[filt_k].append(X)
            print(f"[run {rid}] filter={filt_k}k | {X.size} samples (SIGNED X)")
    except Exception as e:
        print(f"[run {rid}] skipped: {e}")

if missing_filter:
    print("WARNING: filter value not found for runs:", missing_filter)

if not samples_by_filter:
    raise RuntimeError("No data grouped by filter—check name parsing or ranges.")

# ==================== METRICS PER FILTER ====================
rows = []
for filt_k in sorted(samples_by_filter.keys()):
    samples = np.concatenate(samples_by_filter[filt_k])
    if samples.size == 0:
        continue
    if np.min(samples) >= 0:
        print(f"NOTE: filter {filt_k}k has no negative values—check you are reading signed X.")
    res = fit_two_gaussians_and_metrics(samples, bins=bins, label=f"{filt_k}k")
    rows.append(dict(
        filter_k=filt_k,
        fidelity=res["fidelity"],
        F0=res["F0"], F1=res["F1"],
        snr=res["snr"],
        x0=res["x0"], x1=res["x1"], w1=res["w1"], w2=res["w2"],
        threshold=res["threshold"]
    ))
    print(f"filter {filt_k:>3}k  |  Fidelity={res['fidelity']*100:6.2f}%  |  SNR={res['snr']:.2f}  "
          f"|  mid={res['threshold']:.2f} µV")

# Sort by filter value
rows = sorted(rows, key=lambda d: d["filter_k"])
filters = np.array([r["filter_k"] for r in rows], dtype=float)
fids    = np.array([r["fidelity"] for r in rows], dtype=float)
snrs    = np.array([r["snr"] for r in rows], dtype=float)

# ==================== PLOTS ====================
plt.figure(figsize=(7.4, 4.6))
plt.plot(filters, fids*100, "o-", lw=1.5, ms=5)
plt.xlabel("filter (k)", fontsize=13)
plt.ylabel("Fidelity (%)", fontsize=13)
plt.title(f"Fidelity vs filter (runs {run_first}–{run_last})", fontsize=13)
plt.grid(alpha=0.3)
plt.tight_layout()
if save_png:
    plt.savefig(os.path.join(out_dir, f"fidelity_vs_filter_{run_first}_{run_last}.png"), dpi=220)
plt.show()

plt.figure(figsize=(7.4, 4.6))
plt.plot(filters, snrs, "o-", lw=1.5, ms=5)
plt.xlabel("filter (k)", fontsize=13)
plt.ylabel("SNR (2·Δμ²/(σ₁²+σ₂²))", fontsize=13)
plt.title(f"SNR vs filter (runs {run_first}–{run_last})", fontsize=13)
plt.grid(alpha=0.3)
plt.tight_layout()
if save_png:
    plt.savefig(os.path.join(out_dir, f"SNR_vs_filter_{run_first}_{run_last}.png"), dpi=220)
plt.show()

# ==================== CSV SUMMARY ====================
if save_csv:
    csv_path = os.path.join(out_dir, f"fidelity_snr_vs_filter_{run_first}_{run_last}.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["filter_k", "fidelity", "F0", "F1", "snr", "x0", "x1", "w1", "w2", "threshold"])
        for r in rows:
            w.writerow([r["filter_k"], r["fidelity"], r["F0"], r["F1"], r["snr"],
                        r["x0"], r["x1"], r["w1"], r["w2"], r["threshold"]])
    print("Saved CSV ->", csv_path)




