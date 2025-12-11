import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit
from typing import Optional

# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v38_28_11_25.db"
bins      = 160
use_calibri = True
units_in_volts = True
out_dir   = r"C:\Users\Public\Fidelity_AGG_X_signed_only"
save_png  = True
fidelity_decimals = 15

# Scan range of run_ids to auto-build groups
RUN_SCAN_FIRST = 740
RUN_SCAN_LAST  = 781

# Optional manual initial guess for the fit (None = auto)
p_guess = None
# =======================================================

os.makedirs(out_dir, exist_ok=True)
if use_calibri:
    mpl.rcParams['font.family'] = 'Calibri'


# ==================== AUTO GROUPING (from run name) ====================
# matches "...C3.350us..." or "...TC3.350us..."
TC_RE = re.compile(r"(?:\bTC|\bC)\s*([0-9]+(?:\.[0-9]+)?)\s*us", re.IGNORECASE)

def _get_run_label_string(ds) -> str:
    # Try a few common attributes
    for attr in ("exp_name", "name", "run_name"):
        if hasattr(ds, attr):
            v = getattr(ds, attr)
            try:
                v = v() if callable(v) else v
            except Exception:
                pass
            if isinstance(v, str) and v:
                return v

    # Fallback: run description
    try:
        desc = ds.get_run_description()
        if isinstance(desc, dict):
            for k in ("name", "experiment_name", "exp_name"):
                if k in desc and isinstance(desc[k], str) and desc[k]:
                    return desc[k]
            return str(desc)
        return str(desc)
    except Exception:
        return ""

def tc_us_from_dataset(ds) -> Optional[float]:
    s = _get_run_label_string(ds)
    m = TC_RE.search(s)
    if not m:
        return None
    return float(m.group(1))  # microseconds

def build_integration_groups_from_runs(run_first: int, run_last: int):
    """
    Build groups for consecutive runs that share the same TC extracted from name.
    Returns list of dicts: {label, t_ns, run_first, run_last}
    """
    groups = []
    current_tc = None
    current_start = None
    current_end = None

    for rid in range(run_first, run_last + 1):
        try:
            ds = load_by_id(rid)
        except Exception:
            continue

        tc_us = tc_us_from_dataset(ds)
        if tc_us is None:
            # close any open group if we hit a run without TC
            if current_tc is not None:
                groups.append({
                    "label": f"{current_tc:.3f} us",
                    "t_ns": current_tc * 1000.0,
                    "run_first": current_start,
                    "run_last": current_end,
                })
                current_tc = current_start = current_end = None
            continue

        if current_tc is None:
            current_tc = tc_us
            current_start = rid
            current_end = rid
        else:
            # extend group only if same tc and consecutive run_id
            if abs(tc_us - current_tc) < 1e-9 and rid == current_end + 1:
                current_end = rid
            else:
                groups.append({
                    "label": f"{current_tc:.3f} us",
                    "t_ns": current_tc * 1000.0,
                    "run_first": current_start,
                    "run_last": current_end,
                })
                current_tc = tc_us
                current_start = rid
                current_end = rid

    if current_tc is not None:
        groups.append({
            "label": f"{current_tc:.3f} us",
            "t_ns": current_tc * 1000.0,
            "run_first": current_start,
            "run_last": current_end,
        })

    return groups


# ==================== HELPERS: DATA + FIT ====================
def _ravel(a):
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()

def extract_X_signed_uV(ds, units_in_volts=True, debug=False):
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

    raise RuntimeError("Signed X channel not found. Print keys and adjust patterns in extract_X_signed_uV().")

def make_hist(y, bins=120):
    counts, edges = np.histogram(y, bins=bins)
    centers = 0.5*(edges[:-1] + edges[1:])
    bw = np.diff(edges).mean()
    return centers, counts.astype(float), bw

def two_means_threshold(v):
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

def f(x, *p):
    x0, x1, w1, h1, w2, h2 = p
    return (h1*np.exp(-(x-x0)**2/(2*w1**2)) +
            h2*np.exp(-(x-x1)**2/(2*w2**2)))

def g1(x, x0, w, h):
    return h*np.exp(-(x-x0)**2/(2*w**2))

def optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2):
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
    p = max(0.0, min(1.0, float(p)))
    factor = 10**decimals
    v = np.floor(p * 100.0 * factor) / factor
    return f"{v:.{decimals}f}"

def format_scientific_prob(p, decimals=6):
    p = max(0.0, float(p))
    if p == 0.0:
        return f"0.{'0'*decimals}e+00"
    return f"{p:.{decimals}e}"


# ==================== ANALYSIS FOR ONE GROUP ====================
def analyze_group(run_first, run_last, group_label, t_ns):
    print("\n" + "="*70)
    print(f"Integration time: {t_ns:.3f} ns  |  group '{group_label}'")
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

    ticks, counts, bw = make_hist(X_all, bins=bins)

    # initial guess
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
            i = np.argmin(np.abs(ticks - xc))
            return max(counts[i], 1.0)

        h1 = peak_height(x0); h2 = peak_height(x1)
        p0 = np.array([x0, x1, max(w1, 1e-6), h1, max(w2, 1e-6), h2], dtype=float)
    else:
        p0 = np.array(p_guess, dtype=float)

    lb = [-np.inf, -np.inf, 1e-12, 0.0, 1e-12, 0.0]
    ub = [ np.inf,  np.inf,  np.inf, np.inf,  np.inf, np.inf]
    popt, _ = curve_fit(f, ticks, counts, p0=p0, bounds=(lb, ub), maxfev=40000)

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
    print(f"Fidelity  = {format_percent_floor(Fidelity, decimals=fidelity_decimals)}%   |   SNR = {SNR:.2f}")

    # Plot histogram (linear scale as in your latest version)
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
    plt.title(f"Histogram & two-Gaussian fit ({group_label})", fontsize=13)
    plt.legend(fontsize=8)
    plt.tight_layout()

    if save_png:
        fig_path = os.path.join(
            out_dir,
            f"agg_hist_2G_Xsigned_{group_label.replace(' ','_')}_runs_{run_first}_{run_last}.png"
        )
        plt.savefig(fig_path, dpi=220)
        print("Saved figure ->", fig_path)

    plt.show()

    return Fidelity, SNR


# ==================== MAIN ====================
initialise_or_create_database_at(db_path)

# auto-build groups
integration_groups = build_integration_groups_from_runs(RUN_SCAN_FIRST, RUN_SCAN_LAST)

print("\nAUTO integration_groups found:")
for g in integration_groups:
    print(g)

if not integration_groups:
    raise RuntimeError("No groups found. Check RUN_SCAN_FIRST/LAST and the name pattern (C...us / TC...us).")

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

order = np.argsort(all_times_ns)
all_times_ns = all_times_ns[order]
all_fidelities = all_fidelities[order]
all_snrs = all_snrs[order]

# error vs integration time
errors = 1.0 - all_fidelities

plt.figure(figsize=(6.0, 4.2))
plt.plot(all_times_ns, errors, "o-", lw=1.8, ms=5)
plt.xlabel("Integration time (ns)", fontsize=13)
plt.ylabel("1 - Fidelity", fontsize=13)
plt.title("Readout error vs integration time", fontsize=13)
plt.grid(True, which="both", ls="--", alpha=0.4)
plt.tight_layout()

if save_png:
    fig_path = os.path.join(out_dir, "error_vs_integration_time.png")
    plt.savefig(fig_path, dpi=220)
    print("Saved figure ->", fig_path)

plt.show()

print("\n=== SUMMARY BY INTEGRATION TIME ===")
for t, F, snr in zip(all_times_ns, all_fidelities, all_snrs):
    err = 1.0 - F
    err_str = format_scientific_prob(err, decimals=fidelity_decimals)
    print(
        f"t = {t:9.1f} ns | "
        f"Fidelity = {format_percent_floor(F, decimals=fidelity_decimals)}% "
        f"| 1 - F = {err_str} | SNR = {snr:.2f}"
    )
