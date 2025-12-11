import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit
from scipy.special import log_ndtr  # stable log(CDF)
# =======================================================
# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v47_08_12_25.db"

out_dir = r"C:\Users\Public\Fidelity_AGG_X_signed_only"
save_png = True
show_group_plots = False

use_calibri = True
units_in_volts = True  # True if dataset in V; False if already in µV

bins_max = 500
fidelity_decimals = 18

# Runs scheme
start_run  = 241
end_run    = 860
block_size = 10
t0_ns      = 100.0   # first block 241–250
dt_ns      = 150.0   # increment per block

# If run 241 corresponds to prepared state "0", set True. If it is state "1", set False.
first_run_is_state0 = True

# Numerical/plot floors (plot only)
plot_floor_abs = 1e-18         # floor shown in log plot to avoid zeros
stat_floor_factor = 0.5        # ~0.5/N statistical floor (order of magnitude)

# Optional manual fit guess: (x0, x1, w1, h1, w2, h2) in µV, counts
p_guess = None
# =======================================================

os.makedirs(out_dir, exist_ok=True)
if use_calibri:
    mpl.rcParams["font.family"] = "Calibri"

# ------------------ integration groups ------------------
integration_groups = []
k = 0
run_first = start_run
while run_first <= end_run:
    run_last = min(run_first + block_size - 1, end_run)
    t_ns = t0_ns + k * dt_ns
    label = f"{t_ns/1000:.2f} us" if t_ns >= 1000 else f"{t_ns:.0f} ns"
    integration_groups.append(
        {"label": label, "t_ns": float(t_ns), "run_first": int(run_first), "run_last": int(run_last)}
    )
    k += 1
    run_first += block_size

# sanity check for last group (851–860 should be 9250 ns)
if integration_groups and integration_groups[-1]["run_first"] == 851 and integration_groups[-1]["run_last"] == 860:
    assert abs(integration_groups[-1]["t_ns"] - 9250.0) < 1e-9, "Last group is not 9.25 us (9250 ns)!"

# ------------------------ helpers ------------------------
def _ravel(a):
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()

def extract_X_signed_uV(ds, units_in_volts=True, debug=False):
    """
    Strictly extract a SIGNED X (Real) demod channel.
    Returns 1D array in µV.
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
                kl == "x"
                or kl.endswith(" sample x")
                or "sample x" in kl
                or kl.endswith(":x")
                or kl.endswith("_x")
                or " real" in kl
                or "(x" in kl
                or kl.startswith("x ")
            ):
                X = _ravel(block[k]) * scale
                X = X[np.isfinite(X)]
                if X.size:
                    return X

    raise RuntimeError("Signed X channel not found. Print keys and adjust extract_X_signed_uV patterns.")

def make_hist(y, bins=120):
    counts, edges = np.histogram(y, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bw = np.diff(edges).mean()
    return centers, counts.astype(float), bw

def two_means_threshold(v):
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
        t_new = 0.5 * (lo.mean() + hi.mean())
        if abs(t_new - t) < 1e-9:
            t = t_new
            break
        t = t_new
    return float(t)

# Two-Gaussian model (unnormalized heights for histogram counts)
def f2g(x, x0, x1, w1, h1, w2, h2):
    return (h1 * np.exp(-(x - x0) ** 2 / (2 * w1 ** 2)) +
            h2 * np.exp(-(x - x1) ** 2 / (2 * w2 ** 2)))

def optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2):
    """
    Solve h1*exp(-(t-x0)^2/(2w1^2)) = h2*exp(-(t-x1)^2/(2w2^2))
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

def run_is_state0(run_id: int) -> bool:
    """
    Alternating labeling anchored at start_run:
    if first_run_is_state0: start_run->0, start_run+1->1, ...
    else: start_run->1, start_run+1->0, ...
    """
    offset = run_id - start_run
    if first_run_is_state0:
        return (offset % 2 == 0)
    return (offset % 2 == 1)

def format_percent_floor(p, decimals=6):
    p = max(0.0, min(1.0, float(p)))
    factor = 10**decimals
    v = np.floor(p * 100.0 * factor) / factor
    return f"{v:.{decimals}f}"

def format_scientific_prob(p, decimals=6):
    p = float(p)
    if p <= 0:
        return f"0.{'0'*decimals}e+00"
    return f"{p:.{decimals}e}"

# ==================== ANALYSIS FOR ONE GROUP ====================
def analyze_group(run_first, run_last, group_label, t_ns, show_plot=False):
    print("\n" + "=" * 78)
    print(f"Integration time: {t_ns/1000.0:.3f} us | group '{group_label}' | runs {run_first}-{run_last}")
    print("=" * 78)

    X0_list, X1_list = [], []

    for rid in range(run_first, run_last + 1):
        try:
            ds = load_by_id(rid)
            debug = (rid in (run_first, run_first + 1))
            X = extract_X_signed_uV(ds, units_in_volts=units_in_volts, debug=debug)
            if X.size == 0:
                continue
            if run_is_state0(rid):
                X0_list.append(X)
            else:
                X1_list.append(X)
        except Exception as e:
            print(f"[run {rid}] skipped: {e}")

    if not X0_list or not X1_list:
        raise RuntimeError(
            f"Missing state0 or state1 data in runs {run_first}-{run_last}. "
            f"(state0? {bool(X0_list)} | state1? {bool(X1_list)})"
        )

    X0 = np.concatenate(X0_list)
    X1 = np.concatenate(X1_list)
    X0 = X0[np.isfinite(X0)]
    X1 = X1[np.isfinite(X1)]

    N0, N1 = int(X0.size), int(X1.size)
    X_all = np.concatenate([X0, X1])
    X_all = X_all[np.isfinite(X_all)]
    N = int(X_all.size)

    print(f"N0={N0}, N1={N1}, Ntot={N} | range=[{X_all.min():.2f},{X_all.max():.2f}] µV")

    # histogram for fit
    bins_eff = min(bins_max, max(80, int(np.sqrt(N))))
    ticks, counts, bw = make_hist(X_all, bins=bins_eff)

    # initial guess
    if p_guess is None:
        thr0 = two_means_threshold(X_all)
        lo = X_all[X_all <= thr0]
        hi = X_all[X_all > thr0]
        if lo.size < 10 or hi.size < 10:
            p30, p70 = np.percentile(X_all, [30, 70])
            lo = X_all[X_all <= p30]
            hi = X_all[X_all >= p70]

        x0 = float(np.mean(lo) if lo.size else np.min(X_all))
        x1 = float(np.mean(hi) if hi.size else np.max(X_all))
        if x0 > x1:
            x0, x1 = x1, x0

        w1 = float(np.std(lo, ddof=1) if lo.size > 1 else np.std(X_all, ddof=1))
        w2 = float(np.std(hi, ddof=1) if hi.size > 1 else np.std(X_all, ddof=1))

        def peak_height(xc):
            i = int(np.argmin(np.abs(ticks - xc)))
            return float(max(counts[i], 1.0))

        h1 = peak_height(x0)
        h2 = peak_height(x1)
        p0 = np.array([x0, x1, max(w1, 1e-6), h1, max(w2, 1e-6), h2], dtype=float)
    else:
        p0 = np.array(p_guess, dtype=float)

    lb = [-np.inf, -np.inf, 1e-12, 0.0, 1e-12, 0.0]
    ub = [ np.inf,  np.inf,  np.inf, np.inf,  np.inf, np.inf]

    # Poisson weighting on histogram counts
    sigma = np.sqrt(np.maximum(counts, 1.0))

    popt, _ = curve_fit(
        f2g, ticks, counts,
        p0=p0, bounds=(lb, ub),
        sigma=sigma, absolute_sigma=True,
        maxfev=120000
    )

    x0, x1, w1, h1, w2, h2 = [float(v) for v in popt]

    # enforce x0<x1
    if x0 > x1:
        x0, x1 = x1, x0
        w1, w2 = w2, w1
        h1, h2 = h2, h1

    thr = optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2)

    # SNR (your definition)
    SNR = 2.0 * (x1 - x0) ** 2 / (w1**2 + w2**2)

    # =================== MODEL ERROR (stable, log-domain) ===================
    # decision rule: decide 0 if X < thr else 1
    z0 = (thr - x0) / w1
    z1 = (thr - x1) / w2

    # err0_model = P0(X > thr) = 1 - Phi(z0) = Phi(-z0)
    # err1_model = P1(X < thr) = Phi(z1)
    log_err0_model = log_ndtr(-z0)
    log_err1_model = log_ndtr(z1)

    # err_model = 0.5*(err0+err1) in log domain
    m = max(log_err0_model, log_err1_model)
    log_sum = m + np.log(np.exp(log_err0_model - m) + np.exp(log_err1_model - m))
    log_err_model = np.log(0.5) + log_sum

    err0_model = float(np.exp(log_err0_model))
    err1_model = float(np.exp(log_err1_model))
    err_model  = float(np.exp(log_err_model))
    F_model = 1.0 - err_model  # may print as 1.0, but err_model remains accurate

    # =================== EMPIRICAL ERROR (using labels) ===================
    err0_emp = float(np.mean(X0 > thr))   # prepared 0 but threshold says 1
    err1_emp = float(np.mean(X1 < thr))   # prepared 1 but threshold says 0
    err_emp  = 0.5 * (err0_emp + err1_emp)
    F_emp    = 1.0 - err_emp

    print(f"fit: x0={x0:.2f} w1={w1:.2f} | x1={x1:.2f} w2={w2:.2f} | thr={thr:.2f} µV | SNR={SNR:.2f}")
    print(f"MODEL: err0={err0_model:.3e} err1={err1_model:.3e} err={err_model:.3e}")
    print(f"EMP  : err0={err0_emp:.3e} err1={err1_emp:.3e} err={err_emp:.3e}")

    # ---------------- plot per-group ----------------
    if show_plot or save_png:
        xplot = np.linspace(ticks.min(), ticks.max(), 4000)
        fit_total = f2g(xplot, x0, x1, w1, h1, w2, h2)

        fig = plt.figure(figsize=(7.4, 4.8))
        plt.bar(ticks, counts, width=bw, color="0.75", edgecolor="none", alpha=0.65, label=f"N={N}")
        plt.yscale("log")
        plt.ylim(bottom=1)
        plt.plot(xplot, fit_total, "-", lw=2.0, label="Fit (G1+G2)")
        plt.axvline(thr, ls="--", lw=1.4, color="tab:red", label=f"thr={thr:.1f} µV")
        plt.xlabel("Amplitude X (µV)")
        plt.ylabel("Counts")
        plt.title(f"{group_label}  ({t_ns/1000:.3f} µs)")
        plt.legend(fontsize=9)
        plt.tight_layout()

        if save_png:
            fig_path = os.path.join(
                out_dir,
                f"hist_2G_{group_label.replace(' ','_')}_{int(t_ns)}ns_runs_{run_first}_{run_last}.png"
            )
            fig.savefig(fig_path, dpi=220)

        if show_plot:
            plt.show()
        else:
            plt.close(fig)

    return {
        "t_ns": float(t_ns),
        "F_model": float(F_model),
        "err_model": float(err_model),
        "F_emp": float(F_emp),
        "err_emp": float(err_emp),
        "err0_emp": float(err0_emp),
        "err1_emp": float(err1_emp),
        "SNR": float(SNR),
        "N": int(N),
    }

# ==================== MAIN ====================
initialise_or_create_database_at(db_path)

rows = []
for grp in integration_groups:
    rows.append(
        analyze_group(
            run_first=grp["run_first"],
            run_last=grp["run_last"],
            group_label=grp["label"],
            t_ns=grp["t_ns"],
            show_plot=show_group_plots,
        )
    )

# sort by time
rows = sorted(rows, key=lambda d: d["t_ns"])

t_ns = np.array([r["t_ns"] for r in rows], float)
t_us = t_ns / 1000.0

err_model = np.array([r["err_model"] for r in rows], float)
err_emp   = np.array([r["err_emp"] for r in rows], float)
N         = np.array([r["N"] for r in rows], float)

# plot floors (only for plot)
stat_floor = stat_floor_factor / np.maximum(N, 1.0)
floor = np.maximum(stat_floor, plot_floor_abs)

err_model_plot = np.maximum(err_model, floor)
err_emp_plot   = np.maximum(err_emp, floor)

plt.figure(figsize=(6.8, 4.8))
plt.plot(t_us, err_emp_plot, "o-", lw=1.8, ms=4, label="Empirical (from alternating 0/1)")
plt.plot(t_us, err_model_plot, "o--", lw=1.4, ms=3.5, label="Gaussian-tail model (stable)")
plt.yscale("log")
plt.xlabel("Integration time (µs)", fontsize=13)
plt.ylabel("1 - Fidelity", fontsize=13)
plt.title("Readout error vs integration time", fontsize=13)
plt.grid(True, which="both", ls="--", alpha=0.4)
plt.legend(fontsize=9)
plt.tight_layout()

if save_png:
    fig_path = os.path.join(out_dir, "error_vs_integration_time.png")
    plt.savefig(fig_path, dpi=220)
    print("Saved figure ->", fig_path)

plt.show()

print("\n=== SUMMARY BY INTEGRATION TIME ===")
for r in rows:
    t_us_i = r["t_ns"] / 1000.0
    Fm = r["F_model"]
    Fe = r["F_emp"]
    em = r["err_model"]
    ee = r["err_emp"]
    print(
        f"t = {t_us_i:7.3f} us | "
        f"F_model = {format_percent_floor(Fm, fidelity_decimals)}% | 1-F_model = {format_scientific_prob(em, fidelity_decimals)} | "
        f"F_emp = {format_percent_floor(Fe, fidelity_decimals)}% | 1-F_emp = {format_scientific_prob(ee, fidelity_decimals)} | "
        f"err0={r['err0_emp']:.3e}, err1={r['err1_emp']:.3e} | SNR={r['SNR']:7.2f} | N={r['N']}"
    )
