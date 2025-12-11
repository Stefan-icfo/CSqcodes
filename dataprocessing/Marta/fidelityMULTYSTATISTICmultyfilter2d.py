import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit
from scipy.stats import multivariate_normal
from sklearn.mixture import GaussianMixture


# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v38_28_11_25.db"

bins = 160
use_calibri = True
units_in_volts = True        # True if dataset stores V; convert to µV internally
out_dir = r"C:\Users\Public\Fidelity_AGG_XY_2D_pretty_Bayes"
save_png = True

# formatting
fidelity_decimals = 15

# 1D numerical integration grid
integration_grid_points = 2_000_001

# 2D density grid size (higher = smoother, slower)
density_gridsize = 420

# 2D plot limits: use percentiles to avoid crazy outliers
xy_percentiles = (0.2, 99.8)

# Colormap for 2D density (pick the one you like)
density_cmap = "viridis"   # e.g. "viridis", "magma", "inferno", "plasma"

# Integration-time groups
integration_groups = [
    {"label": "29.9 ns", "t_ns": 29.9, "run_first": 528, "run_last": 534},
    {"label": "49.9 ns", "t_ns": 49.9, "run_first": 535, "run_last": 541},
    {"label": "102.6 ns", "t_ns": 102.6, "run_first": 542, "run_last": 548},
    {"label": "120 ns", "t_ns": 120.0, "run_first": 549, "run_last": 555},
    {"label": "130 ns", "t_ns": 130.0, "run_first": 562, "run_last": 570},
    {"label": "150 ns", "t_ns": 150.0, "run_first": 571, "run_last": 582},
    {"label": "200 ns", "t_ns": 200.0, "run_first": 584, "run_last": 590},
    {"label": "250 ns", "t_ns": 250.0, "run_first": 591, "run_last": 597},
    {"label": "300 ns", "t_ns": 300.0, "run_first": 598, "run_last": 604},
    {"label": "350 ns", "t_ns": 350.0, "run_first": 605, "run_last": 611},
    {"label": "400 ns", "t_ns": 400.0, "run_first": 612, "run_last": 618},
    {"label": "500 ns", "t_ns": 500.0, "run_first": 619, "run_last": 625},
    {"label": "600 ns", "t_ns": 600.0, "run_first": 626, "run_last": 632},
    {"label": "700 ns", "t_ns": 700.0, "run_first": 633, "run_last": 639},
    {"label": "800 ns", "t_ns": 800.0, "run_first": 648, "run_last": 654},
    {"label": "1.0 us", "t_ns": 1000.0, "run_first": 655, "run_last": 661},
    {"label": "1.2 us", "t_ns": 1200.0, "run_first": 665, "run_last": 672},
    {"label": "1.4 us", "t_ns": 1400.0, "run_first": 803, "run_last": 809},
    {"label": "1.6 us", "t_ns": 1600.0, "run_first": 679, "run_last": 684},
    {"label": "1.8 us", "t_ns": 1800.0, "run_first": 685, "run_last": 691},
    {"label": "2.0 us", "t_ns": 2000.0, "run_first": 810, "run_last": 816},
    {"label": "2.2 us", "t_ns": 2200.0, "run_first": 817, "run_last": 823},
    {"label": "2.5 us", "t_ns": 2500.0, "run_first": 698, "run_last": 704},
    {"label": "2.7 us", "t_ns": 2700.0, "run_first": 824, "run_last": 830},
    {"label": "3.0 us", "t_ns": 3000.0, "run_first": 705, "run_last": 711},
    {"label": "3.3 us", "t_ns": 3300.0, "run_first": 831, "run_last": 837},
    {"label": "3.5 us", "t_ns": 3500.0, "run_first": 712, "run_last": 718},
    {"label": "3.8 us", "t_ns": 3800.0, "run_first": 838, "run_last": 844},
    {"label": "4.0 us", "t_ns": 4000.0, "run_first": 719, "run_last": 725},
    {"label": "4.3 us", "t_ns": 4300.0, "run_first": 851, "run_last": 857},
    {"label": "4.5 us", "t_ns": 4500.0, "run_first": 726, "run_last": 732},
]

# Optional manual initial guesses for 1D fits (None = auto)
p_guess_X = None
p_guess_y = None
# =======================================================

os.makedirs(out_dir, exist_ok=True)
if use_calibri:
    mpl.rcParams["font.family"] = "Calibri"


# ------------------------ Helpers ------------------------
def _ravel(a):
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()


def extract_param_by_name(ds, name_candidates, units_in_volts=True, debug=False):
    """
    Search ds.get_parameter_data() for a parameter key matching candidates.
    Candidates can be exact names or regex with prefix 're:'.
    Returns array (µV).
    """
    scale = 1e6 if units_in_volts else 1.0
    pdata = ds.get_parameter_data()
    cand = [c.lower() for c in name_candidates]

    for dep_key, block in pdata.items():
        keys = list(block.keys())
        if debug:
            print(f"[dep='{dep_key}'] keys: {keys}")

        for k in keys:
            kl = k.lower().strip()

            if kl in cand:
                v = _ravel(block[k]) * scale
                v = v[np.isfinite(v)]
                if v.size:
                    return v

            for cc in cand:
                if cc.startswith("re:"):
                    pat = cc[3:]
                    if re.search(pat, kl):
                        v = _ravel(block[k]) * scale
                        v = v[np.isfinite(v)]
                        if v.size:
                            return v

    raise RuntimeError(f"Parameter not found. Tried candidates: {name_candidates}")


def extract_X_signed_uV(ds, units_in_volts=True, debug=False):
    # adjust patterns if needed
    candidates = [
        "x",
        "re:sample\\s*x",
        "re:.*:x$",
        "re:.*_x$",
        "re:real",
    ]
    return extract_param_by_name(ds, candidates, units_in_volts=units_in_volts, debug=debug)


def extract_y_uV(ds, units_in_volts=True, debug=False):
    # adjust patterns if needed
    candidates = [
        "y",
        "re:.*:y$",
        "re:.*_y$",
        "re:\\by\\b",
    ]
    return extract_param_by_name(ds, candidates, units_in_volts=units_in_volts, debug=debug)


def make_hist(v, bins=120):
    counts, edges = np.histogram(v, bins=bins)
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


# ------------------------ 1D model ------------------------
def two_gauss(x, *p):
    x0, x1, w1, h1, w2, h2 = p
    return (h1 * np.exp(-(x - x0) ** 2 / (2 * w1**2)) +
            h2 * np.exp(-(x - x1) ** 2 / (2 * w2**2)))


def gauss1(x, x0, w, h):
    return h * np.exp(-(x - x0) ** 2 / (2 * w**2))


def optimal_threshold_equal_priors(x0, x1, w1, h1, w2, h2):
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


def fit_1d_two_gaussian(v, bins, p_guess=None):
    ticks, counts, bw = make_hist(v, bins=bins)

    if p_guess is None:
        thr0 = two_means_threshold(v)
        lo = v[v <= thr0]
        hi = v[v > thr0]
        if lo.size < 10 or hi.size < 10:
            p30, p70 = np.percentile(v, [30, 70])
            lo = v[v <= p30]
            hi = v[v >= p70]

        x0 = np.mean(lo) if lo.size else np.min(v)
        x1 = np.mean(hi) if hi.size else np.max(v)
        if x0 > x1:
            x0, x1 = x1, x0

        w1 = np.std(lo, ddof=1) if lo.size > 1 else np.std(v, ddof=1)
        w2 = np.std(hi, ddof=1) if hi.size > 1 else np.std(v, ddof=1)

        def peak_height(xc):
            i = np.argmin(np.abs(ticks - xc))
            return max(counts[i], 1.0) if counts.size else 1.0

        h1 = peak_height(x0)
        h2 = peak_height(x1)
        p0 = np.array([x0, x1, max(w1, 1e-6), h1, max(w2, 1e-6), h2], dtype=float)
    else:
        p0 = np.array(p_guess, dtype=float)

    lb = [-np.inf, -np.inf, 1e-12, 0.0, 1e-12, 0.0]
    ub = [ np.inf,  np.inf,  np.inf, np.inf,  np.inf, np.inf]
    popt, _ = curve_fit(two_gauss, ticks, counts, p0=p0, bounds=(lb, ub), maxfev=40000)

    x0, x1, w1, h1, w2, h2 = popt
    if x0 > x1:
        x0, x1 = x1, x0
        w1, w2 = w2, w1
        h1, h2 = h2, h1

    return (x0, x1, w1, h1, w2, h2), ticks, counts, bw


def fidelity_from_1d_fit(ticks, popt, grid_points=2_000_001):
    x0, x1, w1, h1, w2, h2 = popt
    thr = optimal_threshold_equal_priors(x0, x1, w1, h1, w2, h2)

    arr = np.linspace(ticks.min(), ticks.max(), grid_points)
    F0 = np.sum(gauss1(arr[arr < thr], x0, w1, h1)) / np.sum(gauss1(arr, x0, w1, h1))
    F1 = np.sum(gauss1(arr[arr > thr], x1, w2, h2)) / np.sum(gauss1(arr, x1, w2, h2))
    return float(0.5 * (F0 + F1)), float(F0), float(F1), float(thr)


def plot_1d_hist_fit_linear(ticks, counts, bw, popt, thr, xlabel, title, out_path=None):
    x0, x1, w1, h1, w2, h2 = popt
    xplot = np.linspace(ticks.min(), ticks.max(), 2000)

    fit_total = two_gauss(xplot, x0, x1, w1, h1, w2, h2)
    g_low = gauss1(xplot, x0, w1, h1)
    g_high = gauss1(xplot, x1, w2, h2)

    plt.figure(figsize=(7.4, 4.8))
    plt.bar(ticks, counts, width=bw, color="0.86", edgecolor="none", alpha=0.95, label="Histogram")
    plt.plot(xplot, fit_total, "-", lw=2.2, label="Fit (G1+G2)")
    plt.plot(xplot, g_low, "--", lw=1.6, label="G1")
    plt.plot(xplot, g_high, "--", lw=1.6, label="G2")
    plt.axvline(thr, ls="--", lw=1.6, color="tab:red", label=f"thr = {thr:.2f} µV")

    plt.xlabel(xlabel, fontsize=13)
    plt.ylabel("Counts", fontsize=13)
    plt.title(title, fontsize=13)
    plt.grid(True, alpha=0.25)
    plt.legend(fontsize=9)
    plt.tight_layout()

    if out_path is not None:
        plt.savefig(out_path, dpi=260)
        print("Saved figure ->", out_path)

    plt.show()


# ------------------------ 2D Bayes fidelity + plotting ------------------------
def gmm_bayes_fidelity_grid(gmm, XY, gridsize=420, percentiles=(0.2, 99.8)):
    """
    Bayes-optimal fidelity (model-based) for a 2-component 2D GMM.
    Computes on a grid:
      decide 0 if w0 p0 >= w1 p1
      F0 = ∫_{R0} p0
      F1 = ∫_{R1} p1
      Fidelity = 0.5*(F0+F1)
    Returns fidelity + grid objects for plotting.
    """
    assert gmm.n_components == 2
    weights = gmm.weights_
    means = gmm.means_
    covs = gmm.covariances_

    x = XY[:, 0]
    y = XY[:, 1]
    xmin, xmax = np.percentile(x, percentiles)
    ymin, ymax = np.percentile(y, percentiles)

    xx = np.linspace(xmin, xmax, gridsize)
    yy = np.linspace(ymin, ymax, gridsize)
    Xg, Yg = np.meshgrid(xx, yy)
    pos = np.dstack((Xg, Yg))

    p0 = multivariate_normal(mean=means[0], cov=covs[0]).pdf(pos)
    p1 = multivariate_normal(mean=means[1], cov=covs[1]).pdf(pos)

    decide0 = (weights[0] * p0) >= (weights[1] * p1)

    dx = (xmax - xmin) / (gridsize - 1)
    dy = (ymax - ymin) / (gridsize - 1)
    dA = dx * dy

    F0 = np.sum(p0[decide0]) * dA
    F1 = np.sum(p1[~decide0]) * dA
    Fidelity = 0.5 * (F0 + F1)

    grid_pack = (Xg, Yg, p0, p1, decide0, xmin, xmax, ymin, ymax)
    return Fidelity, F0, F1, grid_pack


def plot_2d_intensity_gaussian_max_like(gmm, grid_pack, title, out_path=None, cmap="viridis"):
    """
    Produces a continuous intensity image (full axes, no white),
    intensity ~ max component density (or you can swap to mixture density).
    Overlays:
      - white contours for each Gaussian component
      - white dashed Bayes boundary (where w0 p0 == w1 p1)
      - colorbar on the right
    """
    Xg, Yg, p0, p1, decide0, xmin, xmax, ymin, ymax = grid_pack
    w = gmm.weights_

    # Intensity "like your map": use the maximum of the two (weighted) Gaussians
    Z0 = w[0] * p0
    Z1 = w[1] * p1
    Z = np.maximum(Z0, Z1)  # peak/intensity map emphasizing Gaussian maxima

    # Bayes boundary field
    boundary_field = (w[0] * p0) - (w[1] * p1)

    fig, ax = plt.subplots(figsize=(7.2, 5.8))

    im = ax.imshow(
        Z,
        origin="lower",
        extent=[xmin, xmax, ymin, ymax],
        aspect="auto",
        cmap=cmap
    )

    # White contours of each Gaussian "ball"
    for pk in (Z0, Z1):
        # choose a few contour levels as fractions of the peak
        lv = np.array([0.15, 0.35, 0.60, 0.85]) * np.max(pk)
        ax.contour(Xg, Yg, pk, levels=lv, colors="white", linewidths=1.2, alpha=0.95)

    # Bayes boundary in white dashed
    ax.contour(Xg, Yg, boundary_field, levels=[0.0], colors="white", linewidths=2.0, linestyles="--")

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("X (µV)")
    ax.set_ylabel("y (µV)")
    ax.set_title(title)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Intensity (max of weighted 2D Gaussians)")

    fig.tight_layout()

    if out_path is not None:
        fig.savefig(out_path, dpi=260)
        print("Saved figure ->", out_path)

    plt.show()


# ==================== GROUP ANALYSIS ====================
def analyze_group(run_first, run_last, group_label, t_ns):
    """
    Returns: FidX_1D, SNR_X, FidY_1D, Fid2D_Bayes
    """
    print("\n" + "=" * 86)
    print(f"Integration time: {t_ns} ns | group '{group_label}' | runs {run_first}–{run_last}")
    print("=" * 86)

    all_X = []
    all_y = []

    for rid in range(run_first, run_last + 1):
        try:
            ds = load_by_id(rid)
            debug = (rid in (run_first, run_first + 1))

            X = extract_X_signed_uV(ds, units_in_volts=units_in_volts, debug=debug)
            y = extract_y_uV(ds, units_in_volts=units_in_volts, debug=debug)

            # If y is scalar per run, repeat to match X length.
            if y.size == 1 and X.size > 1:
                y = np.full_like(X, float(y.item()))

            n = min(X.size, y.size)
            if n == 0:
                print(f"[run {rid}] empty (skipped)")
                continue

            all_X.append(X[:n])
            all_y.append(y[:n])
            print(f"[run {rid}] collected {n} paired samples (X,y)")

        except Exception as e:
            print(f"[run {rid}] skipped: {e}")

    if not all_X or not all_y:
        raise RuntimeError(f"No paired (X,y) data found for runs {run_first}-{run_last}.")

    X_all = np.concatenate(all_X)
    y_all = np.concatenate(all_y)
    mask = np.isfinite(X_all) & np.isfinite(y_all)
    X_all = X_all[mask]
    y_all = y_all[mask]
    XY = np.column_stack([X_all, y_all])

    print(f"\nAggregated: N = {XY.shape[0]}")
    print(f"X range: [{X_all.min():.2f}, {X_all.max():.2f}] µV")
    print(f"y range: [{y_all.min():.2f}, {y_all.max():.2f}] µV")

    # -------- 1D X
    poptX, ticksX, countsX, bwX = fit_1d_two_gaussian(X_all, bins=bins, p_guess=p_guess_X)
    FidX, _, _, thrX = fidelity_from_1d_fit(ticksX, poptX, grid_points=integration_grid_points)
    x0, x1, w1, h1, w2, h2 = poptX
    SNR_X = 2.0 * (x1 - x0) ** 2 / (w1 ** 2 + w2 ** 2)

    # -------- 1D y
    poptY, ticksY, countsY, bwY = fit_1d_two_gaussian(y_all, bins=bins, p_guess=p_guess_y)
    FidY, _, _, thrY = fidelity_from_1d_fit(ticksY, poptY, grid_points=integration_grid_points)

    # -------- 2D Bayes / overlap fidelity from fitted 2D Gaussians
    gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=0)
    gmm.fit(XY)

    Fid2D, F0_2D, F1_2D, grid_pack = gmm_bayes_fidelity_grid(
        gmm, XY, gridsize=density_gridsize, percentiles=xy_percentiles
    )

    print("\n=== Results ===")
    print(f"Fidelity_X (1D)  = {format_percent_floor(FidX, decimals=fidelity_decimals)}% | SNR_X = {SNR_X:.2f}")
    print(f"Fidelity_y (1D)  = {format_percent_floor(FidY, decimals=fidelity_decimals)}%")
    print(f"Fidelity_2D Bayes(model) = {format_percent_floor(Fid2D, decimals=fidelity_decimals)}%")
    print(f"  (F0_2D = {format_percent_floor(F0_2D, decimals=fidelity_decimals)}%, "
          f"F1_2D = {format_percent_floor(F1_2D, decimals=fidelity_decimals)}%)")

    # -------- Plots
    if save_png:
        tag = f"{group_label.replace(' ','_')}_{int(t_ns)}ns_runs_{run_first}_{run_last}"
        path_X = os.path.join(out_dir, f"hist_X_linear_{tag}.png")
        path_y = os.path.join(out_dir, f"hist_y_linear_{tag}.png")
        path_2D = os.path.join(out_dir, f"map_2D_intensity_{tag}.png")
    else:
        path_X = path_y = path_2D = None

    plot_1d_hist_fit_linear(
        ticksX, countsX, bwX, poptX, thrX,
        xlabel="Amplitude X (µV)",
        title=f"X histogram & 2-Gaussian fit (linear) ({group_label}, {t_ns:.0f} ns)",
        out_path=path_X
    )

    plot_1d_hist_fit_linear(
        ticksY, countsY, bwY, poptY, thrY,
        xlabel="Amplitude y (µV)",
        title=f"y histogram & 2-Gaussian fit (linear) ({group_label}, {t_ns:.0f} ns)",
        out_path=path_y
    )

    plot_2d_intensity_gaussian_max_like(
        gmm, grid_pack,
        title=f"2D intensity (max Gaussian) + white contours + Bayes boundary ({group_label}, {t_ns:.0f} ns)",
        out_path=path_2D,
        cmap=density_cmap
    )

    return FidX, SNR_X, FidY, Fid2D


# ==================== MAIN LOOP ====================
initialise_or_create_database_at(db_path)

all_times_ns = []
all_fidX = []
all_snrX = []
all_fidY = []
all_fid2D = []

for grp in integration_groups:
    FidX, snrX, FidY, Fid2D = analyze_group(
        run_first=grp["run_first"],
        run_last=grp["run_last"],
        group_label=grp.get("label", f"{grp['t_ns']} ns"),
        t_ns=grp["t_ns"],
    )
    all_times_ns.append(grp["t_ns"])
    all_fidX.append(FidX)
    all_snrX.append(snrX)
    all_fidY.append(FidY)
    all_fid2D.append(Fid2D)

all_times_ns = np.array(all_times_ns, float)
all_fidX = np.array(all_fidX, float)
all_snrX = np.array(all_snrX, float)
all_fidY = np.array(all_fidY, float)
all_fid2D = np.array(all_fid2D, float)

# sort by time
order = np.argsort(all_times_ns)
all_times_ns = all_times_ns[order]
all_fidX = all_fidX[order]
all_snrX = all_snrX[order]
all_fidY = all_fidY[order]
all_fid2D = all_fid2D[order]

# ==================== ERROR PLOT (log y) ====================
plt.figure(figsize=(6.4, 4.4))
plt.plot(all_times_ns, 1.0 - all_fidX, "o-", lw=1.8, ms=5, label="1 - Fidelity (X, 1D)")
plt.plot(all_times_ns, 1.0 - all_fidY, "o-", lw=1.8, ms=5, label="1 - Fidelity (y, 1D)")
plt.plot(all_times_ns, 1.0 - all_fid2D, "o-", lw=1.8, ms=5, label="1 - Fidelity (2D Bayes)")
plt.yscale("log")
plt.xlabel("Integration time (ns)")
plt.ylabel("1 - Fidelity")
plt.title("Readout error vs integration time")
plt.grid(True, which="both", ls="--", alpha=0.45)
plt.legend()
plt.tight_layout()

if save_png:
    fig_path = os.path.join(out_dir, "error_vs_integration_time_X_y_2D_Bayes.png")
    plt.savefig(fig_path, dpi=260)
    print("Saved figure ->", fig_path)

plt.show()

# ==================== TEXT SUMMARY ====================
print("\n=== SUMMARY BY INTEGRATION TIME ===")
for t, Fx, Fy, F2, snr in zip(all_times_ns, all_fidX, all_fidY, all_fid2D, all_snrX):
    errx = 1.0 - Fx
    erry = 1.0 - Fy
    err2 = 1.0 - F2
    print(
        f"t = {t:7.1f} ns | "
        f"FidX = {format_percent_floor(Fx, decimals=fidelity_decimals)}% "
        f"(1-F={format_scientific_prob(errx, decimals=fidelity_decimals)}) | "
        f"FidY = {format_percent_floor(Fy, decimals=fidelity_decimals)}% "
        f"(1-F={format_scientific_prob(erry, decimals=fidelity_decimals)}) | "
        f"Fid2D(Bayes) = {format_percent_floor(F2, decimals=fidelity_decimals)}% "
        f"(1-F={format_scientific_prob(err2, decimals=fidelity_decimals)}) | "
        f"SNR_X = {snr:.2f}"
    )
