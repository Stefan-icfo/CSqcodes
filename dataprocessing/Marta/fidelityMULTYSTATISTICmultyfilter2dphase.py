
import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from qcodes.dataset import load_by_id, initialise_or_create_database_at
from sklearn.mixture import GaussianMixture
from scipy.stats import multivariate_normal


# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v38_28_11_25.db"

use_calibri = True
units_in_volts = True  # True if stored in V; convert to µV
out_dir = r"C:\Users\Public\Fidelity_2D_FIRST"
save_png = True

# 2D map settings
density_gridsize = 420                  # higher = smoother, slower
xy_percentiles = (0.2, 99.8)            # clipping to avoid outliers
density_cmap = "viridis"                # "magma"/"inferno" are also nice

# Fit settings
gmm_covariance_type = "full"            # keep "full" for real ellipses
gmm_random_state = 0

# Integration-time groups
integration_groups = [
    {"label": "29.9 ns", "t_ns": 29.9, "run_first": 528, "run_last": 534},
    {"label": "49.9 ns", "t_ns": 49.9, "run_first": 535, "run_last": 541},
    {"label": "102.6 ns", "t_ns": 102.6, "run_first": 542, "run_last": 548},
    # ... add the rest like before ...
]
# =======================================================

os.makedirs(out_dir, exist_ok=True)
if use_calibri:
    mpl.rcParams["font.family"] = "Calibri"


# ------------------------ Extract helpers ------------------------
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
                    if re.search(cc[3:], kl):
                        v = _ravel(block[k]) * scale
                        v = v[np.isfinite(v)]
                        if v.size:
                            return v

    raise RuntimeError(f"Parameter not found. Tried candidates: {name_candidates}")


def extract_X_signed_uV(ds, units_in_volts=True, debug=False):
    candidates = [
        "x",
        "re:sample\\s*x",
        "re:.*:x$",
        "re:.*_x$",
        "re:real",
    ]
    return extract_param_by_name(ds, candidates, units_in_volts=units_in_volts, debug=debug)


def extract_y_uV(ds, units_in_volts=True, debug=False):
    candidates = [
        "y",
        "re:.*:y$",
        "re:.*_y$",
        "re:\\by\\b",
    ]
    return extract_param_by_name(ds, candidates, units_in_volts=units_in_volts, debug=debug)


# ------------------------ 2D model + fidelity ------------------------
def gaussian2d_pdf_grid(Xg, Yg, mean, cov):
    pos = np.dstack((Xg, Yg))
    return multivariate_normal(mean=mean, cov=cov).pdf(pos)


def compute_2d_model_fidelity(gmm, XY, gridsize=420, percentiles=(0.2, 99.8)):
    """
    Bayes-optimal model fidelity for fitted 2-component 2D GMM.
    Defines decision regions by: decide 0 if w0*p0 >= w1*p1.
    Then:
        F0 = ∫_{R0} p0
        F1 = ∫_{R1} p1
        Fidelity = 0.5*(F0+F1)
    """
    assert gmm.n_components == 2
    w = gmm.weights_
    mu = gmm.means_
    cov = gmm.covariances_

    x = XY[:, 0]
    y = XY[:, 1]
    xmin, xmax = np.percentile(x, percentiles)
    ymin, ymax = np.percentile(y, percentiles)

    xx = np.linspace(xmin, xmax, gridsize)
    yy = np.linspace(ymin, ymax, gridsize)
    Xg, Yg = np.meshgrid(xx, yy)

    p0 = gaussian2d_pdf_grid(Xg, Yg, mu[0], cov[0])
    p1 = gaussian2d_pdf_grid(Xg, Yg, mu[1], cov[1])

    decide0 = (w[0] * p0) >= (w[1] * p1)

    dx = (xmax - xmin) / (gridsize - 1)
    dy = (ymax - ymin) / (gridsize - 1)
    dA = dx * dy

    F0 = np.sum(p0[decide0]) * dA
    F1 = np.sum(p1[~decide0]) * dA
    Fid = 0.5 * (F0 + F1)

    grid_pack = (Xg, Yg, p0, p1, decide0, xmin, xmax, ymin, ymax)
    return Fid, F0, F1, grid_pack


def plot_2d_continuous_intensity(gmm, grid_pack, title, out_path=None, cmap="viridis"):
    """
    Continuous 2D intensity map WITHOUT white gaps:
      intensity = max(weighted p0, weighted p1)  (peak-map look)
      contours (white) of both components
      Bayes boundary (white dashed)
      colorbar on right
    """
    Xg, Yg, p0, p1, decide0, xmin, xmax, ymin, ymax = grid_pack
    w = gmm.weights_

    Z0 = w[0] * p0
    Z1 = w[1] * p1
    Z = np.maximum(Z0, Z1)                # peak-map look (like your first figure)

    boundary_field = (w[0] * p0) - (w[1] * p1)

    fig, ax = plt.subplots(figsize=(7.4, 5.8))
    im = ax.imshow(
        Z,
        origin="lower",
        extent=[xmin, xmax, ymin, ymax],
        aspect="auto",
        cmap=cmap,
        interpolation="nearest"
    )

    # white contours of each Gaussian ("balls")
    for Zk in (Z0, Z1):
        levels = np.array([0.15, 0.35, 0.60, 0.85]) * np.max(Zk)
        ax.contour(Xg, Yg, Zk, levels=levels, colors="white", linewidths=1.2, alpha=0.95)

    # Bayes boundary w0*p0 = w1*p1
    ax.contour(Xg, Yg, boundary_field, levels=[0.0], colors="white",
               linewidths=2.0, linestyles="--")

    ax.set_xlabel("X (µV)")
    ax.set_ylabel("y (µV)")
    ax.set_title(title)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Intensity = max(weighted 2D Gaussians)")

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=260)
        print("Saved figure ->", out_path)
    plt.show()


def format_percent(p, decimals=12):
    p = max(0.0, min(1.0, float(p)))
    return f"{p*100:.{decimals}f}"


def format_sci(p, decimals=6):
    p = max(0.0, float(p))
    return f"{p:.{decimals}e}"


# ==================== ANALYSIS FOR ONE GROUP (2D FIRST) ====================
def analyze_group_2d_first(run_first, run_last, group_label, t_ns):
    print("\n" + "=" * 86)
    print(f"Integration time: {t_ns} ns  |  group '{group_label}'  |  runs {run_first}–{run_last}")
    print("=" * 86)

    all_X = []
    all_y = []

    for rid in range(run_first, run_last + 1):
        try:
            ds = load_by_id(rid)
            debug = (rid in (run_first, run_first + 1))

            X = extract_X_signed_uV(ds, units_in_volts=units_in_volts, debug=debug)
            y = extract_y_uV(ds, units_in_volts=units_in_volts, debug=debug)

            # If y is scalar per run, repeat it to match X length
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

    if not all_X:
        raise RuntimeError("No paired data collected.")

    X_all = np.concatenate(all_X)
    y_all = np.concatenate(all_y)
    mask = np.isfinite(X_all) & np.isfinite(y_all)
    X_all = X_all[mask]
    y_all = y_all[mask]
    XY = np.column_stack([X_all, y_all])

    print(f"\nAggregated N = {XY.shape[0]}")
    print(f"X range: [{X_all.min():.2f}, {X_all.max():.2f}] µV")
    print(f"y range: [{y_all.min():.2f}, {y_all.max():.2f}] µV")

    # ---- Fit two 2D Gaussians directly (GMM)
    gmm = GaussianMixture(
        n_components=2,
        covariance_type=gmm_covariance_type,
        random_state=gmm_random_state
    )
    gmm.fit(XY)

    # ---- Compute 2D Bayes threshold + model fidelity (integral on grid)
    Fid2D, F0_2D, F1_2D, grid_pack = compute_2d_model_fidelity(
        gmm, XY, gridsize=density_gridsize, percentiles=xy_percentiles
    )

    print("\n=== 2D-FIRST Results (from 2D Gaussian model) ===")
    print(f"F0_2D = {format_percent(F0_2D, decimals=15)}%")
    print(f"F1_2D = {format_percent(F1_2D, decimals=15)}%")
    print(f"Fidelity_2D = {format_percent(Fid2D, decimals=15)}%")
    print(f"1 - F = {format_sci(1.0 - Fid2D, decimals=15)}")

    # ---- Plot 2D continuous intensity map + white Gaussian contours + boundary
    if save_png:
        tag = f"{group_label.replace(' ', '_')}_{int(t_ns)}ns_runs_{run_first}_{run_last}"
        out_path = os.path.join(out_dir, f"2D_first_map_{tag}.png")
    else:
        out_path = None

    plot_2d_continuous_intensity(
        gmm,
        grid_pack,
        title=f"2D map + 2D Gaussian fit + Bayes boundary ({group_label}, {t_ns:.1f} ns)",
        out_path=out_path,
        cmap=density_cmap
    )

    return float(Fid2D)


# ==================== MAIN LOOP ====================
initialise_or_create_database_at(db_path)

all_times = []
all_F2D = []

for grp in integration_groups:
    F2D = analyze_group_2d_first(
        run_first=grp["run_first"],
        run_last=grp["run_last"],
        group_label=grp["label"],
        t_ns=grp["t_ns"],
    )
    all_times.append(grp["t_ns"])
    all_F2D.append(F2D)

all_times = np.array(all_times, float)
all_F2D = np.array(all_F2D, float)

order = np.argsort(all_times)
all_times = all_times[order]
all_F2D = all_F2D[order]

# Plot 1-F vs integration time
plt.figure(figsize=(6.6, 4.5))
plt.plot(all_times, 1.0 - all_F2D, "o-", lw=1.8, ms=6)
plt.yscale("log")
plt.xlabel("Integration time (ns)")
plt.ylabel("1 - Fidelity (2D)")
plt.title("2D-first model error vs integration time")
plt.grid(True, which="both", ls="--", alpha=0.45)
plt.tight_layout()
if save_png:
    p = os.path.join(out_dir, "error_vs_time_2D_first.png")
    plt.savefig(p, dpi=260)
    print("Saved figure ->", p)
plt.show()

print("\n=== SUMMARY BY INTEGRATION TIME (2D-first) ===")
for t, F in zip(all_times, all_F2D):
    print(f"t = {t:7.1f} ns | Fidelity_2D = {format_percent(F, decimals=15)}% | 1-F = {format_sci(1.0-F, decimals=15)}")
