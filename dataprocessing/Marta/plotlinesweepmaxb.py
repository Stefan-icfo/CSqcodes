# -*- coding: utf-8 -*-
# Ridge stitching + jump detection (x>0) + piecewise-linear smoothing
# Reports jump spacing along x (mV) and jump amplitude along y (mV).

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import qcodes as qc

# ========================== USER SETTINGS ==========================
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"
RUN_ID_1 = 1946          # first dataset id
RUN_ID_2 = None          # optional second dataset id (or None)
PREFERRED_PARAM = "G"    # parameter to pivot to a grid

# Jump-detector knobs
SMOOTH_WINDOW = 3        # odd integer; lower -> more sensitive to small wiggles
MAD_K = 0.8              # robust threshold (in MAD units of residual); 0.6–1.2 is typical
MIN_SEP_MV_X = None      # min spacing between jumps along x, in mV; None = auto-estimate
MIN_PEAK_AMP_MV_Y = 0.18 # min absolute jump amplitude along y, in mV (you said ~0.20–0.25 mV)
# ===================================================================


# ------------------------- grid & ridge helpers -------------------------
def maxima_trace(x_vec, y_vec, G2d):
    """Row-wise maxima: for each row in G2d, return x_at_max, the row's y, and the max value."""
    if G2d.ndim != 2:
        raise ValueError("G must be 2D.")
    valid_rows = np.isfinite(G2d).any(axis=1)
    y_valid = y_vec[valid_rows]
    G_valid = G2d[valid_rows, :]
    G_safe  = np.where(np.isfinite(G_valid), G_valid, -np.inf)
    idx_max = np.argmax(G_safe, axis=1)
    x_at_max = x_vec[idx_max]
    gmax = G_safe[np.arange(G_safe.shape[0]), idx_max]
    return x_at_max, y_valid, gmax

def grid_from_dataset(ds, preferred_param="G"):
    """Build a regular grid (pivot) from a QCoDeS dataset, guessing x/y setpoints."""
    groups = ds.get_parameter_data()
    if preferred_param in groups:
        pname = preferred_param
    else:
        if not groups:
            raise ValueError("No parameters found in dataset.")
        pname = next(iter(groups.keys()))
    grp = groups[pname]

    setpoints = [k for k in grp.keys() if k != pname]
    if len(setpoints) < 1:
        raise ValueError("Expected >=1 setpoint.")
    if len(setpoints) == 1:
        x_name = y_name = setpoints[0]
    else:
        sp1, sp2 = setpoints[:2]
        v1 = np.asarray(grp[sp1]); v2 = np.asarray(grp[sp2])
        if len(np.unique(v1)) >= len(np.unique(v2)):
            x_name, y_name = sp1, sp2
        else:
            x_name, y_name = sp2, sp1

    df = pd.DataFrame({
        "x":   np.asarray(grp[x_name]),
        "y":   np.asarray(grp[y_name]),
        "val": np.asarray(grp[pname]),
    })
    piv = df.pivot_table(index="y", columns="x", values="val", aggfunc="mean", sort=True)
    x = piv.columns.values
    y = piv.index.values
    G = piv.to_numpy()

    # ensure y increases upward
    if len(y) > 1 and y[0] > y[-1]:
        y = y[::-1]
        G = np.flipud(G)
    return pname, x_name, y_name, x, y, G

def stitch_ridges(y_list, x_list, tol=1e-9):
    """Concatenate multiple (y, x_at_max) series, sort by y, remove nearly-duplicate y's."""
    y_all = np.concatenate(y_list)
    x_all = np.concatenate(x_list)
    idx = np.argsort(y_all)
    y_all = y_all[idx]
    x_all = x_all[idx]
    keep = np.ones_like(y_all, dtype=bool)
    if len(y_all) > 1:
        dy = np.diff(y_all)
        keep[1:] = np.abs(dy) > tol
    return y_all[keep], x_all[keep]


# ------------------------- jump detection helpers ------------------------
def moving_average(a, w):
    """Centered moving average with odd window length."""
    w = int(max(1, w))
    if w % 2 == 0:
        w += 1
    k = np.ones(w) / w
    return np.convolve(a, k, mode="same")

def _local_maxima(arr):
    return np.where((arr[1:-1] > arr[:-2]) & (arr[1:-1] > arr[2:]))[0] + 1
def _local_minima(arr):
    return np.where((arr[1:-1] < arr[:-2]) & (arr[1:-1] < arr[2:]))[0] + 1

def _autocorr_period_mV(x, r, min_points=4):
    """Coarse period estimate of residual r(x) in mV along x."""
    r = r - np.mean(r)
    ac = np.correlate(r, r, mode="full")
    ac = ac[ac.size // 2:]
    peaks = _local_maxima(ac)
    peaks = peaks[peaks >= min_points]
    if peaks.size == 0:
        return None
    dx = np.median(np.diff(x))
    return 1000.0 * peaks[0] * dx  # mV

def detect_jumps(x, y, *,
                 x_positive_only=True,
                 smooth_window=3,
                 mad_k=0.8,
                 min_sep_mV_x=None,
                 min_peak_amp_mV_y=0.18):
    """
    Detect jumps as local extremes (peaks + troughs) of the detrended residual.
    - min_sep_mV_x: min spacing between jumps along x (mV). If None, auto-estimate.
    - min_peak_amp_mV_y: min absolute amplitude along y (mV) to accept a jump.
    Returns: jump_x (V), jump_y_amp_mV (absolute residual amplitude in mV), and info dict.
    """
    mask = np.isfinite(x) & np.isfinite(y)
    if x_positive_only:
        mask &= (x > 0)
    x = x[mask]; y = y[mask]
    if x.size < 5:
        return np.array([]), np.array([]), {"x": x, "y": y}

    # linear detrend
    p = np.polyfit(x, y, 1)
    y_base = np.polyval(p, x)
    resid = y - y_base
    r = moving_average(resid, smooth_window)

    # candidates: peaks AND troughs
    peaks  = _local_maxima(r)
    trough = _local_minima(r)

    # robust symmetric threshold (MAD)
    mad = np.median(np.abs(r - np.median(r))) + 1e-12
    thr_mad = mad_k * mad

    # amplitude threshold on y (convert mV -> V)
    thr_amp = (min_peak_amp_mV_y / 1000.0)

    pos_ok = [i for i in peaks  if (r[i] - np.median(r)) >  max(thr_mad, thr_amp)]
    neg_ok = [i for i in trough if (np.median(r) - r[i]) > max(thr_mad, thr_amp)]
    cand = np.array(sorted(pos_ok + neg_ok), dtype=int)

    # minimum spacing along x
    if min_sep_mV_x is None:
        est = _autocorr_period_mV(x, r)
        min_sep_mV_x = 0.6*est if est else 30.0
    min_sep_V = min_sep_mV_x / 1000.0

    jump_idx = []
    for i in cand:
        if not jump_idx or (x[i] - x[jump_idx[-1]] >= min_sep_V):
            jump_idx.append(i)
    jump_idx = np.array(jump_idx, dtype=int)

    jump_x = x[jump_idx]
    jump_amp_mV = np.abs(r[jump_idx]) * 1000.0  # amplitude along y, in mV

    info = {"x": x, "y": y, "resid": r, "y_base": y_base, "jump_idx": jump_idx}
    return jump_x, jump_amp_mV, info

def piecewise_linear_between_jumps(x, y, jump_x):
    """Piecewise linear fit of y(x) within each interval between jumps."""
    if x.size < 2:
        return np.copy(y)
    if jump_x.size == 0:
        c = np.polyfit(x, y, 1)
        return np.polyval(c, x)

    # segment edges = midpoints between consecutive jumps + endpoints
    edges = [x.min()] + list(0.5 * (jump_x[:-1] + jump_x[1:])) + [x.max()]
    y_fit = np.empty_like(y)
    for a, b in zip(edges[:-1], edges[1:]):
        sel = (x >= a) & (x <= b)
        if sel.sum() >= 2:
            c = np.polyfit(x[sel], y[sel], 1)
            y_fit[sel] = np.polyval(c, x[sel])
        else:
            y_fit[sel] = y[sel]
    return y_fit


# =============================== MAIN ====================================
def main():
    qc.config["core"]["db_location"] = DB_PATH
    _ = qc.experiments()

    # --- Load first run and build ridge ---
    ds1 = qc.load_by_id(RUN_ID_1)
    p1, x1_name, y1_name, x1, y1, G1 = grid_from_dataset(ds1, preferred_param=PREFERRED_PARAM)
    x1_max, y1_valid, _ = maxima_trace(x1, y1, G1)

    ys_list = [y1_valid]
    xs_list = [x1_max]
    x_label = y1_name   # x-axis (sweep)
    y_label = x1_name   # y-axis (channel at the crest, e.g., ch06)

    # --- Optional second run ---
    if RUN_ID_2 is not None:
        try:
            ds2 = qc.load_by_id(RUN_ID_2)
            p2, x2_name, y2_name, x2, y2, G2 = grid_from_dataset(ds2, preferred_param=PREFERRED_PARAM)
            x2_max, y2_valid, _ = maxima_trace(x2, y2, G2)
            ys_list.append(y2_valid)
            xs_list.append(x2_max)
            if y2_name != y1_name:
                print(f"[Warn] Sweep axis differs: '{y1_name}' vs '{y2_name}'. Proceeding anyway.")
            if x2_name != x1_name:
                print(f"[Warn] Transverse axis differs: '{x1_name}' vs '{x2_name}'. Proceeding anyway.")
        except Exception as e:
            print(f"[Info] Could not load RUN_ID_2={RUN_ID_2}: {e}")

    # --- Stitch into a single curve (x = sweep; y = ridge value) ---
    x_sweep, y_ridge = stitch_ridges(ys_list, xs_list, tol=1e-9)

    # --- Detect jumps (x>0) using your amplitude hint on ch06 ---
    jump_x, jump_amp_mV, info = detect_jumps(
        x_sweep, y_ridge,
        x_positive_only=True,
        smooth_window=SMOOTH_WINDOW,
        mad_k=MAD_K,
        min_sep_mV_x=MIN_SEP_MV_X,
        min_peak_amp_mV_y=MIN_PEAK_AMP_MV_Y,
    )

    # --- Report spacing along x and amplitude along y ---
    if jump_x.size >= 2:
        spacing_mV = np.diff(jump_x) * 1000.0
        print("Jump positions on x (mV):", np.round(jump_x * 1000.0, 3))
        print("Jump spacing on x  (mV):", np.round(spacing_mV, 3))
        print(f"Mean spacing = {np.mean(spacing_mV):.3f} mV | Median = {np.median(spacing_mV):.3f} mV | Std = {np.std(spacing_mV, ddof=1):.3f} mV")
        print("Jump amplitudes on y (mV):", np.round(jump_amp_mV, 3))
    elif jump_x.size == 1:
        print(f"Only one jump detected at x = {jump_x[0]*1000.0:.3f} mV (amplitude on y = {jump_amp_mV[0]:.3f} mV).")
    else:
        print("No jumps detected. Try lowering MAD_K or MIN_PEAK_AMP_MV_Y, or set MIN_SEP_MV_X manually.")

    # --- Piecewise-linear reconstruction between jumps ---
    y_piecewise = piecewise_linear_between_jumps(info["x"], info["y"], jump_x)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(x_sweep, y_ridge, color="k", lw=1.2, label="Stitched crest (raw)")
    ax.plot(info["x"], y_piecewise, color="tab:orange", lw=2.2, label="Piecewise-linear between jumps")
    # mark jumps
    y_at_jumps = np.interp(jump_x, x_sweep, y_ridge) if jump_x.size else []
    for xj, yj in zip(jump_x, y_at_jumps):
        ax.plot([xj], [yj], "o", ms=8, mfc="none", mec="red", mew=2)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title("Crest of G — stitched; jumps (x>0) and piecewise-linear smoothing")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()

