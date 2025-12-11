import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ===================== USER SETTINGS =====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v44_05_12_25.db"
run_id = 20
SIGNAL_PARAM = "x"

# Units
UNITS_IN_VOLTS = True  # True if x stored in V (common), False if already in µV

# Start index: look backwards from edge for the first sample x in this band (µV)
START_BAND_UV = (-1450.0, -1200.0)

# Fit window timing
W_PRE_US  = 8.0     # how far back we search (and can include) before the edge
W_POST_US = 60.0    # how long after the edge we include in the fit window

# Minimum points for a fit window
MIN_POINTS_FIT = 6

# Keep/reject fits
MIN_R2 = 0.30  # with few points, keep loose

# Quality filter: keep only edges with enough points on the transition (avoid "vertical" edges)
TRANSITION_LOW_UV  = -1200.0
TRANSITION_HIGH_UV = +1200.0
MIN_TRANSITION_PTS = 4   # increase -> stricter selection

# Reject windows that are mostly already on the high plateau
PLATEAU_HIGH_UV = 1300.0
MAX_PLATEAU_FRAC_IN_WINDOW = 0.75

# Plots
PLOT_DEBUG_EDGES = True
PLOT_ALL_FITS_ON_FULL_TRACE = True
PLOT_ONE_ZOOM = True
ZOOM_EDGE_INDEX = 10  # 0..N-1
# =========================================================


def load_time_trace_from_db(db_path, run_id, signal_param):
    initialise_or_create_database_at(db_path)
    ds = load_by_id(run_id)
    print(f"Loaded run_id = {run_id} from DB:\n{ds}\n")

    data_dict = ds.get_parameter_data(signal_param)
    inner = data_dict[signal_param]

    time_keys = [k for k in inner.keys() if k != signal_param]
    if not time_keys:
        raise RuntimeError("Could not find time parameter in dataset.")
    time_name = time_keys[0]

    t = np.asarray(inner[time_name]).ravel()
    v = np.asarray(inner[signal_param]).ravel()

    o = np.argsort(t)
    t, v = t[o], v[o]

    print(f"Time parameter:   {time_name}")
    print(f"Signal parameter: {signal_param}")
    print(f"Number of points: {len(t)}")
    return t, v


def exp_rise(t, A, tau, C):
    """Exponential rise: C + A*(1-exp(-t/tau))"""
    return C + A * (1.0 - np.exp(-t / tau))


def r2_score(y, yfit):
    ss_res = np.sum((y - yfit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 1.0


def find_rising_edges(t, v):
    """
    Robust rising-edge finder for square-like traces:
    - low/high from percentiles
    - threshold halfway
    - rising edges where boolean crosses 0->1
    """
    v_low = np.percentile(v, 5)
    v_high = np.percentile(v, 95)
    thr = 0.5 * (v_low + v_high)

    high = v > thr
    edges = np.where((~high[:-1]) & (high[1:]))[0] + 1

    dt = np.median(np.diff(t)) if len(t) > 2 else np.nan

    print("\n--- Edge finding ---")
    print(f"v_low(5%)   = {v_low*1e6:.1f} µV")
    print(f"v_high(95%) = {v_high*1e6:.1f} µV")
    print(f"threshold   = {thr*1e6:.1f} µV")
    print(f"dt          = {dt*1e6:.3f} µs")
    print(f"rising edges found = {len(edges)}")

    return edges, thr, dt


def choose_start_index_by_band(t, v, edge_idx, dt, w_pre_us, band_uv, units_in_volts=True):
    """
    Search backwards from edge_idx over W_PRE_US and pick the earliest sample
    that falls inside band_uv (in µV). If none found, fallback to edge_idx - pre.
    """
    scale = 1e6 if units_in_volts else 1.0
    lo, hi = band_uv
    pre = int((w_pre_us * 1e-6) / dt)
    i0 = max(0, edge_idx - pre)

    vv_uv = v[i0:edge_idx + 1] * scale
    inside = np.where((vv_uv >= lo) & (vv_uv <= hi))[0]

    if inside.size > 0:
        # earliest index within this pre-window that hits the band
        return i0 + int(inside[0])

    return i0


def edge_is_good(v_segment, units_in_volts=True):
    """
    Quality filter:
    - requires enough points in a wide "transition band"
    - rejects windows that are mostly already on the high plateau
    """
    scale = 1e6 if units_in_volts else 1.0
    vv = v_segment * scale

    in_trans = np.sum((vv >= TRANSITION_LOW_UV) & (vv <= TRANSITION_HIGH_UV))
    if in_trans < MIN_TRANSITION_PTS:
        return False

    frac_plateau = np.mean(vv >= PLATEAU_HIGH_UV)
    if frac_plateau > MAX_PLATEAU_FRAC_IN_WINDOW:
        return False

    return True


def fit_one_edge_with_custom_start(t, v, edge_idx, dt, i_start, w_post_us):
    """
    Fit from i_start to edge_idx + W_POST_US.
    Time is referenced to i_start (t=0 at start).
    """
    post = int((w_post_us * 1e-6) / dt)
    i_end = min(len(v), edge_idx + post)

    if i_end - i_start < MIN_POINTS_FIT:
        return None

    tt = t[i_start:i_end] - t[i_start]
    vv = v[i_start:i_end]

    # Initial guesses
    n0 = min(5, len(vv))
    C0 = float(np.median(vv[:n0]))
    plateau = float(np.median(vv[-max(3, len(vv)//5):]))
    A0 = plateau - C0
    tau0 = max(dt, 0.2 * (tt.max() - tt.min()))

    bounds = ([-np.inf, 0.0, -np.inf], [np.inf, np.inf, np.inf])  # tau >= 0

    popt, _ = curve_fit(
        exp_rise, tt, vv,
        p0=[A0, tau0, C0],
        bounds=bounds,
        maxfev=50000
    )

    vfit = exp_rise(tt, *popt)
    r2 = r2_score(vv, vfit)

    return {
        "i_start": i_start,
        "i_end": i_end,
        "tt": tt,
        "vv": vv,
        "vfit": vfit,
        "popt": popt,
        "r2": r2
    }


def main():
    t, v = load_time_trace_from_db(db_path, run_id, SIGNAL_PARAM)
    edges, thr, dt = find_rising_edges(t, v)

    if not np.isfinite(dt) or dt <= 0:
        raise RuntimeError("Bad dt; check time array.")

    if PLOT_DEBUG_EDGES:
        plt.figure(figsize=(10, 4))
        plt.plot(t * 1e6, v * 1e6, "k-", lw=1)
        plt.axhline(thr * 1e6, ls="--", lw=1)
        plt.plot(t[edges] * 1e6, v[edges] * 1e6, "ro", ms=5)
        plt.xlabel("time (µs)")
        plt.ylabel("x (µV)")
        plt.title("Debug: rising edges (red) and threshold (dashed)")
        plt.tight_layout()
        plt.show()

    kept_taus = []
    kept_r2 = []
    kept_results = []

    if PLOT_ALL_FITS_ON_FULL_TRACE:
        plt.figure(figsize=(10, 4))
        plt.plot(t * 1e6, v * 1e6, "k-", lw=1, label="data")

    print("\n--- Fitting ---")
    for k, eidx in enumerate(edges):
        # choose earlier start based on your requested low-band
        i_start = choose_start_index_by_band(
            t, v, eidx, dt,
            w_pre_us=W_PRE_US,
            band_uv=START_BAND_UV,
            units_in_volts=UNITS_IN_VOLTS
        )

        # build candidate window to apply quality filter
        post = int((W_POST_US * 1e-6) / dt)
        i_end = min(len(v), eidx + post)

        if i_end - i_start < MIN_POINTS_FIT:
            continue

        if not edge_is_good(v[i_start:i_end], units_in_volts=UNITS_IN_VOLTS):
            continue

        res = fit_one_edge_with_custom_start(t, v, eidx, dt, i_start, W_POST_US)
        if res is None:
            continue

        A, tau, C = res["popt"]
        r2 = res["r2"]

        print(f"edge {k:02d}: "
              f"start@{t[i_start]*1e6:.2f} µs, end@{t[res['i_end']-1]*1e6:.2f} µs | "
              f"tau={tau*1e6:.3f} µs | R2={r2:.3f} | N={res['i_end']-res['i_start']}")

        if r2 < MIN_R2:
            continue

        kept_taus.append(tau)
        kept_r2.append(r2)
        kept_results.append((k, eidx, res))

        if PLOT_ALL_FITS_ON_FULL_TRACE:
            t_abs = (t[res["i_start"]] + res["tt"]) * 1e6
            plt.plot(t_abs, res["vfit"] * 1e6, lw=2)

    if PLOT_ALL_FITS_ON_FULL_TRACE:
        plt.xlabel("time (µs)")
        plt.ylabel("x (µV)")
        plt.title(f"Custom-start exp-rise fits | kept {len(kept_taus)} (R2>{MIN_R2})")
        plt.tight_layout()
        plt.show()

    if PLOT_ONE_ZOOM and kept_results:
        idx = min(ZOOM_EDGE_INDEX, len(kept_results) - 1)
        k, eidx, res = kept_results[idx]

        t_abs = (t[res["i_start"]] + res["tt"]) * 1e6
        plt.figure(figsize=(9, 3.5))
        plt.plot(t_abs, res["vv"] * 1e6, "k-", lw=1, label="data (window)")
        plt.plot(t_abs, res["vfit"] * 1e6, lw=2, label=f"fit (tau={res['popt'][1]*1e6:.3f} µs)")
        plt.xlabel("time (µs)")
        plt.ylabel("x (µV)")
        plt.title(f"Zoom edge {k:02d}: {t[res['i_start']]*1e6:.1f}→{t[res['i_end']-1]*1e6:.1f} µs")
        plt.legend()
        plt.tight_layout()
        plt.show()

    if kept_taus:
        taus_us = np.array(kept_taus) * 1e6
        print("\n=== SUMMARY (kept fits) ===")
        print(f"Mean tau  = {taus_us.mean():.3f} µs")
        print(f"Std tau   = {taus_us.std():.3f} µs")
        print(f"N fits    = {len(taus_us)}")
        print(f"Median R2 = {np.median(kept_r2):.3f}")
    else:
        print("\nNo fits kept.")
        print("Try:")
        print(" - decrease MIN_TRANSITION_PTS (e.g. 3)")
        print(" - increase W_POST_US (e.g. 100 us)")
        print(" - lower PLATEAU_HIGH_UV (e.g. 1200)")
        print(" - or check UNITS_IN_VOLTS setting")


if __name__ == "__main__":
    main()
