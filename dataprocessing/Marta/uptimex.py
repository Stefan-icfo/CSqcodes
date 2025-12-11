import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ====== DATABASE PATH AND RUN ID ======
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v44_05_12_25.db"
run_id = 20

# Use the demod X signal stored in the dataset
SIGNAL_PARAM = "x"


# ---------- Load time trace from QCoDes database ----------
def load_time_trace_from_db(db_path, run_id, signal_param=SIGNAL_PARAM):
    """
    Load time and signal arrays from a QCoDeS run.

    Parameters
    ----------
    db_path : str
        Path to the .db file.
    run_id : int
        Run ID of the measurement.
    signal_param : str
        Name of the signal parameter ('x', 'y', 'v_r', 'v_r_avg', ...).

    Returns
    -------
    t : np.ndarray
        Time array (seconds).
    v : np.ndarray
        Signal array (same units as stored in the dataset).
    """
    initialise_or_create_database_at(db_path)
    ds = load_by_id(run_id)
    print(f"Loaded run_id = {run_id} from DB:")
    print(ds)

    data_dict = ds.get_parameter_data(signal_param)
    inner = data_dict[signal_param]

    # Time parameter is any setpoint that is not the signal itself
    time_keys = [k for k in inner.keys() if k != signal_param]
    if len(time_keys) == 0:
        raise RuntimeError("Could not find a time parameter in the dataset.")
    time_name = time_keys[0]

    t = np.asarray(inner[time_name])
    v = np.asarray(inner[signal_param])

    print(f"Time parameter:   {time_name}")
    print(f"Signal parameter: {signal_param}")
    print(f"Number of points: {len(t)}")

    return t, v


# ---------- Exponential rise model ----------
def exp_rise(t, A, tau, C):
    """
    Exponential rise: C + A * (1 - exp(-t / tau)).
    """
    return C + A * (1.0 - np.exp(-t / tau))


# ---------- Find all peaks (high plateaus) ----------
def find_all_peaks(t, v, peak_height_frac=0.4):
    """
    Detect peaks based on a threshold that is a fraction between
    a low baseline and the maximum signal.

    Returns indices of peaks (high plateaus).
    """
    baseline = np.percentile(v, 5)   # low plateau
    v_max = np.max(v)

    height_threshold = baseline + peak_height_frac * (v_max - baseline)

    if len(t) > 1:
        dt = np.median(np.diff(t))
        distance_samples = max(1, int(200e-6 / dt))  # ~30 µs
    else:
        distance_samples = 50

    peaks, _ = find_peaks(v, height=height_threshold,
                          distance=distance_samples)
    print(f"Found {len(peaks)} peaks (high plateaus).")
    return peaks


# ---------- Extract rise segments (low plateau → middle of high plateau) ----------
def extract_rise_segments(t, v,
                          peak_height_frac=0.6,
                          low_frac=0.01,
                          plateau_frac=0.9,
                          min_points=2):
    """
    For each interior peak, define a rise window that starts in the low plateau
    and ends in the middle (in time) of the high plateau.

    Steps:
      - Find peaks with find_all_peaks().
      - For each peak pk:
          * prev_pk = previous peak (or extrapolated)
          * next_pk = next peak (or extrapolated)
          * win_start = midpoint between prev_pk and pk
          * win_end   = midpoint between pk and next_pk
            -> window covers: low plateau → rise → high plateau.
          * In this window, use region BEFORE pk as low plateau.
            - baseline = median around the minimum there.
          * amplitude = v[pk] - baseline
          * level_low     = baseline + low_frac     * amplitude
          * plateau_level = baseline + plateau_frac * amplitude
          * i_start = first index after valley s.t. v >= level_low
          * plateau_indices = all j after pk with v >= plateau_level
            (continuous region on the high plateau)
          * i_end = middle index of plateau_indices
    Returns list of (i_start, i_end) with i_end exclusive.
    """
    peaks = find_all_peaks(t, v, peak_height_frac=peak_height_frac)
    if len(peaks) == 0:
        return []

    rise_segments = []

    for i, pk in enumerate(peaks):
        # Previous peak
        if i == 0:
            if len(peaks) > 1:
                period = peaks[1] - peaks[0]
            else:
                period = 100
            prev_pk = max(0, pk - period)
        else:
            prev_pk = peaks[i - 1]

        # Next peak
        if i < len(peaks) - 1:
            next_pk = peaks[i + 1]
        else:
            if len(peaks) > 1:
                period = peaks[i] - peaks[i - 1]
            else:
                period = 100
            next_pk = min(len(v) - 1, pk + period)

        # Window covering low plateau → rise → high plateau
        win_start = (prev_pk + pk) // 2
        win_end = (pk + next_pk) // 2

        if win_end - win_start < min_points + 4:
            continue

        v_win = v[win_start:win_end + 1]
        idx_pk_rel = pk - win_start

        # Low plateau region BEFORE the peak
        v_low = v_win[:idx_pk_rel]
        if len(v_low) < 5:
            continue

        # Baseline = median around minimum in low plateau
        i_min_rel_center = np.argmin(v_low)
        i_base_start = max(0, i_min_rel_center - 2)
        i_base_end = min(len(v_low), i_min_rel_center + 3)
        baseline = np.median(v_low[i_base_start:i_base_end])

        amp = v[pk] - baseline
        if amp <= 0:
            continue

        level_low = baseline + low_frac * amp
        plateau_level = baseline + plateau_frac * amp

        # ---- Start index: first point >= level_low (after valley) ----
        i_start_rel = None
        for j in range(i_min_rel_center, len(v_win)):
            if v_win[j] >= level_low:
                i_start_rel = j
                break
        if i_start_rel is None:
            continue

        # ---- High plateau region: indices >= plateau_level after the peak ----
        plateau_indices_rel = []
        for j in range(idx_pk_rel, len(v_win)):
            if v_win[j] >= plateau_level:
                plateau_indices_rel.append(j)
            else:
                # once we drop below plateau_level, we stop (likely the falling edge)
                if plateau_indices_rel:
                    break

        if len(plateau_indices_rel) < min_points:
            continue

        # Middle of the plateau region (in time)
        plateau_mid_rel = (plateau_indices_rel[0] + plateau_indices_rel[-1]) // 2

        # Require enough points between start and end
        if plateau_mid_rel - i_start_rel < min_points:
            continue

        i_start = win_start + i_start_rel
        i_end = win_start + plateau_mid_rel  # end will be exclusive

        rise_segments.append((i_start, i_end))

    print(f"Rise segments found: {len(rise_segments)}")
    return rise_segments


# ---------- Fit exponential rise ----------
def fit_rise(t, v, i1, i2):
    """
    Fit an exponential rise within indices [i1, i2).

    Returns best-fit parameters and R^2.
    """
    t_win = t[i1:i2] - t[i1]
    v_win = v[i1:i2]

    C0 = v_win[0]
    A0 = v_win[-1] - v_win[0]
    tau0 = (t_win[-1] - t_win[0]) / 2 if len(t_win) > 1 else 1e-6

    popt, _ = curve_fit(exp_rise, t_win, v_win,
                        p0=[A0, tau0, C0], maxfev=5000)

    v_fit = exp_rise(t_win, *popt)
    ss_res = np.sum((v_win - v_fit) ** 2)
    ss_tot = np.sum((v_win - np.mean(v_win)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return popt, t_win, v_win, r2


# ---------- MAIN ----------
def main():
    # Load time trace (time_param, x) from database
    t, v = load_time_trace_from_db(db_path, run_id, SIGNAL_PARAM)

    print("Data loaded for RISE analysis (using X signal).")

    # Extract windows for exponential rise fits (up to middle of plateau)
    rise_windows = extract_rise_segments(
    t, v,
    peak_height_frac=0.2,
    low_frac=0.0005,
    plateau_frac=0.7,   # <-- abbassa, plateau più facile da trovare
    min_points=2       # <-- abbassa, perché edge velocissimo
)


    all_taus_up = []

    plt.figure(figsize=(10, 5))
    # x is stored in volts → plot in µV
    plt.plot(t * 1e6, v * 1e6, 'k-', label="data")

    for (i1, i2) in rise_windows:
        popt, tw, vw, r2 = fit_rise(t, v, i1, i2)

        # Reject poor fits based on R^2
        if r2 < 0.9:
            
            continue

        A, tau, C = popt
        all_taus_up.append(tau)

        v_fit = exp_rise(tw, A, tau, C)
        plt.plot(t[i1:i2] * 1e6, v_fit * 1e6, lw=2)

    plt.xlabel("time (µs)")
    plt.ylabel("x (µV)")
    plt.title("Exponential rise fits on demod X (low plateau → mid high plateau)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print("\n=== τ_rise results on X (only good fits) ===")
    for i, tau in enumerate(all_taus_up):
        print(f"Rise {i}: tau_up = {tau * 1e6:.3f} µs")

    if all_taus_up:
        mean_tau = np.mean(all_taus_up) * 1e6
        std_tau = np.std(all_taus_up) * 1e6
        print(f"\nMean tau_up: {mean_tau:.3f} µs")
        print(f"Std  tau_up: {std_tau:.3f} µs")
    else:
        print("\nNo good fits (check thresholds or parameters).")


if __name__ == "__main__":
    main()
