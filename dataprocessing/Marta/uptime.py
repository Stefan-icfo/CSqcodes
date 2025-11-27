import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ====== DATABASE PATH AND RUN ID ======
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v35_26_11_25.db"
run_id = 37

# Name of the signal parameter in the dataset
# Change to "v_r_avg", "x", etc. if needed
SIGNAL_PARAM = "v_r"


# ---------- Load time trace from QCoDeS database ----------
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
        Name of the signal parameter (e.g. 'v_r', 'v_r_avg', 'x').

    Returns
    -------
    t : np.ndarray
        Time array (seconds).
    v : np.ndarray
        Signal array (same units as stored in the dataset).
    """
    # Point QCoDeS to the correct database
    initialise_or_create_database_at(db_path)

    # Load dataset
    ds = load_by_id(run_id)
    print(f"Loaded run_id = {run_id} from DB:")
    print(ds)

    # get_parameter_data returns a nested dict:
    # data_dict[signal_param][signal_param]          -> values
    # data_dict[signal_param][time_parameter_name]   -> setpoints (time)
    data_dict = ds.get_parameter_data(signal_param)
    inner = data_dict[signal_param]

    # Find the time parameter name (anything that is not the signal itself)
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
    Exponential rise: C + A * (1 - exp(-t / tau))

    Parameters
    ----------
    t : np.ndarray
        Time (relative).
    A : float
        Amplitude.
    tau : float
        Time constant.
    C : float
        Offset (baseline).

    Returns
    -------
    np.ndarray
        Exponential rise evaluated at t.
    """
    return C + A * (1.0 - np.exp(-t / tau))


# ---------- Find all peaks ----------
def find_all_peaks(t, v, peak_height_frac=0.6):
    """
    Detect peaks based on a threshold that is a fraction between
    a low baseline and the maximum signal.

    Parameters
    ----------
    t : np.ndarray
        Time array.
    v : np.ndarray
        Signal array.
    peak_height_frac : float
        Fraction between baseline and max used as height threshold.

    Returns
    -------
    peaks : np.ndarray
        Indices of detected peaks.
    """
    # Estimate low baseline from the 5th percentile
    baseline = np.percentile(v, 5)
    v_max = np.max(v)

    height_threshold = baseline + peak_height_frac * (v_max - baseline)

    # Minimum distance between peaks, here ~30 microseconds
    if len(t) > 1:
        dt = np.median(np.diff(t))
        distance_samples = max(1, int(30e-6 / dt))
    else:
        distance_samples = 50

    peaks, _ = find_peaks(v, height=height_threshold,
                          distance=distance_samples)
    print(f"Found {len(peaks)} peaks.")
    return peaks


# ---------- Extract rise segments (plateau → upper plateau) ----------
def extract_rise_segments(t, v,
                          peak_height_frac=0.6,
                          low_frac=0.0001,
                          high_frac=1.1,
                          min_points=8):
    """
    For each peak, define a rise window that starts in the low plateau
    and ends inside the upper plateau.

    Steps:
      - Use peaks from find_all_peaks().
      - For each peak pk:
          * Define prev_pk = previous peak (or extrapolated)
          * Define next_pk = next peak (or extrapolated)
          * Define a window [win_start, win_end] where:
                win_start = midpoint between prev_pk and pk
                win_end   = midpoint between pk and next_pk
            So the window covers: low plateau → rise → high plateau.
          * Inside this window, define a "valley" (minimum) BEFORE the peak:
                valley is argmin(v_win[:idx_pk_rel])
            This is used to estimate the lower plateau baseline.
          * amplitude = v[pk] - baseline
          * low level  = baseline + low_frac  * amplitude
          * high level = baseline + high_frac * amplitude
          * start index = first point after valley >= low level
          * end index   = first point after start >= high level

    Returns
    -------
    rise_segments : list of (i_start, i_end)
        Index intervals [i_start, i_end) for the rise fits.
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

        # Window for the rise:
        # from midpoint of previous plateau to midpoint of next plateau
        win_start = (prev_pk + pk) // 2
        win_end = (pk + next_pk) // 2

        if win_end - win_start < min_points + 4:
            continue

        v_win = v[win_start:win_end + 1]

        # Index of the peak within the window
        idx_pk_rel = pk - win_start

        # Take the minimum BEFORE the peak as the low plateau
        v_low_region = v_win[:idx_pk_rel]
        if len(v_low_region) == 0:
            continue
        i_min_rel = np.argmin(v_low_region)

        # Use a small neighborhood around the valley to estimate baseline
        i_base_start = max(0, i_min_rel - 3)
        i_base_end = min(len(v_low_region), i_min_rel + 4)
        baseline = np.median(v_low_region[i_base_start:i_base_end])

        amp = v[pk] - baseline
        if amp <= 0:
            continue

        level_low = baseline + low_frac * amp
        level_high = baseline + high_frac * amp

        # Start: first point after the valley that reaches level_low
        i_start_rel = None
        for j in range(i_min_rel, len(v_win)):
            if v_win[j] >= level_low:
                i_start_rel = j
                break
        if i_start_rel is None:
            continue

        # End: first point after start that reaches level_high
        i_end_rel = None
        for j in range(i_start_rel, len(v_win)):
            if v_win[j] >= level_high:
                i_end_rel = j
                break
        if i_end_rel is None:
            continue

        if i_end_rel - i_start_rel < min_points:
            continue

        i_start = win_start + i_start_rel
        i_end = win_start + i_end_rel + 1  # end is exclusive

        rise_segments.append((i_start, i_end))

    print(f"Rise segments found: {len(rise_segments)}")
    return rise_segments


# ---------- Fit exponential rise ----------
def fit_rise(t, v, i1, i2):
    """
    Fit an exponential rise within indices [i1, i2).

    Parameters
    ----------
    t : np.ndarray
        Time array.
    v : np.ndarray
        Signal array.
    i1, i2 : int
        Start and end indices (end is exclusive).

    Returns
    -------
    popt : (A, tau, C)
        Best-fit parameters.
    t_win : np.ndarray
        Time window (relative, t - t[i1]).
    v_win : np.ndarray
        Signal window.
    r2 : float
        Coefficient of determination (goodness of fit).
    """
    t_win = t[i1:i2] - t[i1]
    v_win = v[i1:i2]

    C0 = v_win[0]                       # initial offset ~ low plateau
    A0 = v_win[-1] - v_win[0]           # amplitude guess
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
    # Load time trace from database
    t, v = load_time_trace_from_db(db_path, run_id, SIGNAL_PARAM)

    print("Data loaded for RISE analysis.")

    # Extract windows for exponential rise fits
    rise_windows = extract_rise_segments(
        t, v,
        peak_height_frac=0.4,
        low_frac=0.0001,   # start a bit above low plateau
        high_frac=1,  # stop deep into the high plateau
        min_points=8
    )

    all_taus_up = []

    plt.figure(figsize=(10, 5))
    plt.plot(t * 1e6, v * 1e3, 'k-', label="data")

    for (i1, i2) in rise_windows:
        popt, tw, vw, r2 = fit_rise(t, v, i1, i2)

        # Reject poor fits based on R^2
        if r2 < 0.98:
            continue

        A, tau, C = popt
        all_taus_up.append(tau)

        v_fit = exp_rise(tw, A, tau, C)
        plt.plot(t[i1:i2] * 1e6, v_fit * 1e3, lw=2)

    plt.xlabel("time (µs)")
    plt.ylabel(f"{SIGNAL_PARAM} (mV)")
    plt.title("Exponential rise fits (low plateau → high plateau)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print("\n=== τ_rise results (only good fits) ===")
    for i, tau in enumerate(all_taus_up):
        print(f"Rise {i}: tau_up = {tau * 1e6:.3f} µs")

    if all_taus_up:
        print("\nMean tau_up:", np.mean(all_taus_up) * 1e6, "µs")
        print("Std  tau_up:", np.std(all_taus_up) * 1e6, "µs")
    else:
        print("\nNo good fits (check thresholds or R²).")


if __name__ == "__main__":
    main()



