import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ====== DATABASE PATH AND RUN ID ======
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v35_26_11_25.db"
run_id = 37

# Use the demod X signal stored in the dataset
SIGNAL_PARAM = "x"


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
        Name of the signal parameter ('x', 'y', 'v_r', 'v_r_avg', ...).

    Returns
    -------
    t : np.ndarray
        Time array (seconds).
    v : np.ndarray
        Signal array.
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
def find_all_peaks(t, v):
    """
    Detect peaks corresponding to the high plateaus.

    Strategy:
      - Use a high percentile (e.g. 80th) as an amplitude threshold.
      - Optionally require some prominence so we don't pick pure noise.
    """
    if len(t) > 1:
        dt = np.median(np.diff(t))
        distance_samples = max(1, int(30e-6 / dt))  # ~30 µs
    else:
        distance_samples = 50

    # High level threshold: upper 20% of the signal distribution
    height_threshold = np.percentile(v, 80)
    amp = np.max(v) - np.min(v)
    prominence = 0.1 * amp  # require at least ~10% of full swing

    peaks, _ = find_peaks(
        v,
        height=height_threshold,
        prominence=prominence,
        distance=distance_samples,
    )
    print(f"Found {len(peaks)} peaks (high plateaus).")
    return peaks


# ---------- Extract rise segments (low plateau → middle of high plateau) ----------
def extract_rise_segments(t, v,
                          min_points=20,
                          valley_margin_samples=3):
    """
    For each peak, define a rise window that starts near the bottom
    of the low plateau and ends in the middle of the high plateau.

    Geometry-based:
      - peaks: indices of high plateaus.
      - For each peak pk:
          * prev_pk = previous peak (or extrapolated)
          * next_pk = next peak (or extrapolated)
          * win_start = midpoint between prev_pk and pk
          * win_end   = midpoint between pk and next_pk
            -> window covers: low plateau → rise → high plateau.
          * In window, search for the minimum BEFORE pk: valley index.
            Start index = valley index - valley_margin_samples.
          * End index = midpoint between pk and win_end
            (approx middle of high plateau in time).

    Returns
    -------
    rise_segments : list of (i_start, i_end)
        Index intervals [i_start, i_end) for the rise fits.
    """
    peaks = find_all_peaks(t, v)
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

        # Window that covers low plateau → rise → high plateau
        win_start = (prev_pk + pk) // 2
        win_end = (pk + next_pk) // 2

        if win_end - win_start < min_points + 4:
            continue

        v_win = v[win_start:win_end + 1]
        idx_pk_rel = pk - win_start

        # --- valley in low plateau BEFORE the peak ---
        v_low = v_win[:idx_pk_rel]
        if len(v_low) < 5:
            continue

        valley_rel = int(np.argmin(v_low))
        valley_rel = max(0, valley_rel - valley_margin_samples)

        # --- middle of high plateau in time ---
        plateau_mid_abs = (pk + win_end) // 2
        plateau_mid_rel = plateau_mid_abs - win_start

        if plateau_mid_rel - valley_rel < min_points:
            continue

        i_start = win_start + valley_rel
        i_end = win_start + plateau_mid_rel  # exclusive in fit

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

    C0 = v_win[0]                     # offset ~ low plateau
    A0 = v_win[-1] - v_win[0]         # amplitude guess
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

    # Extract windows for exponential rise fits
    rise_windows = extract_rise_segments(
        t, v,
        min_points=30,          # ensure we have a long window
        valley_margin_samples=5 # start a bit BEFORE valley → more low plateau
    )

    all_taus_up = []

    plt.figure(figsize=(10, 5))
    # x is stored in volts → plot in µV
    plt.plot(t * 1e6, v * 1e6, 'k-', label="data")

    for (i1, i2) in rise_windows:
        popt, tw, vw, r2 = fit_rise(t, v, i1, i2)

        # Reject poor fits based on R^2
        if r2 < 0.98:
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
