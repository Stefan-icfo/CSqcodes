import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ====== DATABASE PATH AND RUN ID ======
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v44_05_12_25.db"
run_id = 20
# Use demod X as signal instead of v_r
SIGNAL_PARAM = "x"          # could also be "v_r", "v_r_avg", ...

# ---------- Exponential decay model ----------
def exp_decay(t, A, tau, C):
    """Exponential decay: C + A * exp(-t/tau)."""
    return C + A * np.exp(-t / tau)

# ---------- Load time trace from QCoDeS database ----------
def load_time_trace_from_db(db_path, run_id, signal_param=SIGNAL_PARAM):
    """
    Load time and signal arrays from a QCoDeS run.

    Returns
    -------
    t : np.ndarray
        Time array (s).
    v : np.ndarray
        Signal array.
    """
    initialise_or_create_database_at(db_path)
    ds = load_by_id(run_id)
    print(f"Loaded run_id = {run_id} from DB:")
    print(ds)

    data_dict = ds.get_parameter_data(signal_param)
    inner = data_dict[signal_param]

    # time parameter = any key that is not the signal itself
    time_keys = [k for k in inner.keys() if k != signal_param]
    if len(time_keys) == 0:
        raise RuntimeError("Could not find time parameter in dataset.")
    time_name = time_keys[0]

    t = np.asarray(inner[time_name])
    v = np.asarray(inner[signal_param])

    print(f"Time parameter:   {time_name}")
    print(f"Signal parameter: {signal_param}")
    print(f"Number of points: {len(t)}")

    return t, v

# ---------- Find decay segments ----------
def extract_decay_segments(t, v,
                           peak_height_frac=0.57,
                           start_fraction=0.999,
                           min_points=10):
    """
    For each peak:
      - estimate a global low baseline from the 5th percentile
      - find peaks above a threshold between baseline and max
      - for each peak, define:
          * end of fit = midpoint (in index) between this peak and next peak
          * start of fit = first index after the peak where
              v <= baseline + start_fraction*(peak - baseline)
        (so the fit starts slightly below the peak, along the decay)
    Returns
    -------
    decay_segments : list of (i_start, i_end) with i_end exclusive.
    """

    # global low baseline (low plateau)
    baseline = np.percentile(v, 5)
    v_max = np.max(v)

    # threshold for finding peaks
    height_threshold = baseline + peak_height_frac * (v_max - baseline)

    # minimal distance between peaks (~30 µs)
    if len(t) > 1:
        dt = np.median(np.diff(t))
        distance_samples = max(1, int(30e-6 / dt))
    else:
        distance_samples = 20

    peaks, props = find_peaks(v, height=height_threshold,
                              distance=distance_samples)
    print(f"Found {len(peaks)} peaks.")

    decay_segments = []

    for i, pk in enumerate(peaks):
        # next peak index (for last one, extrapolate using previous period)
        if i < len(peaks) - 1:
            next_pk = peaks[i + 1]
        else:
            if len(peaks) > 1:
                period = peaks[i] - peaks[i - 1]
            else:
                period = distance_samples
            next_pk = min(len(v) - 1, pk + period)

        # end of fit = midpoint between two peaks
        mid_idx = (pk + next_pk) // 2

        # amplitude relative to baseline
        amp = v[pk] - baseline
        if amp <= 0:
            continue

        # level where we start the decay fit
        start_level = baseline + start_fraction * amp

        # search for first point after peak that drops below start_level
        i_start = None
        for j in range(pk, mid_idx):
            if v[j] <= start_level:
                i_start = j
                break
        if i_start is None:
            continue

        i_end = mid_idx

        if i_end - i_start < min_points:
            continue

        decay_segments.append((i_start, i_end))

    print(f"Decay segments found: {len(decay_segments)}")
    return decay_segments

# ---------- Fit one decay ----------
def fit_decay(t, v, i1, i2):
    """
    Fit exponential decay between indices [i1, i2).
    """
    t_win = t[i1:i2] - t[i1]
    v_win = v[i1:i2]

    # initial guesses
    A0 = v_win[0] - v_win[-1]
    tau0 = (t_win[-1] - t_win[0]) / 2
    C0 = v_win[-1]

    popt, _ = curve_fit(exp_decay, t_win, v_win,
                        p0=[A0, tau0, C0], maxfev=5000)
    return popt, t_win, v_win

# ---------- MAIN ----------
def main():
    # load t and x from DB
    t, v = load_time_trace_from_db(db_path, run_id, SIGNAL_PARAM)

    print("Data loaded from database.")

    decay_windows = extract_decay_segments(
        t, v,
        peak_height_frac=0.6,
        start_fraction=0.98,   # move closer/further from the peak if needed
        min_points=10
    )

    all_taus = []

    plt.figure(figsize=(10, 5))
    # x is in V → plot in µV
    plt.plot(t * 1e6, v * 1e6, 'k-', label="data")

    for (i1, i2) in decay_windows:
        popt, tw, vw = fit_decay(t, v, i1, i2)
        A, tau, C = popt
        all_taus.append(tau)

        v_fit = exp_decay(tw, A, tau, C)
        plt.plot(t[i1:i2] * 1e6, v_fit * 1e6, lw=2)

    plt.xlabel("time (µs)")
    plt.ylabel("x (µV)")
    plt.title("Exponential decay fits on X (start near peak, end mid-plateau)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print("\n=== τ_decay results on X ===")
    for i, tau in enumerate(all_taus):
        print(f"Decay {i}: tau = {tau*1e6:.3f} µs")

    if all_taus:
        print("\nMean:", np.mean(all_taus) * 1e6, "µs")
        print("Std :",  np.std(all_taus) * 1e6, "µs")
    else:
        print("\nNo decay segments found – check peak_height_frac/start_fraction.")

if __name__ == "__main__":
    main()







