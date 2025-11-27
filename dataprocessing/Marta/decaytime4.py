

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ====== PATH DATABASE E RUN ID ======
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v35_26_11_25.db"
run_id = 37

# Nome del parametro di segnale nel dataset
# Cambia qui se vuoi usare 'v_r_avg', 'x', ecc.
SIGNAL_PARAM = "v_r"          # oppure "v_r_avg"

# ---------- Modello esponenziale ----------
def exp_decay(t, A, tau, C):
    return C + A * np.exp(-t / tau)

# ---------- CARICAMENTO DATI DAL DATABASE ----------
def load_time_trace_from_db(db_path, run_id, signal_param=SIGNAL_PARAM):
    """
    Carica tempo e segnale da un run QCoDeS:
    - db_path: percorso al file .db
    - run_id: ID del run
    - signal_param: nome del parametro di segnale (es. 'v_r', 'v_r_avg', 'x')

    Ritorna:
        t, v  (numpy array)
    """
    # punta al database giusto
    initialise_or_create_database_at(db_path)

    # carica il dataset
    ds = load_by_id(run_id)
    print(f"Caricato run_id = {run_id} dal DB:")
    print(ds)

    # get_parameter_data restituisce un dict annidato:
    # data_dict[signal_param][signal_param] -> valori
    # data_dict[signal_param][time_param_name] -> setpoint (tempo)
    data_dict = ds.get_parameter_data(signal_param)

    inner = data_dict[signal_param]

    # trova il nome del parametro tempo (tutto tranne il segnale stesso)
    time_keys = [k for k in inner.keys() if k != signal_param]
    if len(time_keys) == 0:
        raise RuntimeError("Non riesco a trovare il parametro tempo nel dataset.")
    time_name = time_keys[0]

    t = np.asarray(inner[time_name])
    v = np.asarray(inner[signal_param])

    print(f"Param tempo: {time_name}")
    print(f"Param segnale: {signal_param}")
    print(f"Numero punti: {len(t)}")

    return t, v

# ---------- MAIN: trova tutti i segmenti di decay ----------
def extract_decay_segments(t, v,
                           peak_height_frac=0.6,
                           start_fraction=0.9,
                           min_points=10):
    """
    Per ogni picco:
      - stima la baseline bassa (plateau)
      - trova il punto dove il segnale è sceso a
            baseline + start_fraction * (peak - baseline)
        (inizio fit, più in basso lungo il decay)
      - fissa la fine del fit al punto a metà tra questo picco e il successivo
    Ritorna lista di (i_start, i_end) in indici assoluti (i_end esclusivo).
    """

    # baseline globale (plateau basso)
    baseline = np.percentile(v, 5)
    v_max = np.max(v)

    # soglia per i picchi: frazione tra baseline e massimo
    height_threshold = baseline + peak_height_frac * (v_max - baseline)

    # distanza minima tra picchi (~30 µs)
    if len(t) > 1:
        dt = np.median(np.diff(t))
        distance_samples = max(1, int(30e-6 / dt))
    else:
        distance_samples = 20

    peaks, props = find_peaks(v, height=height_threshold,
                              distance=distance_samples)
    print(f"Trovati {len(peaks)} picchi.")

    decay_segments = []

    for i, pk in enumerate(peaks):
        # indice del prossimo picco
        if i < len(peaks) - 1:
            next_pk = peaks[i + 1]
        else:
            # per l'ultimo picco, usa periodo medio
            if len(peaks) > 1:
                period = peaks[i] - peaks[i - 1]
            else:
                period = distance_samples
            next_pk = min(len(v) - 1, pk + period)

        # fine fit = metà tra i due picchi
        mid_idx = (pk + next_pk) // 2

        # ampiezza del picco rispetto alla baseline
        amp = v[pk] - baseline
        if amp <= 0:
            continue

        # livello dove vogliamo iniziare il fit (es. 97% del decay)
        start_level = baseline + start_fraction * amp

        # cerca il primo punto dopo il picco che scende sotto start_level
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

    print(f"Segmenti di decay trovati: {len(decay_segments)}")
    return decay_segments

# ---------- FIT ----------
def fit_decay(t, v, i1, i2):
    """
    Fit esponenziale tra indici [i1, i2) (i2 esclusivo).
    """
    t_win = t[i1:i2] - t[i1]
    v_win = v[i1:i2]

    # guess iniziali
    A0 = v_win[0] - v_win[-1]
    tau0 = (t_win[-1] - t_win[0]) / 2
    C0 = v_win[-1]

    popt, _ = curve_fit(exp_decay, t_win, v_win,
                        p0=[A0, tau0, C0], maxfev=5000)
    return popt, t_win, v_win

# ---------- MAIN ----------
def main():
    # carica t e v dal DB
    t, v = load_time_trace_from_db(db_path, run_id, SIGNAL_PARAM)

    print("File caricato dal database.")

    decay_windows = extract_decay_segments(
        t, v,
        peak_height_frac=0.6,
        start_fraction=0.9,   # prova 0.8 o 0.9 se vuoi spostarlo su/giù
        min_points=10
    )

    all_taus = []

    plt.figure(figsize=(10, 5))
    plt.plot(t * 1e6, v * 1e3, 'k-', label="data")

    for (i1, i2) in decay_windows:
        popt, tw, vw = fit_decay(t, v, i1, i2)
        A, tau, C = popt
        all_taus.append(tau)

        v_fit = exp_decay(tw, A, tau, C)
        plt.plot(t[i1:i2] * 1e6, v_fit * 1e3, lw=2)

    plt.xlabel("time (µs)")
    plt.ylabel(f"{SIGNAL_PARAM} (mV)")
    plt.title("Exponential decay fits (start dal 97% del picco, fine a metà plateau)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print("\n=== Risultati τ_decay ===")
    for i, tau in enumerate(all_taus):
        print(f"Decay {i}: tau = {tau*1e6:.3f} µs")

    print("\nMedia:", np.mean(all_taus) * 1e6, "µs")
    print("Std  :", np.std(all_taus) * 1e6, "µs")


if __name__ == "__main__":
    main()





