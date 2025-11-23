import os
import numpy as np
import matplotlib.pyplot as plt
from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

# ==================== USER SETTINGS ====================

db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v32_19_11_25.db"
run_id = 592           # run che contiene la demod time trace
which_pulse = 1        # quale impulso usare per il fit (0=primo, 1=secondo, ...)

# =======================================================


def exp_rise(t, A, tau, C):
    """C + A * (1 - exp(-t/tau))"""
    return C + A * (1 - np.exp(-t / tau))


def exp_decay(t, A, tau, C):
    """C + A * exp(-t/tau)"""
    return C + A * np.exp(-t / tau)


def fit_pulse_rise_decay(t, v_r, which_pulse=1, plot=True):
    """
    Fit esponenziale su salita (uptime) e discesa (decay time) di un impulso.

    Parameters
    ----------
    t : array (s)
    v_r : array (V)
    which_pulse : int
        Indice del picco da usare (0 = primo picco, 1 = secondo, ecc.).
    plot : bool
        Se True, mostra i plot con fit.

    Returns
    -------
    results : dict con tau_rise, tau_decay, ecc.
    """
    t = np.asarray(t)
    v_r = np.asarray(v_r)

    # Trova i picchi (impulsi) – soglia circa sopra il livello medio
    thr = np.mean(v_r) + 0.3 * (np.max(v_r) - np.mean(v_r))
    peaks, props = find_peaks(v_r, height=thr)

    if len(peaks) == 0:
        raise RuntimeError("Non trovo nessun picco: controlla i dati o la soglia.")

    if which_pulse >= len(peaks):
        raise ValueError(f"which_pulse={which_pulse}, ma ci sono solo {len(peaks)} picchi.")

    pk = peaks[which_pulse]

    # Minimi prima e dopo il picco (per delimitare finestra di fit)
    mins, _ = find_peaks(-v_r)
    mins_before = mins[mins < pk]
    mins_after = mins[mins > pk]

    if len(mins_before) == 0 or len(mins_after) == 0:
        raise RuntimeError("Non trovo minimo prima/dopo il picco. Prova un altro which_pulse.")

    min_before = mins_before[-1]   # ultimo minimo prima del picco
    min_after = mins_after[0]      # primo minimo dopo il picco

    # -------- FIT RISE (uptime) --------
    t_rise = t[min_before:pk + 1]
    v_rise = v_r[min_before:pk + 1]
    t_rise_rel = t_rise - t_rise[0]

    C0_r = v_rise[0]
    A0_r = v_rise[-1] - v_rise[0]
    tau0_r = (t_rise_rel[-1] - t_rise_rel[0]) / 5.0

    popt_rise, pcov_rise = curve_fit(
        exp_rise, t_rise_rel, v_rise,
        p0=[A0_r, tau0_r, C0_r]
    )
    A_r, tau_r, C_r = popt_rise

    # -------- FIT DECAY --------
    t_dec = t[pk:min_after + 1]
    v_dec = v_r[pk:min_after + 1]
    t_dec_rel = t_dec - t_dec[0]

    C0_d = v_dec[-1]
    A0_d = v_dec[0] - v_dec[-1]
    tau0_d = (t_dec_rel[-1] - t_dec_rel[0]) / 5.0

    popt_dec, pcov_dec = curve_fit(
        exp_decay, t_dec_rel, v_dec,
        p0=[A0_d, tau0_d, C0_d]
    )
    A_d, tau_d, C_d = popt_dec

    if plot:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(t * 1e6, v_r * 1e3, 'o', ms=2, label='dati')

        ax.plot(t_rise * 1e6,
                exp_rise(t_rise_rel, *popt_rise) * 1e3,
                '-', label=f'rise fit, τ_rise = {tau_r*1e6:.1f} µs')

        ax.plot(t_dec * 1e6,
                exp_decay(t_dec_rel, *popt_dec) * 1e3,
                '-', label=f'decay fit, τ_decay = {tau_d*1e6:.1f} µs')

        ax.set_xlabel('time (µs)')
        ax.set_ylabel('v_r (mV)')
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        plt.show()

    return {
        "tau_rise_s": tau_r,
        "tau_rise_us": tau_r * 1e6,
        "tau_decay_s": tau_d,
        "tau_decay_us": tau_d * 1e6,
        "rise_params": popt_rise,
        "decay_params": popt_dec,
        "indices": {
            "min_before": int(min_before),
            "peak": int(pk),
            "min_after": int(min_after),
        }
    }


def main():
    # inizializza DB e carica dataset
    initialise_or_create_database_at(db_path)
    ds = load_by_id(run_id)

    df = ds.to_pandas_dataframe()

    # estrai colonne (adatta i nomi se nel tuo dataset sono diversi)
    t = df['time_param'].to_numpy()   # in secondi
    v_r = df['v_r'].to_numpy()        # in Volt

    # per sicurezza, ordiniamo per tempo crescente
    idx_sort = np.argsort(t)
    t = t[idx_sort]
    v_r = v_r[idx_sort]

    # Fit di un impulso
    res = fit_pulse_rise_decay(t, v_r, which_pulse=which_pulse, plot=True)

    print(f"=== Risultati fit run {run_id} (impulso {which_pulse}) ===")
    print(f"tau_rise  = {res['tau_rise_us']:.2f} µs")
    print(f"tau_decay = {res['tau_decay_us']:.2f} µs")


if __name__ == "__main__":
    main()





