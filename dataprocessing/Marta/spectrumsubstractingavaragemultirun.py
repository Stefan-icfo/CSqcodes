#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# ===================== USER CONFIG =====================
DATA_RUNS       = [3942, 3946, 3948]   # i tre run "data"
BACKGROUND_RUNS = [3952]               # uno (o più) run di background
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"
out_dir  = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
os.makedirs(out_dir, exist_ok=True)

# Subtraction mode:
DO_VOLTAGE_DOMAIN = True   # True: V_sig = sqrt(P_data) - sqrt(P_bg);  P_sig = V_sig^2
USE_INDEX_AXIS    = False  # True: x = indice; False: x = frequenza del primo data-run (troncata)

# Averaging (smoothing per il plot):
AVG_BLOCK = 5    # media a blocchi NON sovrapposti di 5 punti
# =======================================================

qc.config["core"]["db_location"] = db_path

# --------------- helpers ---------------
CANDIDATE_SIGNAL_KEYS = [
    "avg_avg_psd_nodrive", "avg_psd_nodrive",
    "avg_avg_psd", "avg_psd",
    "V_fft_avg_avg", "PSD", "psd"
]
FREQ_KEY_CANDIDATES = ["freq_param", "frequency", "f", "freq"]

def load_1d_in_order(run_id):
    """
    Carica una traccia 1D nell'ordine di acquisizione (senza sorting).
    Ritorna (x_MHz, y, dep_key, freq_key).
    """
    ds = qc.load_by_id(run_id)
    pd = ds.get_parameter_data()

    # dependent
    dep = next((k for k in CANDIDATE_SIGNAL_KEYS if k in pd), next(iter(pd.keys())))
    blk = pd[dep]

    # frequency key
    freq_key = None
    for k in blk.keys():
        if k == dep:
            continue
        if any(t in k.lower() for t in FREQ_KEY_CANDIDATES):
            freq_key = k
            break
    if freq_key is None:
        freq_key = next((k for k in blk.keys() if k != dep), None)

    # flatten (ordine originale)
    x = np.asarray(blk[freq_key]).ravel() / 1e6 if freq_key else np.arange(len(np.asarray(blk[dep]).ravel()), dtype=float)
    y = np.asarray(blk[dep]).ravel()

    # pulizia NaN
    m = np.isfinite(y)
    if freq_key is not None:
        m &= np.isfinite(x)
        x = x[m]
    else:
        x = np.arange(np.count_nonzero(m), dtype=float)
    y = y[m]

    print(f"[run {run_id}] dep='{dep}', f_key='{freq_key}', N={len(y)}")
    return x, y, dep, freq_key

def block_average(x, y, n):
    """
    Media a blocchi non sovrapposti di ampiezza n.
    I punti in coda che non formano un blocco completo vengono scartati.
    """
    if n <= 1:
        return x.copy(), y.copy()
    m = (len(y) // n) * n
    if m == 0:
        return np.array([]), np.array([])
    xb = x[:m].reshape(-1, n).mean(axis=1)
    yb = y[:m].reshape(-1, n).mean(axis=1)
    return xb, yb

# --------------- carica tutti i run ---------------
data_traces = [load_1d_in_order(r) for r in DATA_RUNS]
bg_traces   = [load_1d_in_order(r) for r in BACKGROUND_RUNS]

# determina la lunghezza minima comune per tagliare punto-per-punto
all_lengths = [len(y) for (_, y, _, _) in data_traces + bg_traces]
n = min(all_lengths) if all_lengths else 0
if n == 0:
    raise RuntimeError("Nessun punto da elaborare (una o più tracce sono vuote).")

# asse x per il plot: indice o frequenza del primo data-run
x_ref = data_traces[0][0][:n]
x_plot = (np.arange(n, dtype=float) if USE_INDEX_AXIS else x_ref)

# estrai e taglia gli array y
data_Ys = [y[:n] for (_, y, _, _) in data_traces]
bg_Ys   = [y[:n] for (_, y, _, _) in bg_traces]

# --------------- media fra i run e sottrazione ---------------
if DO_VOLTAGE_DOMAIN:
    # media in dominio di VOLTAGGIO: V_mean = mean(sqrt(P_i))
    V_datas = [np.sqrt(np.clip(y, 0.0, None)) for y in data_Ys]
    V_bg    = [np.sqrt(np.clip(y, 0.0, None)) for y in bg_Ys]

    V_data_mean = np.mean(np.stack(V_datas, axis=0), axis=0)
    V_bg_mean   = np.mean(np.stack(V_bg, axis=0),    axis=0)

    V_sig = V_data_mean - V_bg_mean
    y_sig = V_sig**2
    y_label = "PSD (pW/Hz)"
else:
    # media in dominio di POTENZA: P_sig = mean(P_data) - mean(P_bg)
    P_data_mean = np.mean(np.stack(data_Ys, axis=0), axis=0)
    P_bg_mean   = np.mean(np.stack(bg_Ys,   axis=0), axis=0)
    y_sig = P_data_mean - P_bg_mean
    y_label = "PSD (W/Hz)  [power subtraction]"

# diagnostica
print(f"[process] usando i primi {n} punti (allineamento per indice).")
d_min, d_max = np.nanmin(np.stack(data_Ys)), np.nanmax(np.stack(data_Ys))
b_min, b_max = np.nanmin(np.stack(bg_Ys)),   np.nanmax(np.stack(bg_Ys))
r_min, r_max = np.nanmin(y_sig), np.nanmax(y_sig)
print(f"  data range: [{d_min:.3e}, {d_max:.3e}]")
print(f"  back range: [{b_min:.3e}, {b_max:.3e}]")
print(f"  result rng: [{r_min:.3e}, {r_max:.3e}]")

# converte in pW/Hz per il plot (tuo convenzione)
y_sig_plot = y_sig / 1e-12 if DO_VOLTAGE_DOMAIN else y_sig  # se power-subtraction, lasciare in W/Hz

# --------------- media a blocchi di 5 punti per il plot ---------------
x_avg, y_avg = block_average(x_plot, y_sig_plot, AVG_BLOCK)

# --------------- plot ---------------
plt.figure(figsize=(10, 6))

# curva media (smussata per il plot)
if x_avg.size and y_avg.size:
    plt.plot(x_avg, y_avg, '-o', ms=3.5, lw=1.8,
             label=f"{AVG_BLOCK}-point block average")

xlabel = "Point index" if USE_INDEX_AXIS else "Frequency (MHz)"
plt.xlabel(xlabel, fontsize=18, fontname="Calibri")
plt.ylabel(y_label, fontsize=18, fontname="Calibri")

title_runs = "_".join(str(r) for r in DATA_RUNS)
title_bg   = "_".join(str(r) for r in BACKGROUND_RUNS)
plt.title(f"Mean of runs [{title_runs}] − background [{title_bg}]  (block {AVG_BLOCK})",
          fontsize=18, fontname="Calibri")

plt.tick_params(axis='both', which='major', labelsize=16)
for lab in (plt.gca().get_xticklabels() + plt.gca().get_yticklabels()):
    lab.set_fontname("Calibri")

out_png = os.path.join(out_dir, f"psd_mean_runs_{title_runs}_minus_bg_{title_bg}_avg{AVG_BLOCK}.png")
plt.savefig(out_png, dpi=200)
print("Saved figure to:", out_png)

try:
    plt.show()
except Exception as e:
    print("plt.show() failed:", e)


