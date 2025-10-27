# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import re
from typing import Dict, Tuple

# =========================================================
# CONFIG (tuoi valori)
# =========================================================
qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"
)
background_id = 2885
all_data_ids = list(range(2948, 3033, 2))  # 2948..3032

# --- Analisi/plot ---
BASELINE_PERCENTILE = 5.0    # baseline robusta dopo sottrazione background
CLIP_NEGATIVES = True        # clip <0 per il plot/fill
SAVE_FIGS = False            # salva png locali
PLOT_PER_RUN = True          # mostra anche un grafico per ogni run

# =========================================================
# FASE 2: INSERISCI QUI le finestre (in MHz) dopo aver visto i plot
# esempio: {2948: {"left": (136.760,136.767), "right": (136.786,136.797)}, ...}
# lascia vuoto al primo giro
# =========================================================
MANUAL_WINDOWS: Dict[int, Dict[str, Tuple[float,float]]] = {
    # 2948: {"left": (136.760, 136.767), "right": (136.786, 136.797)},
}

# =========================================================
# HELPERS
# =========================================================
def find_psd_key(param_data: dict) -> str:
    for k, v in param_data.items():
        if "freq_param" in v:
            return k
    return list(param_data.keys())[0]

def load_psd(run_id):
    ds = qc.load_by_id(run_id)
    p = ds.get_parameter_data()
    key = find_psd_key(p)
    y = np.asarray(p[key][key]).flatten()
    x = np.asarray(p[key]["freq_param"]).flatten() / 1e6  # MHz
    idx = np.argsort(x)
    return x[idx], y[idx], ds

def extract_drive_amplitude(exp_name: str):
    pats = [
        r"[dD]rive\s*([0-9]+)\s*u[V]?",
        r"_(\d+)u_(?=1Dcs)",
        r"_(\d+)u(?:_|$)",
        r"\b(\d+)u\b",
    ]
    for pat in pats:
        m = re.search(pat, exp_name)
        if m:
            return int(m.group(1))
    return None

def integrate_between_freqs(x, y_corr, fmin, fmax, baseline_mode="local"):
    """Integra (y - baseline) tra fmin e fmax (MHz). baseline locale = percentile basso nella finestra."""
    i0 = int(np.searchsorted(x, fmin, side="left"))
    i1 = int(np.searchsorted(x, fmax, side="right")) - 1
    i0 = max(0, min(i0, len(x)-2))
    i1 = max(i0+1, min(i1, len(x)-1))

    ywin = y_corr[i0:i1+1].copy()
    if baseline_mode == "local":
        baseline = np.percentile(ywin, BASELINE_PERCENTILE)
    else:
        baseline = np.percentile(y_corr, BASELINE_PERCENTILE)

    yb = ywin - baseline
    if CLIP_NEGATIVES:
        yb = np.clip(yb, 0, None)
    area = np.trapz(yb, x[i0:i1+1])  # (W/Hz) * MHz
    return area, i0, i1, float(baseline)

# =========================================================
# CARICA BACKGROUND
# =========================================================
x_bg, y_bg, _ = load_psd(background_id)

# =========================================================
# FASE 1: CARICA TUTTI I RUN, SOTTRAI BACKGROUND E PLOTTA
# =========================================================
runs = []  # lista di dict: {"run_id", "drive", "x", "y_corr", "exp_name"}
for run_id in all_data_ids:
    try:
        x, y, ds = load_psd(run_id)
        n = min(len(y), len(y_bg))
        x = x[:n]; y = y[:n]; y_bg_trim = y_bg[:n]
        y_corr = y - y_bg_trim
        drive = extract_drive_amplitude(ds.exp_name)
        runs.append({"run_id": run_id, "drive": drive, "x": x, "y_corr": y_corr, "exp_name": ds.exp_name})
    except Exception as e:
        print(f"⚠️ Run {run_id}: skipped → {e}")

if not runs:
    raise RuntimeError("❌ Nessuno spettro caricato.")

# overlay di tutti i run
plt.figure(figsize=(11.5, 5))
for r in runs:
    yplot = np.clip(r["y_corr"], 0, None) if CLIP_NEGATIVES else r["y_corr"]
    lbl = f"{r['run_id']}" + (f" ({r['drive']} µV)" if r['drive'] is not None else "")
    plt.plot(r["x"], yplot, lw=0.9, alpha=0.8, label=lbl)
plt.xlabel("Frequenza (MHz)")
plt.ylabel("PSD (W/Hz)")
plt.title("Tutti i run (background sottratto)")
plt.legend(ncol=4, fontsize=7, framealpha=0.5)
plt.tight_layout()
if SAVE_FIGS: plt.savefig("all_runs_overlay.png", dpi=180)
plt.show()

# tabella riassunto utile per scegliere le finestre
print("\n=== RIEPILOGO RUN per scegliere finestre (MHz) ===")
for r in runs:
    x = r["x"]
    print(f"Run {r['run_id']:4d}  drive={str(r['drive'])+' µV':>6}  freq_range=[{x[0]:.6f}, {x[-1]:.6f}]  Npts={len(x)}")

# plot per-run
if PLOT_PER_RUN:
    for r in runs:
        x = r["x"]; y_corr = r["y_corr"]
        yplot = np.clip(y_corr, 0, None) if CLIP_NEGATIVES else y_corr
        baseline_glob = np.percentile(y_corr, BASELINE_PERCENTILE)
        plt.figure(figsize=(11, 4.2))
        plt.plot(x, yplot, lw=1.0, label=f"Run {r['run_id']} (drive={r['drive']} µV)")
        plt.axhline(max(baseline_glob,0) if CLIP_NEGATIVES else baseline_glob,
                    ls="--", lw=0.8, color="gray", alpha=0.8, label="baseline (globale)")
        plt.xlabel("Frequenza (MHz)")
        plt.ylabel("PSD (W/Hz)")
        plt.title("Spettro background-sottratto")
        plt.legend(loc="upper left", fontsize=8)
        plt.tight_layout()
        if SAVE_FIGS: plt.savefig(f"run_{r['run_id']}.png", dpi=180)
        plt.show()

# =========================================================
# FASE 2: Se hai compilato MANUAL_WINDOWS, calcola aree e colora verde
# =========================================================
if MANUAL_WINDOWS:
    print("\n=== CALCOLO AREE (baseline locale in ciascuna finestra) ===")
    summary = []
    for r in runs:
        rid = r["run_id"]
        if rid not in MANUAL_WINDOWS:
            continue
        x = r["x"]; y_corr = r["y_corr"]
        spec = MANUAL_WINDOWS[rid]

        areaL, l0, l1, baseL = integrate_between_freqs(x, y_corr, *spec["left"], baseline_mode="local")
        areaR, r0, r1, baseR = integrate_between_freqs(x, y_corr, *spec["right"], baseline_mode="local")

        # plot con riempimento verde
        yplot = np.clip(y_corr, 0, None) if CLIP_NEGATIVES else y_corr
        plt.figure(figsize=(11, 4.2))
        lbl = f"Run {rid} (drive={r['drive']} µV)"
        plt.plot(x, yplot, lw=1.0, label=lbl)
        # baseline locali (tratteggiate sottili)
        plt.hlines([max(baseL,0) if CLIP_NEGATIVES else baseL,
                    max(baseR,0) if CLIP_NEGATIVES else baseR],
                   xmin=[x[l0], x[r0]], xmax=[x[l1], x[r1]],
                   colors=["gray","gray"], linestyles="--", linewidths=0.8, alpha=0.8)
        # fill verde
        plt.fill_between(x[l0:l1+1],
                         max(baseL,0) if CLIP_NEGATIVES else baseL,
                         yplot[l0:l1+1], alpha=0.35, color="green",
                         label=f"Area L = {areaL:.3e}")
        plt.fill_between(x[r0:r1+1],
                         max(baseR,0) if CLIP_NEGATIVES else baseR,
                         yplot[r0:r1+1], alpha=0.35, color="green",
                         label=f"Area R = {areaR:.3e}")
        plt.xlabel("Frequenza (MHz)")
        plt.ylabel("PSD (W/Hz)")
        plt.title(f"Sideband integrati (Run {rid})")
        plt.legend(loc="upper left", fontsize=8)
        plt.tight_layout()
        if SAVE_FIGS: plt.savefig(f"run_{rid}_areas.png", dpi=180)
        plt.show()

        summary.append((rid, r["drive"], areaL, areaR,
                        spec["left"][0], spec["left"][1], spec["right"][0], spec["right"][1]))

        print(f"Run {rid}:  Area_L={areaL:.4e}  Area_R={areaR:.4e}  "
              f"L=({spec['left'][0]:.6f},{spec['left'][1]:.6f})  "
              f"R=({spec['right'][0]:.6f},{spec['right'][1]:.6f})")

    # riepilogo ordinato per run_id
    if summary:
        summary.sort(key=lambda t: t[0])
        print("\nRUN  Drive[µV]    Area_L (W/Hz·MHz)    Area_R (W/Hz·MHz)     L_win[MHz]                R_win[MHz]")
        for rid, drv, aL, aR, fL0, fL1, fR0, fR1 in summary:
            print(f"{rid:4d} {str(drv):>10}   {aL:>14.6e}     {aR:>14.6e}    "
                  f"({fL0:.6f},{fL1:.6f})   ({fR0:.6f},{fR1:.6f})")
else:
    print("\nℹ️  Adesso dimmi per ogni run le finestre sideband in MHz (sinistra/destra), "
          "oppure inseriscile nel dict MANUAL_WINDOWS e rilancia lo script.")



