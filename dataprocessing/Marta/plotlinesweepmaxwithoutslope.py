# -*- coding: utf-8 -*-
# CUCITURA DI DUE RUN IN UN'UNICA LINEA NERA CONTINUA
# - ricostruisce G(x,y) via pivot (robusto all'ordine inner/outer)
# - estrae la cresta dei massimi per riga
# - concatena le due run su un unico asse di sweep e plottà una sola linea nera
# - SCAMBIA GLI ASSI: ch02 → X, ch06 → Y
# - STIMA LA PENDENZA, POI RUOTA I DATI (default: appiattisce la slope)

import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
import numpy as np

# ---------- Helper: massimi riga-per-riga (robusto ai NaN) ----------
def maxima_trace(x_vec, y_vec, G2d):
    if G2d.ndim != 2:
        raise ValueError("G deve essere 2D")
    valid_rows = np.isfinite(G2d).any(axis=1)
    y_valid = y_vec[valid_rows]
    G_valid = G2d[valid_rows, :]
    G_safe = np.where(np.isfinite(G_valid), G_valid, -np.inf)
    max_idx = np.argmax(G_safe, axis=1)
    x_at_max = x_vec[max_idx]
    gmax = G_safe[np.arange(G_safe.shape[0]), max_idx]
    return x_at_max, y_valid, gmax

# ---------- Helper: costruzione griglia robusta dal dataset ----------
def grid_from_dataset(ds, preferred_param='G'):
    groups = ds.get_parameter_data()
    if preferred_param in groups:
        pname = preferred_param
    else:
        if not groups:
            raise ValueError("Nessun parametro disponibile nel dataset.")
        pname = next(iter(groups.keys()))
    grp = groups[pname]

    setpoints = [k for k in grp.keys() if k != pname]
    if len(setpoints) < 1:
        raise ValueError("Attesi >=1 setpoint; trovato <1.")
    if len(setpoints) == 1:
        x_name = y_name = setpoints[0]
    else:
        sp1, sp2 = setpoints[:2]
        v1 = np.asarray(grp[sp1]); v2 = np.asarray(grp[sp2])
        if len(np.unique(v1)) >= len(np.unique(v2)):
            x_name, y_name = sp1, sp2
        else:
            x_name, y_name = sp2, sp1

    df = pd.DataFrame({'x': np.asarray(grp[x_name]),
                       'y': np.asarray(grp[y_name]),
                       'val': np.asarray(grp[pname])})

    piv = df.pivot_table(index='y', columns='x', values='val', aggfunc='mean', sort=True)
    x = piv.columns.values
    y = piv.index.values
    G = piv.to_numpy()

    # y crescente verso l'alto
    if len(y) > 1 and y[0] > y[-1]:
        y = y[::-1]
        G = np.flipud(G)
    return pname, x_name, y_name, x, y, G

# ---------- Helper: stitching di due serie (y, x_at_max) ----------
def stitch_ridges(y_list, x_list, tol=1e-9):
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

# ---------- Imposta il database ----------
qc.config["core"]["db_location"] = (
    "C:\\" + "Users\\" + "LAB-nanooptomechanic\\" + "Documents\\" +
    "MartaStefan\\" + "CSqcodes\\" + "Data\\" + "Raw_data\\" + "CD12_B5_F4v19_211025.db"
)
_ = qc.experiments()

# ---------- IDs dei dataset ----------
ID1 = 35
ID2 = None  # cambia o metti None se hai solo una run

# ---------- Caricamento ----------
ds1 = qc.load_by_id(ID1)
ds2 = None
if ID2 is not None:
    try:
        ds2 = qc.load_by_id(ID2)
    except Exception as e:
        print(f"[Info] Non carico DS {ID2}: {e}")

# ---------- Ricostruzione + cresta ----------
p1, x1_name, y1_name, x1, y1, G1 = grid_from_dataset(ds1, preferred_param='G')
x1_max, y1_valid, g1_max = maxima_trace(x1, y1, G1)
print(f"DS {ID1} → '{p1}', x='{x1_name}', y='{y1_name}', range y: [{y1_valid[0]}, {y1_valid[-1]}]")

ys_list = [y1_valid]
xs_list = [x1_max]

if ds2 is not None:
    p2, x2_name, y2_name, x2, y2, G2 = grid_from_dataset(ds2, preferred_param='G')
    x2_max, y2_valid, g2_max = maxima_trace(x2, y2, G2)
    print(f"DS {ID2} → '{p2}', x='{x2_name}', y='{y2_name}', range y: [{y2_valid[0]}, {y2_valid[-1]}]")
    if y2_name != y1_name:
        print(f"[Attenzione] L'asse di sweep (y) differisce tra run: '{y1_name}' vs '{y2_name}'. Procedo comunque.")
    if x2_name != x1_name:
        print(f"[Attenzione] L'asse trasverso (x) differisce tra run: '{x1_name}' vs '{x2_name}'. Procedo comunque.")
    ys_list.append(y2_valid)
    xs_list.append(x2_max)

# ---------- Stitching ----------
# ---------- STITCH ----------
y_st, x_st = stitch_ridges(ys_list, xs_list, tol=1e-9)

REQ_LO, REQ_HI = 0.79, 0.86

def has(s, key): return key.lower() in str(s).lower()

# Metto esplicitamente: X = ch02, Y = ch06 (grezzi)
if has(y1_name, 'ch06'):
    ch06 = y_st.astype(float)
    ch02 = x_st.astype(float)
elif has(x1_name, 'ch06'):
    ch06 = x_st.astype(float)
    ch02 = y_st.astype(float)
else:
    raise ValueError(f"Non trovo 'ch06' nei nomi assi: x='{x1_name}', y='{y1_name}'")

# Pulisci NaN/inf
mfin = np.isfinite(ch02) & np.isfinite(ch06)
ch02, ch06 = ch02[mfin], ch06[mfin]

# --- Stima della pendenza m usando (ch06=0.79, 0.86) -> interpolo ch02(ch06) ---
# ordina per ch06 per l'interpolazione inversa X(Y)
ord_y = np.argsort(ch06)
yy, xx = ch06[ord_y], ch02[ord_y]

# finestra effettiva disponibile (clamp se fuori range)
LO = max(REQ_LO, float(yy.min()))
HI = min(REQ_HI, float(yy.max()))
if LO >= HI:
    raise ValueError(f"Nessuna sovrapposizione con [{REQ_LO},{REQ_HI}] V su ch06. Range disponibile: [{yy.min():.5g},{yy.max():.5g}] V")

x_at_lo = np.interp(LO, yy, xx)
x_at_hi = np.interp(HI, yy, xx)

# m = d(ch06)/d(ch02) stimato con i due estremi richiesti (o clampati)
m = (HI - LO) / (x_at_hi - x_at_lo)
x0 = x_at_lo  # pivot così la cresta si appiattisce a LO

# --- Compensazione di cross-capacitance (shear) ---
ch06_comp = ch06 - m * (ch02 - x0)

# Filtro per visualizzare solo la finestra 0.79–0.86 della ch06 compensata
mask = (ch06_comp >= REQ_LO) & (ch06_comp <= REQ_HI)
if not np.any(mask):
    # se non c'è copertura completa, mostra la parte sovrapposta reale
    mask = (ch06_comp >= LO) & (ch06_comp <= HI)

X = ch02[mask]
Y = ch06_comp[mask]

# ordina per X per una linea pulita
ordx = np.argsort(X)
X, Y = X[ordx], Y[ordx]

# ---------- PLOT: X = ch02, Y = ch06 compensata ----------
plt.figure(figsize=(10, 4))
plt.plot(X, Y, 'k-', linewidth=1.8)
plt.xlabel('ch02')
plt.ylabel('ch06 (compensata) [V]')
plt.yticks([0.79, 0.86])
plt.ylim(REQ_LO, REQ_HI)
plt.title(f"Cross-cap compensata (m={m:.4g}, pivot x0={x0:.4g})")
plt.tight_layout()
plt.show()

print(f"Compensazione fatta: m={m:.6g} da (0.79V -> x={x_at_lo:.6g}) e (0.86V -> x={x_at_hi:.6g}).")

