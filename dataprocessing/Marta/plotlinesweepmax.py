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
    """
    Concatena più serie (y, x_at_max), ordina per y e rimuove duplicati vicini entro tol.
    Ritorna: y_stitched, x_stitched
    """
    y_all = np.concatenate(y_list)
    x_all = np.concatenate(x_list)
    idx = np.argsort(y_all)
    y_all = y_all[idx]
    x_all = x_all[idx]

    # rimuovi punti con y troppo vicini (duplicati/overlap tra run)
    keep = np.ones_like(y_all, dtype=bool)
    if len(y_all) > 1:
        dy = np.diff(y_all)
        keep[1:] = np.abs(dy) > tol
    return y_all[keep], x_all[keep]

# ---------- Imposta il database ----------
import sqlite3
db = qc.config["core"]["db_location"]

con = sqlite3.connect(db)
cur = con.cursor()
cur.execute("SELECT run_id, experiment_id, name, run_timestamp FROM runs ORDER BY run_id DESC LIMIT 30;")
rows = cur.fetchall()
con.close()

print("Last 30 runs:")
for r in rows:
    print(r)


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
xlbl = y1_name
ylbl = x1_name

if ds2 is not None:
    p2, x2_name, y2_name, x2, y2, G2 = grid_from_dataset(ds2, preferred_param='G')
    x2_max, y2_valid, g2_max = maxima_trace(x2, y2, G2)
    print(f"DS {ID2} → '{p2}', x='{x2_name}', y='{y2_name}', range y: [{y2_valid[0]}, {y2_valid[-1]}]")

    # Avviso se i nomi dei setpoint non coincidono (assi “diversi”)
    if y2_name != y1_name:
        print(f"[Attenzione] L'asse di sweep (y) differisce tra run: '{y1_name}' vs '{y2_name}'. Procedo comunque.")
    if x2_name != x1_name:
        print(f"[Attenzione] L'asse trasverso (x) differisce tra run: '{x1_name}' vs '{x2_name}'. Procedo comunque.")

    ys_list.append(y2_valid)
    xs_list.append(x2_max)

# ---------- Stitching e plot UNICO in nero ----------
y_st, x_st = stitch_ridges(ys_list, xs_list, tol=1e-9)


plt.figure(figsize=(7, 6))
plt.plot(y_st, x_st, 'k-', linewidth=1.8)  # TUTTO NERO, LINEA CONTINUA
plt.xlabel(xlbl)  # asse orizzontale = setpoint di riga (y_name)
plt.ylabel(ylbl)  # asse verticale   = setpoint di colonna (x_name)
plt.title("Cresta massimi di G — run cucite (linea nera continua)")
plt.tight_layout()
plt.show()

print("Fatto: cucitura completata e plottata in nero.")
