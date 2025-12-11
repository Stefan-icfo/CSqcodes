import os
import numpy as np
import matplotlib.pyplot as plt

from qcodes.dataset import load_by_id, initialise_or_create_database_at


# ==========================
# USER SETTINGS
# ==========================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v47_08_12_25.db"
run_id  = 131

# quale mappa vuoi plottare (come in QCoDeS plottr)
Z_NAME = "signal_shift_Vx_deriv_trc"   # <-- questa è quella che vuoi tu
# Z_NAME = "signal_shift_Vx_deriv"
# Z_NAME = "signal_shift_V"
# Z_NAME = "G"

# assi (setpoints)
Y_NAME = "QDAC_ch02_dc_constant_V"     # y-axis (outer gate)
X_NAME = "QDAC_ch04_dc_constant_V"     # x-axis (inner gate)

# downsample: 1 riga ogni 2 sull'asse Y
Y_STRIDE = 2

# plotting
use_mV_axes = True        # converte assi in mV (solo se sono in V)
z_scale = 1e6             # per avere µV se Z è in V (come nel tuo colorbar "10^-6 Volt")
z_label = f"{Z_NAME} (µV)"


# ==========================
# HELPERS
# ==========================
def ravel(a):
    """Flatten robusto per array/list-of-arrays."""
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()

def find_block_with_param(pdata: dict, target: str):
    """Trova il block di get_parameter_data() che contiene target come chiave."""
    target_l = target.lower().strip()
    for dep_key, block in pdata.items():
        keys = list(block.keys())
        keys_l = [k.lower().strip() for k in keys]
        if target_l in keys_l:
            # ritorna la chiave esatta (case-sensitive) presente nel block
            exact_key = keys[keys_l.index(target_l)]
            return dep_key, block, exact_key
    return None, None, None

def reshape_to_grid(x, y, z):
    """
    Prova a ricostruire Z su griglia (Ny, Nx) usando i valori unici di x e y.
    Funziona se il dataset è una sweep regolare (come il tuo).
    """
    x = np.asarray(x); y = np.asarray(y); z = np.asarray(z)

    # unici ordinati
    x_unique = np.unique(x)
    y_unique = np.unique(y)

    Nx = x_unique.size
    Ny = y_unique.size

    if Nx * Ny != z.size:
        # fallback: prova a inferire Nx dal "run" tipico (x cambia più veloce)
        # prendiamo la lunghezza del primo blocco di x costante y
        # (cioè fino a quando y cambia)
        # se non funziona, alza errore chiaro
        change = np.where(np.diff(y) != 0)[0]
        if change.size == 0:
            raise ValueError("Non riesco a inferire la griglia: y non cambia mai.")
        Nx2 = change[0] + 1
        if z.size % Nx2 != 0:
            raise ValueError(
                f"Non riesco a fare reshape: z.size={z.size}, "
                f"Nx stimato={Nx2} non divide z.size."
            )
        Nx = Nx2
        Ny = z.size // Nx
        x_grid = x.reshape(Ny, Nx)
        y_grid = y.reshape(Ny, Nx)
        z_grid = z.reshape(Ny, Nx)
        return x_grid, y_grid, z_grid

    # mappa z su griglia usando lookup (robusto anche se l'ordine non è perfetto)
    xi = {val: i for i, val in enumerate(x_unique)}
    yi = {val: i for i, val in enumerate(y_unique)}
    z_grid = np.full((Ny, Nx), np.nan, dtype=float)

    for xv, yv, zv in zip(x, y, z):
        z_grid[yi[yv], xi[xv]] = zv

    # costruiamo griglie X,Y con meshgrid degli unici
    Xg, Yg = np.meshgrid(x_unique, y_unique)
    return Xg, Yg, z_grid


# ==========================
# MAIN
# ==========================
initialise_or_create_database_at(db_path)
ds = load_by_id(run_id)

pdata = ds.get_parameter_data()

dep_key, block, Z_key = find_block_with_param(pdata, Z_NAME)
if block is None:
    raise KeyError(
        f"Parametro '{Z_NAME}' non trovato in get_parameter_data(). "
        f"Prova a stampare le keys per vedere come si chiama davvero."
    )

# Estrai Z e (dallo stesso blocco) i setpoints X,Y
# In molti dataset QCoDeS, block[Z_key] è un dict-like con setpoints inclusi,
# ma spesso è direttamente array. Quindi cerchiamo X_NAME e Y_NAME nel block.
keys = list(block.keys())
if X_NAME not in block or Y_NAME not in block:
    raise KeyError(
        f"Nel blocco che contiene '{Z_NAME}' non trovo i setpoint.\n"
        f"Keys nel block: {keys}\n"
        f"Mi servono: '{X_NAME}' e '{Y_NAME}'."
    )

Z = ravel(block[Z_key]).astype(float)
X = ravel(block[X_NAME]).astype(float)
Y = ravel(block[Y_NAME]).astype(float)

# pulizia
m = np.isfinite(X) & np.isfinite(Y) & np.isfinite(Z)
X, Y, Z = X[m], Y[m], Z[m]

# scala unità
if use_mV_axes:
    X_plot = X * 1e3
    Y_plot = Y * 1e3
    x_unit = "mV"
    y_unit = "mV"
else:
    X_plot = X
    Y_plot = Y
    x_unit = "V"
    y_unit = "V"

Z_plot = Z * z_scale

# reshape in griglia
Xg, Yg, Zg = reshape_to_grid(X_plot, Y_plot, Z_plot)

# prendi una riga ogni 2 sull'asse Y
Xg_ds = Xg[::Y_STRIDE, :]
Yg_ds = Yg[::Y_STRIDE, :]
Zg_ds = Zg[::Y_STRIDE, :]

# plot stile “plottr”
plt.figure(figsize=(8.8, 4.8))
im = plt.pcolormesh(Xg_ds, Yg_ds, Zg_ds, shading="auto")
plt.xlabel(f"ch4 ({x_unit})")
plt.ylabel(f"ch2 ({y_unit})")
plt.title(f"run {run_id} | {Z_NAME} (downsample: 1 riga ogni {Y_STRIDE} in Y)")
cbar = plt.colorbar(im)
cbar.set_label(z_label)
plt.tight_layout()
plt.show()
