import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# === Impostazioni matplotlib per usare Calibri ===
mpl.rcParams['font.family'] = 'Calibri'

# --- DB + run ---
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD13_E3_C2.db"
initialise_or_create_database_at(db_path)
run_id = 302 # cambia se vuoi un altro ID
ds = load_by_id(run_id)

# --- Estrazione dati robusta ---
pdata = ds.get_parameter_data()

if "G" in pdata:
    dep_key = "G"
else:
    dep_key = next(iter(pdata.keys()))

block = pdata[dep_key]
setpoint_key = next(k for k in block.keys() if k != dep_key)

def _to_1d(arr):
    if isinstance(arr, (list, tuple)) and len(arr) and hasattr(arr[0], "__len__"):
        return np.concatenate([np.asarray(a).ravel() for a in arr])
    return np.asarray(arr).ravel()

x = _to_1d(block[setpoint_key])     # Gate voltage
y = _to_1d(block[dep_key]) * 1e6    # Convert conductance to µS

# --- Plot ---
plt.figure(figsize=(6,4))
plt.plot(x, y, "-", color="black", lw=1.4)   # solo linea nera
plt.xticks(fontfamily="Calibri", fontsize=16)
plt.yticks(fontfamily="Calibri", fontsize=16)
plt.xlabel("Voltage GCS", fontsize=18)
plt.ylabel(r"G ($\mu$S)", fontsize=18)

plt.tight_layout()
plt.show()







