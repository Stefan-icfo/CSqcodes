import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# === Font Calibri per tutto (se installato sul sistema) ===
mpl.rcParams['font.family'] = 'Calibri'

# === Database ===
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD13_E3_C2.db"
initialise_or_create_database_at(db_path)

# === Utility: estrai (x, y) in modo robusto per un run id ===
def get_xy(run_id, dep_name_prefer="G"):
    ds = load_by_id(run_id)
    pdata = ds.get_parameter_data()
    dep_key = dep_name_prefer if dep_name_prefer in pdata else next(iter(pdata.keys()))
    block = pdata[dep_key]
    setpoint_key = next(k for k in block.keys() if k != dep_key)

    def _to_1d(arr):
        if isinstance(arr, (list, tuple)) and len(arr) and hasattr(arr[0], "__len__"):
            return np.concatenate([np.asarray(a).ravel() for a in arr])
        return np.asarray(arr).ravel()

    x = _to_1d(block[setpoint_key])
    y = _to_1d(block[dep_key]) * 1e6  # converti G in µS
    # ordina per x (così le linee sono pulite anche se lo sweep ha salti)
    order = np.argsort(x)
    return x[order], y[order], setpoint_key, dep_key

# === Specifica i run e le etichette (temperatura stimata) ===
runs = [
    (299, "~30 K", "black"),
    (300, "~20 K", "red"),
    (302, "~10 K", "blue"),
]

# === Plot ===
plt.figure(figsize=(7,4.5))

x_label_from = None
y_label_from = None

for run_id, tlabel, color in runs:
    x, y, xname, yname = get_xy(run_id)
    if x_label_from is None:
        x_label_from = xname
    if y_label_from is None:
        y_label_from = yname
    plt.plot(x, y, "-", lw=1.6, color=color, label=f"Run {run_id} ({tlabel})")

# Etichette assi (Calibri) e tick labels in Calibri
plt.xlabel("Voltage GCS", fontsize=18, fontfamily="Calibri")
plt.ylabel(r"G ($\mu$S)", fontsize=18, fontfamily="Calibri")
plt.xticks(fontfamily="Calibri", fontsize=16)
plt.yticks(fontfamily="Calibri", fontsize=16)

# Legenda
leg = plt.legend(frameon=False, fontsize=12)
for text in leg.get_texts():
    text.set_fontfamily("Calibri")

plt.tight_layout()
plt.show()







