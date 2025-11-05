import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import sqlite3
import re

# ============== CONFIG (v11) ==============
qc.config["core"]["db_location"] = (
    r"W:\\Electromechanics\\Projects\\chargesensor\\triton measurements\\Raw_data\\CD12_B5_F4v11.db"
)
background_id = 2885
all_data_ids = list(range(2948, 3033, 2))   # 2948, 2950, ..., 3032

# ============== HELPERS ==============
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

def run_exists(run_id: int) -> bool:
    db = qc.config["core"]["db_location"]
    with sqlite3.connect(db) as con:
        row = con.execute("SELECT 1 FROM runs WHERE run_id=?", (run_id,)).fetchone()
    return row is not None

# ============== VALIDAZIONI ==============
db_path = qc.config["core"]["db_location"]
print("DB:", db_path)

if not run_exists(background_id):
    raise RuntimeError(f"Background run_id {background_id} NOT found in this DB. "
                       f"Apri il DB giusto o cambia background_id.")

valid_ids = [rid for rid in all_data_ids if run_exists(rid)]
if not valid_ids:
    raise RuntimeError("Nessun run di all_data_ids esiste in questo DB. "
                       "Controlla la lista in base al DB attuale (v11: 2948..3032).")
if len(valid_ids) != len(all_data_ids):
    print(f"⚠️ Alcuni run non esistono e verranno saltati. Userò: {valid_ids}")

# ============== CARICA BACKGROUND ==============
x_bg, y_bg, _ = load_psd(background_id)

# ============== LOOP E PLOT ==============
spectra, drives, valid_runs = [], [], []
x_common = None

for run_id in valid_ids:
    try:
        x, y, ds = load_psd(run_id)
        drive = extract_drive_amplitude(ds.exp_name)
        if drive is None:
            print(f"⚠️ Run {run_id}: drive non trovato in exp_name '{ds.exp_name}' → skip")
            continue

        # subtract bg punto-per-punto
        n = min(len(y), len(y_bg))
        y_corr = y[:n] - y_bg[:n]
        x_corr = x[:n]

        if x_common is None:
            x_common = x_corr

        m = min(len(x_common), len(y_corr))
        x_corr = x_corr[:m]; y_corr = y_corr[:m]

        spectra.append(y_corr)
        drives.append(drive)
        valid_runs.append(run_id)

    except Exception as e:
        print(f"⚠️ Run {run_id}: skipped → {e}")

if not spectra:
    raise RuntimeError("❌ No valid spectra loaded.")

spectra = np.array(spectra)
drives = np.array(drives)

# ordina per drive
order = np.argsort(drives)
spectra = spectra[order]
drives = drives[order]

# heatmap
plt.figure(figsize=(12, 6))
extent = [x_common[0], x_common[-1], drives[0], drives[-1]]
plt.imshow(spectra, aspect="auto", extent=extent, origin="lower", cmap="plasma")
plt.clim(0, 2e-6)
plt.xlabel("Frequency (MHz)")
plt.ylabel("Drive amplitude (µV)")
plt.title("Background-subtracted spectra (point-by-point)")
plt.colorbar(label="PSD (W/Hz)")
plt.tight_layout()
plt.show()

print(f"✅ Loaded {len(valid_runs)} runs: {valid_runs}")
print(f"Drive amplitudes (µV): {drives}")
