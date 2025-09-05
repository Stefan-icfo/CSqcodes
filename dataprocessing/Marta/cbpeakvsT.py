import re
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# -----------------------------
# 0) Database + run selection
# -----------------------------
qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v7.db"
)

# Odd run IDs 683..707, but exclude 693
RUN_IDS = [rid for rid in range(681, 708, 2) if rid != 693]

# Parameter names (setpoint and traces saved by GVG_fun_sensitivity)
g_key = "G"               # conductance for y-axis
vg_key = None             # detect automatically

# Optional fallback mapping if a temperature is not in the title
manual_temp_by_run = {
    683: "550 mK",
    685: "750 mK",
}

# -----------------------------
# 1) Helpers
# -----------------------------
def extract_vg_key(param_data_block):
    for k in param_data_block.keys():
        if k != g_key:
            return k
    raise KeyError("Could not find gate voltage key.")

def get_xy_from_dataset(ds):
    pd = ds.get_parameter_data()
    if g_key not in pd:
        blk_name = None
        for name, blk in pd.items():
            if g_key in blk:
                blk_name = name
                break
        if blk_name is None:
            raise KeyError(f"'{g_key}' not present in dataset.")
        block = pd[blk_name]
    else:
        block = pd[g_key]

    vg_key_local = extract_vg_key(block)
    y = np.asarray(block[g_key]).flatten()
    x = np.asarray(block[vg_key_local]).flatten()
    order = np.argsort(x)
    return x[order], y[order]

def collect_texts(ds):
    texts = []
    for attr in ("name", "run_name", "exp_name", "sample_name", "snapshot_raw", "snapshot"):
        try:
            v = getattr(ds, attr)
            if isinstance(v, str) and v:
                texts.append(v)
        except Exception:
            pass
    try:
        md = ds.get_metadata()
        if isinstance(md, dict):
            for k in ("name", "exp_name", "sample_name", "description"):
                val = md.get(k)
                if isinstance(val, str):
                    texts.append(val)
    except Exception:
        pass
    return "  ".join(texts)

def temperature_label(ds, run_id):
    blob = collect_texts(ds)
    m = re.search(r"(\d{2,4})\s*mk", blob, flags=re.IGNORECASE)
    if m:
        return f"{int(m.group(1))} mK"
    if run_id in manual_temp_by_run:
        return manual_temp_by_run[run_id]
    return f"run {run_id}"

# -----------------------------
# 2) Load & plot
# -----------------------------
plt.figure(figsize=(10, 6))
labels_used = set()
loaded_any = False

for rid in RUN_IDS:
    try:
        ds = qc.load_by_id(rid)
        Vg, G = get_xy_from_dataset(ds)
        label = temperature_label(ds, rid)

        if label in labels_used:
            label = f"{label} (run {rid})"
        labels_used.add(label)

        plt.plot(Vg, G, marker="o", linestyle="-", markersize=2.5, linewidth=1.0, label=label)
        loaded_any = True
    except Exception as e:
        print(f"[WARN] Skipped run {rid}: {e}")

if not loaded_any:
    raise RuntimeError("No runs loaded. Check DB path, run IDs, or parameter names.")

plt.xlabel("Gate voltage (V)", fontsize=13)
plt.ylabel("Conductance G (S)", fontsize=13)
plt.title("GVGvsT", fontsize=15)
plt.grid(True, alpha=0.3)
plt.legend(title="Temperature", fontsize=10)
plt.tight_layout()
plt.show()
