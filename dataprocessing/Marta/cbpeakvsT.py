import re
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import matplotlib as mpl

# -----------------------------
# 0) Database + run selection
# -----------------------------
qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v7.db"
)

# Odd run IDs 683..707, but exclude 693
RUN_IDS = [rid for rid in range(683, 708, 2) if rid != 693]

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
# -----------------------------
# 2) Load, color by temperature, and plot
# -----------------------------
plt.figure(figsize=(10, 6))

curves = []  # (Vg, G, label, T_mK)

def temperature_info(ds, run_id):
    # text label + numeric value in mK (if found)
    blob = collect_texts(ds)
    m = re.search(r"(\d{2,4})\s*mk", blob, flags=re.IGNORECASE)
    if m:
        val = float(m.group(1))
        return f"{int(val)} mK", val
    if run_id in manual_temp_by_run:
        txt = manual_temp_by_run[run_id]
        m2 = re.search(r"(\d{2,4})", txt)
        val = float(m2.group(1)) if m2 else None
        return txt, val
    return f"run {run_id}", None

for rid in RUN_IDS:
    try:
        ds = qc.load_by_id(rid)
        Vg, G = get_xy_from_dataset(ds)
        label, T_mK = temperature_info(ds, rid)
        curves.append((Vg, G, label, T_mK))
    except Exception as e:
        print(f"[WARN] Skipped run {rid}: {e}")

if not curves:
    raise RuntimeError("No runs loaded. Check DB path, run IDs, or parameter names.")

# Colormap: orange -> dark red; higher T = darker red
temps = [T for *_ , T in curves if T is not None]
vmin, vmax = (min(temps), max(temps)) if temps else (0, 1)
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cmap = mpl.cm.get_cmap("OrRd")  # orange to dark red

for Vg, G, label, T_mK in curves:
    color = "gray" if T_mK is None else cmap(norm(T_mK))
    plt.plot(Vg, G*1e6, marker="o", linestyle="-", markersize=2.5, linewidth=1.0,
             color=color, label=label)

# Labels & colorbar (legend optional; colorbar encodes T)

plt.xlabel("Gate voltage (V)", fontsize=18)
plt.ylabel("Conductance G ($\mu$S)", fontsize=18)
plt.title("GVG vs Temperature", fontsize=18)
mpl.rcParams["xtick.labelsize"] = 16         # numeri asse x
mpl.rcParams["ytick.labelsize"] = 16
sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
cbar = plt.colorbar(sm, pad=0.02)
cbar.set_label("Temperature (mK)", fontsize=18)
plt.tight_layout()
plt.show()
