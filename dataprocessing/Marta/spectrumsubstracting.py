import os
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# -----------------------------
# CONFIG
# -----------------------------
run_data = 542
run_background = 540
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v9.db"
out_dir = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Figures"
os.makedirs(out_dir, exist_ok=True)
qc.config["core"]["db_location"] = db_path

# -----------------------------
# Utilities
# -----------------------------
def list_keys(run_id):
    ds = qc.load_by_id(run_id)
    pd = ds.get_parameter_data()
    return list(pd.keys()), pd

def choose_freq_key(block_dict, dep_key):
    """Find the correct frequency setpoint key inside pdata[dep_key]."""
    candidates = [k for k in block_dict.keys() if k != dep_key]
    # Prefer keys containing 'freq'
    freq_like = [k for k in candidates if 'freq' in k.lower()]
    for k in freq_like + candidates:
        arr = np.asarray(block_dict[k]).flatten()
        if arr.size > 0 and np.isfinite(arr).any():
            return k
    # Fallback: just return any candidate
    return candidates[0] if candidates else None

def find_psd_key(pdata):
    """Pick a PSD-like dependent key by priority."""
    priorities = ["avg_avg_psd_nodrive", "avg_psd"]
    for k in priorities:
        if k in pdata:
            return k
    # Any key containing 'psd'
    for k in pdata.keys():
        if "psd" in k.lower():
            return k
    # As a last resort, return the longest array block
    return max(pdata.keys(), key=lambda k: np.size(pdata[k][k])) if pdata else None

def load_1d_psd(run_id):
    """
    Robustly load (freq_MHz, PSD_W/Hz) for a run.
    - If the chosen PSD block is time-resolved, average over freq bins.
    - If the block is empty/NaN, try fallbacks.
    """
    keys, pdata = list_keys(run_id)
    if not pdata:
        raise ValueError(f"[Run {run_id}] No parameter data in dataset.")

    tried = []
    def attempt_with(dep_key):
        tried.append(dep_key)
        block = pdata[dep_key]
        freq_key = choose_freq_key(block, dep_key)
        if freq_key is None:
            raise KeyError(f"[Run {run_id}] Could not find a setpoint key for '{dep_key}'.")
        y_raw = np.asarray(block[dep_key]).flatten()
        x_raw = np.asarray(block[freq_key]).flatten()

        # Drop non-finite early
        mask0 = np.isfinite(x_raw) & np.isfinite(y_raw)
        x_raw, y_raw = x_raw[mask0], y_raw[mask0]

        if x_raw.size == 0 or y_raw.size == 0:
            raise ValueError(f"[Run {run_id}] Empty arrays with dep='{dep_key}', freq='{freq_key}'.")

        # Heuristic: if many repeats in x, average per Hz bin (time-resolved)
        unique_rounded = np.unique(np.round(x_raw))
        if len(x_raw) > 5 and len(unique_rounded) < 0.8 * len(x_raw):
            bins = np.round(x_raw).astype(np.int64)
            sums, counts = {}, {}
            for b, yy in zip(bins, y_raw):
                if np.isfinite(yy):
                    sums[b] = sums.get(b, 0.0) + yy
                    counts[b] = counts.get(b, 0) + 1
            if not sums:
                raise ValueError(f"[Run {run_id}] All NaN in time-resolved block '{dep_key}'.")
            x_hz = np.array(sorted(sums.keys()), dtype=float)
            y = np.array([sums[b] / counts[b] for b in x_hz.astype(int)], dtype=float)
        else:
            x_hz, y = x_raw, y_raw

        # Final clean & sort
        m = np.isfinite(x_hz) & np.isfinite(y)
        x_hz, y = x_hz[m], y[m]
        if x_hz.size == 0 or y.size == 0:
            raise ValueError(f"[Run {run_id}] No finite points after cleaning for '{dep_key}'.")
        idx = np.argsort(x_hz)
        x_mhz = x_hz[idx] / 1e6
        y = y[idx]

        # Quick diag
        print(f"[Run {run_id}] dep='{dep_key}', freq_key='{freq_key}', N={len(x_mhz)}, "
              f"f=[{x_mhz.min():.3f},{x_mhz.max():.3f}] MHz, PSD=[{np.nanmin(y):.3e},{np.nanmax(y):.3e}]")
        return x_mhz, y

    # Try main choice
    dep = find_psd_key(pdata)
    try:
        return attempt_with(dep)
    except Exception as e1:
        # If the first choice failed, try other PSD-like blocks
        fallbacks = [k for k in pdata.keys() if k != dep and ("psd" in k.lower() or k.endswith("_psd"))]
        for fb in fallbacks:
            try:
                return attempt_with(fb)
            except Exception as e2:
                continue
        # As a last resort, try ANY other dependent block
        for k in pdata.keys():
            if k in tried:
                continue
            try:
                return attempt_with(k)
            except Exception:
                pass
        # If all failed, raise with context
        raise RuntimeError(
            f"[Run {run_id}] Could not load a usable PSD. "
            f"Tried: {tried + fallbacks}. Available blocks: {list(pdata.keys())}"
        )

# -----------------------------
# Load data & background
# -----------------------------
x_data, P_total = load_1d_psd(run_data)
x_bkg,  P_noise = load_1d_psd(run_background)

# -----------------------------
# Voltage-domain subtraction (point-by-point)
# -----------------------------
n = min(len(P_total), len(P_noise))
if n == 0:
    raise RuntimeError("No overlapping length after loading. Check the runs and keys.")

# Clip negatives before sqrt for numerical safety
V_total = np.sqrt(np.clip(P_total[:n], 0.0, None))
V_noise = np.sqrt(np.clip(P_noise[:n], 0.0, None))
V_signal = V_total - V_noise
P_signal = V_signal**2

x_plot = x_data[:n]
print(f"[Subtraction] n={n}, P_total[:3]={P_total[:3]}, P_noise[:3]={P_noise[:3]}, P_signal[:3]={P_signal[:3]}")

# -----------------------------
# Plot & Save
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x_data[:n], P_total[:n], label=f"Run {run_data} (raw)", alpha=0.6)
plt.plot(x_bkg[:n],  P_noise[:n], label=f"Run {run_background} (background)", alpha=0.6)
plt.plot(x_plot,     P_signal, "o", markersize=3, label="Corrected spectrum")

plt.xlabel("Frequency (MHz)")
plt.ylabel("PSD (W/Hz)")
plt.title("Corrected PSD after voltage-domain subtraction (point-by-point)")
plt.grid(True)
plt.legend()
plt.tight_layout()

png_path = os.path.join(out_dir, f"corrected_psd_run{run_data}_minus_{run_background}.png")
plt.savefig(png_path, dpi=200)
print(f"Saved figure to: {png_path}")

try:
    plt.show()
except Exception as e:
    print(f"plt.show() failed (headless backend?): {e}")


