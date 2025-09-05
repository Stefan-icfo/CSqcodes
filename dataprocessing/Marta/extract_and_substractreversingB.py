import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc

# -----------------------------
# CONFIG: update run IDs and DB path
# -----------------------------
run_data = 877        # example
run_background = 716  # example
qc.config["core"]["db_location"] = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v2.db"

# -----------------------------
# Helpers
# -----------------------------
def load_sorted(run_id, signal_key="avg_avg_psd_nodrive", freq_key="freq_param"):
    """
    Load (x, y) from QCoDeS for run_id.
    Return frequency in MHz sorted ASCENDING together with PSD.
    """
    ds = qc.load_by_id(run_id)
    pdict = ds.get_parameter_data()
    y = pdict[signal_key][signal_key].flatten()
    x = (pdict[signal_key][freq_key].flatten()) / 1e6  # Hz -> MHz
    idx = np.argsort(x)
    return x[idx], y[idx]

def reverse_x_only(x, y):
    """
    Reverse the frequency array only (high->low) while KEEPING y order unchanged.
    (x1,y1)->(x_last,y1), (x2,y2)->(x_last-1,y2), ...
    """
    return x[::-1].copy(), y.copy()

# -----------------------------
# Load and reverse-x both runs
# -----------------------------
x_d, y_d = load_sorted(run_data)
x_b, y_b = load_sorted(run_background)

x_d_rev, y_d_rev = reverse_x_only(x_d, y_d)
x_b_rev, y_b_rev = reverse_x_only(x_b, y_b)

# -----------------------------
# Mirror TRANSFORM IN Y for the BACKGROUND around its own maximum
# -----------------------------
y_bg_max = float(np.max(y_b_rev))        # horizontal mirror level
y_b_mir  = 2.0 * y_bg_max - y_b_rev      # mirror in Y

# -----------------------------
# Plot: background (x reversed), mirror LINE (horizontal), mirrored background (Y mirrored)
# -----------------------------
plt.figure(figsize=(11, 6))
plt.plot(x_b_rev, y_b_rev, '-', lw=1.8, label=f'Background (x reversed)')
plt.axhline(y_bg_max, linestyle='--', lw=1.6, label=f'Mirror line @ y = {y_bg_max:.3e}')
plt.plot(x_b_rev, y_b_mir, '-', lw=1.8, alpha=0.9, label='Background mirrored IN Y')

plt.xlabel("Frequency (MHz)", fontsize=13)
plt.ylabel("PSD (W/Hz)", fontsize=13)
plt.title("Background: x reversed; Y mirrored around its maximum", fontsize=15)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

