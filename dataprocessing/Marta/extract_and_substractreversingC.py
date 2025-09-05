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
    ds = qc.load_by_id(run_id)
    pdict = ds.get_parameter_data()
    y = pdict[signal_key][signal_key].flatten()
    x = pdict[signal_key][freq_key].flatten() / 1e6  # Hz -> MHz
    idx = np.argsort(x)
    return x[idx], y[idx]

def reverse_x_only(x, y):
    return x[::-1].copy(), y.copy()

# -----------------------------
# Load runs
# -----------------------------
x_d, y_d = load_sorted(run_data)
x_b, y_b = load_sorted(run_background)

# Data reversed (B)
x_d_rev, y_d_rev = reverse_x_only(x_d, y_d)

# Background reversed + Y mirrored (A)
x_b_rev, y_b_rev = reverse_x_only(x_b, y_b)
y_bg_max = float(np.max(y_b_rev))
y_b_mir = 2.0 * y_bg_max - y_b_rev

# -----------------------------
# Sum A + B (index aligned)
# -----------------------------
n = min(len(y_d_rev), len(y_b_mir))
x_plot = x_d_rev[:n]            # frequency from data reversed
y_sum  = y_b_mir[:n]* y_d_rev[:n]

# -----------------------------
# Plot: sum (blue) and data reversed (red)
# -----------------------------
plt.figure(figsize=(11, 6))
plt.plot(x_plot, y_sum, color='blue', lw=1.8, label='A * B ')
#plt.plot(x_plot, y_d_rev[:n], color='red', lw=1.4, label='Data reversed (B)')
plt.xlabel("Frequency (MHz) [using data reversed]", fontsize=13)
plt.ylabel("PSD / Sum (arb. units)", fontsize=13)
plt.title("Sum: (Background reversed & Y-mirrored) + (Data reversed)", fontsize=15)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()




