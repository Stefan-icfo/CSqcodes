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
    Load (x, y) from QCoDeS for a run_id.
    Returns frequency in MHz sorted ASCENDING together with PSD.
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
    So (x1, y1) becomes (x_last, y1), (x2, y2) -> (x_last-1, y2), etc.
    """
    return x[::-1].copy(), y.copy()

def truncate_pair(a, b):
    """Truncate two arrays to the same (shorter) length."""
    n = min(len(a), len(b))
    return a[:n], b[:n]

# -----------------------------
# Load and prepare data
# -----------------------------
x_d, y_d = load_sorted(run_data)
x_b, y_b = load_sorted(run_background)

# Reverse ONLY the frequency axis for each run
x_d_rev, y_d_rev = reverse_x_only(x_d, y_d)
x_b_rev, y_b_rev = reverse_x_only(x_b, y_b)

# For the difference, align STRICTLY by index (first with first, etc.)
y_d_aligned, y_b_aligned = truncate_pair(y_d_rev, y_b_rev)
x_diff = x_d_rev[:len(y_d_aligned)]  # use data run's reversed x for plotting the diff
y_diff = y_d_aligned - y_b_aligned   # point-by-point difference

# -----------------------------
# Plot 1: spectra with ONLY x reversed
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x_d_rev, y_d_rev, '-', lw=1.6, label=f'Run {run_data} (x reversed only)')
plt.plot(x_b_rev, y_b_rev, '-', lw=1.6, alpha=0.85, label=f'Run {run_background} (x reversed only)')
plt.xlabel("Frequency (MHz)", fontsize=13)
plt.ylabel("PSD (W/Hz)", fontsize=13)
plt.title("Spectra with x reversed only (index 1 stays index 1)", fontsize=15)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# -----------------------------
# Plot 2: point-by-point difference (index-aligned) on reversed-x grid
# -----------------------------
plt.figure(figsize=(10, 6))
plt.plot(x_diff, y_diff, '-', lw=1.6, label=f'{run_data} - {run_background} (index-aligned, x reversed)')
plt.xlabel("Frequency (MHz)", fontsize=13)
plt.ylabel("PSD difference (W/Hz)", fontsize=13)
plt.title("Difference spectrum with ONLY x reversed (first↔first, …)", fontsize=15)
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
