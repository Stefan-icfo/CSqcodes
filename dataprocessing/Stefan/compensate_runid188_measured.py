import numpy as np
import matplotlib.pyplot as plt
import scipy as scp
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
qc.config["core"]["db_location"] = r'D:\databases CD12_B5_F4\CD12_B5_F4v49_15_12_25.db'
print("Opening DB:", qc.config["core"]["db_location"])

run_id = 188

import os
path = os.path.join(os.path.expanduser("~"), "Desktop", "compressed_filter.csv")
compressed_filter = np.loadtxt(path)

freq, spec = extract_1d(
    run_id,
    data_1d_name="V_fft_avg_avg",
    setpoint_name="freq_param",
    plot=False, return_exp_name=False,
)
freq = np.asarray(freq, float)
spec = np.asarray(spec, float)

filt = np.asarray(compressed_filter, float)   # the array you extracted from Zurich

# ---- put the filter on the same grid as the spectrum ----
if filt.size == spec.size:
    H = filt.copy()
else:
    # both span the same (centred) band -> stretch the filter onto the spectrum grid
    H = np.interp(np.linspace(0, 1, spec.size),
                  np.linspace(0, 1, filt.size), filt)

H = H / H.max()              # normalise peak to 1
compensated = spec / H       # V_fft is an amplitude -> divide by |H|. Use H**2 if it's a power/PSD filter.

# ---- plot ----
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(freq, spec,        lw=0.7, label="V_fft_avg_avg (raw)")
ax.plot(freq, compensated, lw=0.9, label="compensated")
ax.set_xlabel("freq_param [Hz]"); ax.set_ylabel("PSD [a.u.]")
ax.set_title(f"run_id {run_id} — divided by measured filter")
ax.legend(); plt.tight_layout(); plt.show()

print(f"filter length {filt.size}, spectrum length {spec.size}")
print(f"filter range: {filt.min():.4f} .. {filt.max():.4f}")
