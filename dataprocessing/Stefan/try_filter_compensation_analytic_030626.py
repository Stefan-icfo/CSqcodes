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

# ---- the dataset you actually meant ----
freq, spec = extract_1d(
    run_id,
    data_1d_name="V_fft_avg_avg",
    setpoint_name="freq_param",
    plot=False, return_exp_name=False,
)
freq = np.asarray(freq, dtype=float)
spec = np.asarray(spec, dtype=float)

# ---- locate the narrow signal peak so we exclude it from the envelope fit ----
meta   = get_metadata(run_id - 1, print_it=False, return_data=True)  # your loop's convention
f_sig  = meta["center_freq"]     # narrow-peak location
MASK_HW = 10e3                   # exclude +/- this window around the peak [Hz]; ~ a few * linewidth
fit_mask = np.abs(freq - f_sig) > MASK_HW

# ---- fit the broad envelope (filter/tank shape) on the masked data ----
# generalised Lorentzian/RC form: flexible enough for either origin of the hump
def envelope(f, A, f0, hw, p, c):
    return A / (1.0 + ((f - f0) / hw) ** 2) ** p + c

f_mid = 0.5 * (freq.min() + freq.max())
p0 = [np.nanmax(spec), f_mid, (freq.max() - freq.min()) / 4, 1.0, np.nanmin(spec)]
popt, _ = scp.optimize.curve_fit(
    envelope, freq[fit_mask], spec[fit_mask], p0=p0, maxfev=20000
)
env = envelope(freq, *popt)

compensated = spec / env          # flat background ~1, peak sticks up

# ---- plot: raw + fitted envelope, then the flattened result ----
fig, axs = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
axs[0].plot(freq, spec, lw=0.7, label="V_fft_avg_avg")
axs[0].plot(freq, env, "r", lw=1.2, label="fitted envelope")
axs[0].plot(freq[~fit_mask], spec[~fit_mask], "k.", ms=2, label="excluded (peak)")
axs[0].set_ylabel("PSD [a.u.]"); axs[0].legend()
axs[1].plot(freq, compensated, lw=0.7)
axs[1].axhline(1.0, color="r", lw=0.5)
axs[1].set_ylabel("compensated"); axs[1].set_xlabel("freq_param [Hz]")
axs[0].set_title(f"run_id {run_id} — envelope removed")
plt.tight_layout(); plt.show()


# =====================================================================
# If you'd rather keep the ANALYTIC demod-filter division, the only fix
# is to centre it on the hump (the demod LO), not on 0, and use the BW
# actually set for this run:
#
#   N    = 3
#   F3DB = 60e3                       # NOT 100e3 — match the real width
#   DEMOD_CENTER = freq[np.argmax(envelope(freq, *popt))]  # ~hump centre
#   tau  = np.sqrt(2**(1/N) - 1) / (2*np.pi*F3DB)
#   H2   = (1 + (2*np.pi*(freq - DEMOD_CENTER)*tau)**2)**(-N)
#   compensated = spec / H2
# =====================================================================