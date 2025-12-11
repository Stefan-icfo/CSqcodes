from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from dataprocessing.extract_fkts import *
import time


###########################
# USER SETTINGS
###########################

qc.config["core"]["db_location"] = r"D:\databases CD12_B5_F4\CD12_B5_F4v24_01_11_25.db"

# choose runs
run_ids = list(range(124, 144))   # 5e
# run_ids = list(range(357, 377))   # 10e
#run_ids = list(range(590, 610))     # 15e
run_ids = list(range(1056, 1076)) # 25e
run_ids = list(range(1023, 1043)) # 24e
run_ids = list(range(990, 1009)) # 23e
run_ids = list(range(957, 976)) # 22e
#qc.config["core"]["db_location"] = r"D:\databases CD12_B5_F4\CD12_B5_F4v19_211025.db"#150M mode

#run_ids = list(range(572, 577))


qc.config["core"]["db_location"] = r"D:\databases CD12_B5_F4\CD12_B5_F4v8.db"
run_ids = list(range(1170, 1179)) 

#run_ids = [1472]
#run_ids = list(range(1808, 1812)) 
#run_ids = list(range(1854, 1858)) 


qc.config["core"]["db_location"] = r"C:\Users\sforstner\Desktop\Triton database\CD12_B5_F4v34_25_11_25.db"
run_ids = list(range(544, 563))#step1 
run_ids = list(range(576, 596)) #step2
run_ids = list(range(608, 628)) #step3

qc.config["core"]["db_location"] = ".\Data\Raw_data\CD12_B5_F4v45_05_12_25.db"
run_ids = list(range(796,806))#300mK
#run_ids = list(range(782,792))#400mK
#run_ids = list(range(768,788))#500mK
#run_ids = list(range(754,764))#600mK
run_ids = list(range(740,750))#700mK
#run_ids = [728,731,734]#1K
#run_ids = [713,768,719]#800mK



##for softening
qc.config["core"]["db_location"] = ".\Data\Raw_data\CD12_B5_F4v39_29_11_25.db"
run_ids = [596,597]#right side of cb
run_ids = [540,541]#left side of cb
run_ids = [516,517]#far left side of cb
run_ids = [556,557]#close left side of cb
run_ids = [524,525]#left side of cb


shortening = 4                      # use only first N runs (or None)
lags = 500                          # number of ACF points to compute
max_lag_fit = 200                   # number of lags to use for fitting
n_discard = 5                      # discard these first lags from fit
custom_prefix = "137MHz700mK"


###########################
# FAST FFT AUTOCORR
###########################

def fft_autocorr(x, n_lags=None):
    """
    Fast normalized positive-lag autocorrelation using FFT.
    Returns ACF for lag >= 0.
    """
    x = x - np.mean(x)
    N = len(x)

    # zero-pad to >= 2N
    nfft = 1
    while nfft < 2 * N:
        nfft *= 2

    fx = np.fft.fft(x, nfft)
    ps = fx * np.conjugate(fx)
    acf_full = np.fft.ifft(ps).real

    acf = acf_full[:N]          # positive lags
    acf = acf / acf[0]          # normalize ACF(0) = 1

    if n_lags is not None:
        acf = acf[:n_lags]

    return acf


###########################
# LOAD AND CONCATENATE RUNS
###########################

if shortening is not None:
    run_ids = run_ids[:shortening]

x_full = []
y_full = []
time_full = []
cumulative_offset = 0.0
base_dt = None

for run_id in run_ids:
    print(f"[{run_id}] loading x,y data...")

    t_x, x_data = extract_1d(
        run_id, data_1d_name="x", setpoint_name="time_param", plot=False
    )
    _,  y_data = extract_1d(
        run_id, data_1d_name="y", setpoint_name="time_param", plot=False
    )

    # relative time
    t_rel = t_x - t_x[0]

    # estimate dt
    if len(t_rel) >= 2:
        dt = float(np.median(np.diff(t_rel)))
        if base_dt is None:
            base_dt = dt
        else:
            if abs(dt - base_dt) > base_dt * 1e-3:
                print(f"WARNING: dt mismatch in run {run_id}: "
                      f"this {dt:.6g}s vs base {base_dt:.6g}s")

    # shifted time for concatenation
    t_shift = t_rel + cumulative_offset

    x_full.append(x_data.astype(float))
    y_full.append(y_data.astype(float))
    time_full.append(t_shift)

    cumulative_offset = t_shift[-1] + dt

# concatenate
x_full = np.concatenate(x_full)
y_full = np.concatenate(y_full)
time_full = np.concatenate(time_full)

print(f"Total samples: {len(x_full)}, total duration: {time_full[-1]:.6f} s")
print(f"dt = {base_dt} s")


###########################
# ENERGY-LIKE ACF: R_E = R_x + R_y
###########################

print("Computing ACF via FFT...")

acf_x = fft_autocorr(x_full, n_lags=lags)
acf_y = fft_autocorr(y_full, n_lags=lags)

# This is really Re⟨a(t)a*(t+τ)⟩ with a = x + i y
auto = acf_x + acf_y
auto = auto / auto[0]     # normalize so R_E(0) = 1

###########################
# DAMPED COSINE FIT
###########################

lags_arr = np.arange(max_lag_fit)
t_lags = lags_arr * base_dt          # seconds
t_lags_ms = t_lags * 1e3             # ms

xFit = t_lags_ms[n_discard:max_lag_fit]
yFit = auto[n_discard:max_lag_fit]

def damped_cosine(x_ms, A, tau_ms, C0, f_per_ms, phi):
    """
    Damped cosine model:
    A * exp(-x/tau_ms) * cos(2π f_per_ms x + phi) + C0
    x_ms in ms, f_per_ms in cycles/ms (i.e. kHz).
    """
    return A * np.exp(-x_ms / tau_ms) * np.cos(2*np.pi*f_per_ms*x_ms + phi) + C0

# --- initial guesses ---

# offset from tail
C0_0 = np.mean(yFit[int(0.7 * len(yFit)):])   # mean of last 30% as rough offset

# amplitude from initial point
A0 = yFit[0] - C0_0

# tau: choose about half of fit window as rough guess
tau0_ms = (xFit[-1] - xFit[0]) / 2.0 if xFit[-1] > xFit[0] else 0.1

# rough frequency guess via FFT of (yFit - C0)
y_centered = yFit - C0_0
dt_ms = xFit[1] - xFit[0]
Y = np.fft.rfft(y_centered)
freqs = np.fft.rfftfreq(len(y_centered), d=dt_ms)  # cycles/ms
if len(freqs) > 1:
    idx_max = np.argmax(np.abs(Y[1:])) + 1  # ignore DC
    f0 = freqs[idx_max]
else:
    f0 = 0.0

phi0 = 0.0

p0 = [A0, tau0_ms, C0_0, f0, phi0]

print("Initial guesses:",
      f"A0={A0:.3g}, tau0={tau0_ms:.3g} ms, C0={C0_0:.3g}, f0={f0:.3g} 1/ms")

# bounds: tau>0, f>=0, phase in [-2π,2π]
bounds_lower = [-np.inf, 0.0, -np.inf, 0.0, -2*np.pi]
bounds_upper = [ np.inf, np.inf,  np.inf, np.inf,  2*np.pi]
try:
    popt, _ = curve_fit(
    damped_cosine,
    xFit, yFit,
    p0=p0,
    bounds=(bounds_lower, bounds_upper),
    )


    A_fit, tau_ms_fit, C0_fit, f_per_ms_fit, phi_fit = popt

    print(f"Fit results:")
    print(f"  tau   = {tau_ms_fit:.4g} ms")
    print(f"  fbeat = {f_per_ms_fit:.4g} 1/ms  (~ {f_per_ms_fit:.4g} kHz)")
    print(f"  A     = {A_fit:.4g}")
    print(f"  C0    = {C0_fit:.4g}")
    print(f"  phi   = {phi_fit:.4g} rad")

except:
    print("no fit possible")
###########################
# PLOT
###########################

fig, ax = plt.subplots()
ax.plot(t_lags_ms, auto[:max_lag_fit], "g*", label="ACF (Energy-like)")
try:
    ax.plot(xFit, damped_cosine(xFit, *popt),
        label=rf"Fit $\tau = {tau_ms_fit:.3g}$ ms, $f={f_per_ms_fit:.3g}$ kHz")
except:
    print("no fit")
ax.set_xlabel("Lag time (ms)")
ax.set_ylabel("Autocorrelation")
ax.grid(True)
ax.legend()
plt.title(custom_prefix
          + qc.config["core"]["db_location"][-15:-3]
          + f"_{run_ids[0]}-{run_ids[-1]}_{round(len(x_full)/1e3)}k samples")
plt.show()
