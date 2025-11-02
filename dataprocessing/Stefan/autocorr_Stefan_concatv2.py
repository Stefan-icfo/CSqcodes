
from tqdm import tqdm
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, hbar, e, h
from scipy.constants import Boltzmann as kB
from scipy.optimize import curve_fit
from dataprocessing.extract_fkts import *
import time
qc.config["core"]["db_location"] = ".\Data\Raw_data\CD12_B5_F4v24_01_11_25.db"





run_ids = list(range(324, 344))  # 
lags=500

# Define your fitting function
def fitting_func(x, a, b, c, d, e):
    return a * np.exp(-b * x) + c + d * np.exp(-e * x)

# List to store results
results = []

# List to store figures for later plotting
figures = []
overall_start_time=time.time_ns()
var_data=[]
time_data=[]
# Collect all runs into one continuous time series
var_data_list = []
time_data_list = []

cumulative_offset = 0.0  # seconds
base_dt = None           # reference timestep (seconds)

for run_id in run_ids:
    try:
        # Extract data for each run_id
        t_data, phase_data = extract_1d(
            run_id, data_1d_name="Phase", setpoint_name="time_param", plot=False
        )
        print(f"[{run_id}] extracted phase data")

        # Relative time for this run, starting at 0
        t_rel = (t_data - t_data[0]).astype(float)

        # Estimate timestep for this run
        if len(t_rel) >= 2:
            dts = np.diff(t_rel)
            dt = float(np.median(dts))
        else:
            print(f"[{run_id}] not enough samples; skipping.")
            continue

        # Check/lock a common sampling interval
        if base_dt is None:
            base_dt = dt
        else:
            # warn if mismatch (tolerance: 0.1% of base_dt)
            if abs(dt - base_dt) > (1e-3 * base_dt + 1e-12):
                print(f"[{run_id}] WARNING: dt mismatch (this {dt:.6g}s vs base {base_dt:.6g}s). "
                      "ACF assumes uniform sampling.")

        # Shift this run's time so it starts after the previous run
        t_shift = t_rel + cumulative_offset

        # Append data
        time_data_list.append(t_shift)
        var_data_list.append(phase_data.astype(float))

        # Update offset so next run starts one dt later (avoids duplicate end/start time)
        cumulative_offset = t_shift[-1] + dt

    except RuntimeError as e:
        print(f"[{run_id}] RuntimeError: {e}. Skipping this run and continuing.")
    except Exception as e:
        print(f"[{run_id}] Error: {e}. Skipping this run and continuing.")

# Concatenate across runs
if len(time_data_list) == 0:
    raise RuntimeError("No valid runs to concatenate.")
time_data = np.concatenate(time_data_list)
var_data  = np.concatenate(var_data_list)

print(f"Concatenated {len(time_data)} samples across {len(time_data_list)} runs.")
print(f"Estimated base dt ~ {base_dt:.6g} s, total duration ~ {time_data[-1]-time_data[0]:.6g} s.")

        # Calculate autocorrelation function

num = len(time_data)
auto = np.zeros(num)

extracted_data_abs_time=time.time_ns()
print(f"extracted_data_time_taken={(extracted_data_abs_time-overall_start_time)/1e9} s ")

        
        

var_data_centered = var_data - np.mean(var_data)

for n in tqdm(range(lags)):
        P = var_data_centered[0:num-n] * var_data_centered[n:num]
        auto[n] = np.sum(P) / (num-n)
            
timestep=time_data[1]-time_data[0]
max_lag_time=500*timestep

print(f"timestep {timestep} s, max timelag={max_lag_time} s")

amplitude = np.abs(auto)

autocorr_finished_abs_time=time.time_ns()
print(f"autocorr_time_taken={(autocorr_finished_abs_time-extracted_data_abs_time)/1e9} s ")

        # Time in milliseconds
t_ms = time_data * 1e3
index = np.where((t_ms > -0.1) & (t_ms < 4))[0]

xData = t_ms[index]
yData = auto[index]

        # Initial guesses for curve fitting
begin_fit = index[0]
end_fit = index[0] + len(index)
a0 = (auto[begin_fit] - auto[end_fit])
b0 = 1 / 6  # Initial guess for decay rate
c0 = auto[end_fit]
d0 = a0
e0 = 1 / 0.001  # Initial guess for decay rate

         #Perform curve fitting
popt, _ = curve_fit(
    fitting_func, 
    xData, 
    yData, 
    p0=[a0, b0, c0, d0, e0],
    bounds=([-np.inf, 0, -np.inf, -np.inf, 0],  # lower bounds: b > 0, e > 
            [np.inf, np.inf, np.inf, np.inf, 1/1e-6])  # upper bounds
)
a, b, c, d, e = popt

decay_ms1 = 1 / b
decay_ms2 = 1 / e


         #Store the results
results.append((run_id, decay_ms1, decay_ms2))

fit_finished_abs_time=time.time_ns()
print(f"fit_time_taken={(fit_finished_abs_time-autocorr_finished_abs_time)/1e9} s ")

        # Create the plot but don't display yet
fig, ax = plt.subplots()
ax.plot(t_ms[index], yData, "g*",label='Data',)
ax.plot(t_ms[index], fitting_func(t_ms[index], *popt), label=r'Fit $\tau_1 =$' + f'{decay_ms1} ms\n' + r'$\tau_2 =$' + f'{decay_ms2} ms')
ax.set_xlabel('Time (ms)', fontsize=14)
ax.set_ylabel('Auto Correlation', fontsize=14)
ax.legend()
ax.grid(True)
figures.append(fig)

    

# Show all plots after processing all data
for fig in figures:
    fig.tight_layout()
    plt.show()

# Print the summary table
print("Run ID   | Decay Time 1 (ms) | Decay Time 2 (ms)")
print("---------------------------------------------------")
for run_id, tau1, tau2 in results:
    print(f"{run_id:<8} | {tau1} | {tau2}")

