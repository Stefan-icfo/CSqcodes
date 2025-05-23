# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 16:19:10 2024

@author: rtormo
"""

from tqdm import tqdm
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, hbar, e, h
from scipy.constants import Boltzmann as kB
from scipy.optimize import curve_fit
from dataprocessing.extract_fkts import *

# List of run_id's to analyze
run_ids = [240,242]

# Define your fitting function
def fitting_func(x, a, b, c, d, e):
    return a * np.exp(-b * x) + c + d * np.exp(-e * x)

# List to store results
results = []

for run_id in run_ids:
    # Extract data for each run_id
    t_data, x_data = extract_1d(run_id, data_1d_name="x", setpoint_name='time_param', plot=False)
    _, y_data = extract_1d(run_id, data_1d_name="v_r", setpoint_name='time_param', plot=False)
    _, phase_data = extract_1d(run_id, data_1d_name="Phase", setpoint_name='time_param', plot=False)

    # Choose phase_data as the variable to analyze
    var_data = phase_data
    time_data = t_data - t_data[0]

    # Downsample data
    frac = 1
    time_data = time_data[0:len(time_data)//frac]
    var_data = var_data[0:len(var_data)//frac]

    # Calculate autocorrelation function
    num = len(time_data)
    auto = np.zeros(num)

    for n in tqdm(range(num)):
        P = var_data[0:num-n] * var_data[n:num]  # This ensures both arrays have the same length
        auto[n] = np.sum(P) / (num-n)  # Normalize the product

    amplitude = np.abs(auto)

    # Time in milliseconds
    t_ms = time_data * 1e3
    index = np.where((t_ms > -0.1) & (t_ms < 12))[0]

    xData = t_ms[index]
    yData = auto[index]

    # Initial guesses for curve fitting
    begin_fit = index[0]
    end_fit = index[0] + len(index)
    a0 = (auto[begin_fit] - auto[end_fit])
    b0 = 1 / 6  # Initial guess for decay rate
    c0 = auto[end_fit]
    d0 = a0
    e0 = 1 / 0.02  # Initial guess for decay rate

    # Perform curve fitting
    popt, _ = curve_fit(fitting_func, xData, yData, p0=[a0, b0, c0, d0, e0])
    a, b, c, d, e = popt

    decay_ms1 = 1 / b
    decay_ms2 = 1 / e

    # Store the results
    results.append((run_id, decay_ms1, decay_ms2))

    # Plot the autocorrelation and fitted curve
    plt.figure()
    plt.plot(t_ms[index], yData, label='Data')
    plt.plot(t_ms[index], fitting_func(t_ms[index], *popt), label=r'Fit $\tau_1 =$' + f'{decay_ms1:0.2f} ms\n' + r'$\tau_2 =$' + f'{decay_ms2:0.2f} ms')
    plt.xlabel('Time (ms)', fontsize=14)
    plt.ylabel('Auto Correlation', fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# Print the summary table
print("Run ID   | Decay Time 1 (ms) | Decay Time 2 (ms)")
print("---------------------------------------------------")
for run_id, tau1, tau2 in results:
    print(f"{run_id:<8} | {tau1:<18.2f} | {tau2:<18.2f}")
