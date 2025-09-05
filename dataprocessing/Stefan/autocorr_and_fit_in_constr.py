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

qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v3.db'
#run_ids = [380,381,382,383,389,391]#
run_ids = [391]
qc.config["core"]["db_location"]='.\\Data\\Raw_data\\CD12_B5_F4v7.db'
#run_ids = list(range(45, 94))#will take very long time
run_ids = list(range(45, 50))#shortened
#run_ids=[48]
#filter_run_id=[389]

# Define your fitting function
def fitting_func(x, a, b, c, d, e):
    return a * np.exp(-b * x) + c + d * np.exp(-e * x)

# List to store results
results = []
# List to store figures for later plotting
figures = []

for run_id in run_ids:
    try:
        # Extract data for each run_id
        t_data, x_data = extract_1d(run_id, data_1d_name="x", setpoint_name='time_param', plot=False)
        print("extracted x,t data")
        #_, y_data = extract_1d(run_id, data_1d_name="v_r", setpoint_name='time_param', plot=False)
        #print("extracted y data")
        #t_data, phase_data = extract_1d(run_id, data_1d_name="Phase", setpoint_name='time_param', plot=False)
        #print("extracted phase data")
        
        # Choose phase_data as the variable to analyze
        var_data = x_data
        time_data = t_data - t_data[0]
        
        # Downsample data
        frac = 1
        time_data = time_data[0:len(time_data)//frac]
        var_data = var_data[0:len(var_data)//frac]
        
        # Calculate autocorrelation function
        num = len(time_data)
        auto = np.zeros(num)
        for n in tqdm(range(500)):
            P = var_data[0:num-n] * var_data[n:num]  # This ensures both arrays have the same length
            auto[n] = np.sum(P) / (num-n)  # Normalize the product
            
        amplitude = np.abs(auto)
        
        # Time in milliseconds
        t_ms = time_data * 1e3
        
        # MODIFIED: Filter data for fitting (only < 1 ms)
        fit_mask = (t_ms >= 0) & (t_ms < 1.0)  # Only fit data < 1 ms
        
        # Data for fitting (restricted range)
        xData_fit = t_ms[fit_mask]
        yData_fit = auto[fit_mask]
        
        # Data for plotting (full range for visualization)
        plot_mask = (t_ms > -0.1) & (t_ms < 5)
        xData_plot = t_ms[plot_mask]
        yData_plot = auto[plot_mask]
        
        # Initial guesses for curve fitting (based on fitting range)
        begin_fit = np.where(fit_mask)[0][0]
        end_fit = np.where(fit_mask)[0][-1]
        
        a0 = (auto[begin_fit] - auto[end_fit])
        b0 = 1 / 0.4  # Initial guess for decay rate (tau ~ 0.4 ms)
        c0 = auto[end_fit]
        d0 = a0 * 0.1  # Smaller second component
        e0 = 1 / 0.1   # Faster second decay
        
        # Perform curve fitting on restricted data
        lower_bounds = [1e-10, -np.inf, -np.inf, 1e-10, -np.inf]  # a > 0, d > 0
        upper_bounds = [np.inf, np.inf, np.inf, np.inf, np.inf]    # no upper limits
        
        popt, pcov = curve_fit(fitting_func, xData_fit, yData_fit,
                       p0=[a0, b0, c0, d0, e0],
                       bounds=(lower_bounds, upper_bounds))
        
        a, b, c, d, e = popt
        decay_ms1 = 1 / b
        decay_ms2 = 1 / e
        
        # Store the results
        results.append((run_id, decay_ms1, decay_ms2))
        
        # Create the plot
        plt.plot(xData_plot, yData_plot, label=f'Data_{run_id}')
        # Plot fit over full range for comparison, but fitted only to < 1ms data
        plt.plot(xData_plot, fitting_func(xData_plot, *popt), 
                label=r'Fit $\tau_1 =$' + f'{decay_ms1:0.2f} ms\n' + r'$\tau_2 =$' + f'{decay_ms2:0.2f} ms')
        
        # Add vertical line to show fitting boundary
        plt.axvline(x=1.0, color='red', linestyle='--', alpha=0.5, 
                   label='Fit boundary (1 ms)')
        
        print(f"Run {run_id}: Fitted τ₁ = {decay_ms1:.2f} ms, τ₂ = {decay_ms2:.2f} ms")
        
    except RuntimeError as e:
        print(f"RuntimeError for run_id {run_id}: {e}. Skipping this run and continuing.")

# Show all plots after processing all data
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('Auto Correlation', fontsize=14)
plt.legend()
plt.grid(True)
plt.show()

# Print the summary table
print("Run ID   | Decay Time 1 (ms) | Decay Time 2 (ms)")
print("---------------------------------------------------")
for run_id, tau1, tau2 in results:
    print(f"{run_id:6d} | {tau1:13.2f} | {tau2:13.2f}")