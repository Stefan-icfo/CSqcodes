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
run_ids = [391]

def hann_fourier_transform(t, A, tau, offset):
    """Fourier transform of Hann window (sinc-like with sidelobes)"""
    x = np.pi * t / tau
    x = np.where(x == 0, 1e-10, x)  # Avoid division by zero
    
    # Main sinc component
    sinc_main = np.sin(x) / x
    
    # Side components (shifted)
    x_shift = np.pi * t / tau - np.pi
    x_shift = np.where(x_shift == 0, 1e-10, x_shift)
    sinc_left = np.sin(x_shift) / x_shift
    
    x_shift2 = np.pi * t / tau + np.pi  
    x_shift2 = np.where(x_shift2 == 0, 1e-10, x_shift2)
    sinc_right = np.sin(x_shift2) / x_shift2
    
    return A * (sinc_main - 0.5 * sinc_left - 0.5 * sinc_right)**2 + offset

def sinc_squared_decay(t, A, tau, offset):
    """Simple sinc-squared function"""
    x = np.pi * t / tau
    x = np.where(x == 0, 1e-10, x)  # Avoid division by zero
    return A * (np.sin(x) / x)**2 + offset

def gaussian_decay(t, A, sigma, offset):
    """Gaussian decay function"""
    return A * np.exp(-(t/sigma)**2) + offset

def exponential_double(x, a, b, c, d, e):
    """Original double exponential function"""
    return a * np.exp(-b * x) + c + d * np.exp(-e * x)

# Choose which function to use
fitting_func = hann_fourier_transform  # Try the Fourier transform of Hann

# List to store results
results = []

for run_id in run_ids:
    try:
        # Extract data for each run_id
        t_data, x_data = extract_1d(run_id, data_1d_name="x", setpoint_name='time_param', plot=False)
        print("extracted x,t data")
        
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
        
        # Filter data for fitting (just the initial clean decay)
        fit_mask = (t_ms >= 0) & (t_ms < 0.4)  # Fit only the clean initial decay
        
        # Data for fitting (just the initial clean decay)
        xData_fit = t_ms[fit_mask]
        yData_fit = auto[fit_mask]
        
        # Data for plotting (full range for visualization)
        plot_mask = (t_ms > -0.1) & (t_ms < 5)
        xData_plot = t_ms[plot_mask]
        yData_plot = auto[plot_mask]
        
        # Initial guesses based on the chosen fitting function
        begin_fit = np.where(fit_mask)[0][0]
        end_fit = np.where(fit_mask)[0][-1]
        
        if fitting_func == gaussian_decay:
            # For Gaussian: A, sigma, offset
            A0 = auto[begin_fit] - auto[end_fit]
            sigma0 = 0.3  # Characteristic width
            offset0 = 0
            p0 = [A0, sigma0, offset0]
            bounds = ([0, 0.01, -1e-10], [np.inf, 2.0, 1e-10])
            
        elif fitting_func in [hann_fourier_transform, sinc_squared_decay]:
            # For Fourier transform of Hann or sinc-squared
            A0 = auto[begin_fit] - auto[end_fit]
            tau0 = 0.3  # Characteristic time scale
            offset0 = 0
            p0 = [A0, tau0, offset0]
            bounds = ([0, 0.1, -1e-10], [np.inf, 1.0, 1e-10])
            
        else:
            # Default exponential parameters
            a0 = (auto[begin_fit] - auto[end_fit])
            b0 = 1 / 0.4
            c0 = auto[end_fit]
            d0 = a0 * 0.1
            e0 = 1 / 0.1
            p0 = [a0, b0, c0, d0, e0]
            bounds = ([1e-10, -np.inf, -np.inf, 1e-10, -np.inf], 
                     [np.inf, np.inf, np.inf, np.inf, np.inf])
        
        # Perform curve fitting on restricted data
        popt, pcov = curve_fit(fitting_func, xData_fit, yData_fit, p0=p0, bounds=bounds)
        
        # Extract parameters based on fitting function
        if fitting_func == gaussian_decay:
            A, sigma, offset = popt
            param1, param2 = sigma, 0
            param1_label, param2_label = 'σ', 'N/A'
        elif fitting_func in [hann_fourier_transform, sinc_squared_decay]:
            A, tau, offset = popt
            param1, param2 = tau, 0
            param1_label, param2_label = 'τ', 'N/A'
        else:
            # Default exponential parameters (5-parameter function)
            a, b, c, d, e = popt
            param1 = 1 / b  # decay_ms1
            param2 = 1 / e  # decay_ms2
            param1_label, param2_label = 'τ₁', 'τ₂'
        
        # Store the results
        results.append((run_id, param1, param2))
        
        # Create the plot
        plt.plot(xData_plot, yData_plot, label=f'Data_{run_id}')
        # Plot fit over full range for comparison
        plt.plot(xData_plot, fitting_func(xData_plot, *popt), 
                label=f'Fit {param1_label}={param1:0.2f} ms\n{param2_label}={param2:0.2f} ms')
        
        # Add vertical line to show fitting boundary
        plt.axvline(x=0.4, color='red', linestyle='--', alpha=0.5, 
                   label='Fit boundary (0.4 ms)')
        
        print(f"Run {run_id}: Fitted {param1_label} = {param1:.2f} ms, {param2_label} = {param2:.2f} ms")
        
    except RuntimeError as e:
        print(f"RuntimeError for run_id {run_id}: {e}. Skipping this run and continuing.")

# Show all plots after processing all data
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('Auto Correlation', fontsize=14)
plt.legend()
plt.grid(True)
plt.show()

# Print the summary table  
print("Run ID   | Parameter 1 (ms) | Parameter 2 (ms)")
print("-----------------------------------------------")
for run_id, param1, param2 in results:
    print(f"{run_id:6d} | {param1:13.2f} | {param2:13.2f}")