import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
from scipy.optimize import curve_fit
import os

# Configure QCoDeS database location
qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part2.db"
)

# -------------------- Utility Functions --------------------

# Function to load data from .txt file
def load_data(filename):
    data = np.loadtxt(filename, delimiter='\t', skiprows=1)  # Skip header
    time_array = data[:, 0]  # First column: time
    trace = data[:, 1]       # Second column: trace
    return time_array, trace

# Function to save time and trace data to .txt file
def save_data(time_array, trace):
    output_filename = r"C:\Users\LAB-nanooptomechanic\Desktop\ringdown1\1frequency.txt"
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    with open(output_filename, "w") as f:
        f.write("Time (s)\tTrace (µV)\n")
        for t, v in zip(time_array.flatten(), trace.flatten()):
            f.write(f"{t:.6f}\t{v:.6f}\n")
    print(f"Data saved to {output_filename}")

# Exponential decay model for fitting
def exponential_fit(x, A, tau):
    return A * np.exp(-(x - 3) / tau)

# Fit exponential in a fixed time range
def fit_exponential(time_array, trace):
    x_initial, x_final = 2.998, 3.013
    mask = (time_array >= x_initial) & (time_array <= x_final)
    x_fit = time_array[mask]
    y_fit = trace[mask]

    A_guess = 0.0001
    tau_guess = 0.0009
    popt, _ = curve_fit(exponential_fit, x_fit, y_fit, p0=[A_guess, tau_guess])
    return x_fit, y_fit, popt

# Plot the original data and fit
def plot_fitting_region_with_custom_exp(time_array

