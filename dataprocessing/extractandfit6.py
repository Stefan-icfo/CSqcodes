import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.optimize import curve_fit

# Function to define the exponential decay model
def exponential_fit(x, A, tau):
    return A * np.exp(-x / tau)

# Function to define the custom exponential Y = 0.002 * exp(-(t - 3) / 0.007)
def custom_exponential(x):
    return 0.002 * np.exp(-(x - 3) / 0.007)  # Decays starting from t=3

# Load data from a .txt file
def load_data(filename):
    data = np.loadtxt(filename, delimiter='\t', skiprows=1)  # Skip header
    time_array = data[:, 0]  # First column
    trace = data[:, 1]  # Second column
    return time_array, trace

# Save data to a specified output path
def save_data(time_array, trace):
    output_filename = r"C:\Users\LAB-nanooptomechanic\Desktop\ringdown\time_trace_data.txt"
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)  # Create the directory if it doesn't exist

    with open(output_filename, "w") as f:
        f.write("Time (s)\tTrace (µV)\n")  # Header
        for t, v in zip(time_array, trace.flatten()):  # Use flatten to ensure 1D
            f.write(f"{t:.6f}\t{v:.6f}\n")  # Write each pair of time and trace values
    print(f"Data saved to {output_filename}")

# Fit the exponential model to the selected region (between 3.00 and 3.03 seconds)
def fit_exponential(time_array, trace):
    # Define fitting range
    x_initial, x_final = 3.00, 3.03
    
    # Select x and y data for fitting
    mask = (time_array >= x_initial) & (time_array <= x_final)
    x_fit = time_array[mask]
    y_fit = trace[mask]

    # Estimate initial guess for parameters A and tau
    A_guess = y_fit[0]  # Initial value
    tau_guess = (x_final - x_initial) / 2  # Rough guess for decay time

    # Perform the curve fit
    popt, _ = curve_fit(exponential_fit, x_fit, y_fit, p0=[A_guess, tau_guess])
    return popt

# Function to plot the fitting region and the custom exponential function
def plot_fitting_region_with_custom_exp(time_array, trace, popt):
    plt.figure(figsize=(10, 6))

    # Full range for plotting (2.98 to 3.1 seconds)
    x_full_range = np.linspace(2.98, 3.1, 100)  
    y_full_range = trace[np.logical_and(time_array >= 2.98, time_array <= 3.1)]

    # Plot the original data in the selected range
    plt.plot(time_array, trace, label='Original Data', color='b')

    # Fit only in the defined region (3.00 to 3.03 seconds)
    plt.plot(x_full_range, exponential_fit(x_full_range, *popt), 'r--', label='Exponential Fit', linewidth=2)

    

    # Additional settings for plot
    plt.title('Exponential Fit and Custom Exponential Function')
    plt.xlabel('Time (s)')
    plt.ylabel('Trace (µV)')
    plt.xlim(2.98, 3.1)  # Limit x-axis
    plt.ylim(min(trace), max(trace))  # Limit y-axis based on original data
    plt.legend()
    plt.grid()
    plt.show()

# Main function
def main():
    # Load data from a text file
    filename = r"C:\Users\LAB-nanooptomechanic\Desktop\ringdown\time_trace_data.txt"  # Adjust as needed
    time_array, trace = load_data(filename)

    # Save the data before fitting
    save_data(time_array, trace)

    # Fit the exponential model to the specified region
    popt = fit_exponential(time_array, trace)

    # Print the fit results
    print(f"Fitted A: {popt[0]:.6f}, Decay Time (tau): {popt[1]:.6f} seconds")

    # Plot the fitting region alongside the custom exponential function
    plot_fitting_region_with_custom_exp(time_array, trace, popt)

if __name__ == "__main__":
    main()

