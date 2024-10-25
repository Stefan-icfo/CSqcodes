import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
from scipy.optimize import curve_fit
import os

# Function to define the exponential decay model (for fitting)
def exponential_fit(x, A, tau):
    return A * np.exp(-(x-3) / tau)

# Function to define the specific exponential Y = 0.002 * exp(-x / 0.007)
def custom_exponential(x):
    return 0.0002 * np.exp(-(x-0.001) / 0.007)

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

# Function to fit the exponential model to the selected region (fixed range: 2.99 to 3.03 seconds)
def fit_exponential(time_array, trace):
    # Define the x range for fitting
    x_initial, x_final = 3.0008, 3.03
    
    # Select x and y data for fitting
    mask = (time_array >= x_initial) & (time_array <= x_final)
    x_fit = time_array[mask]
    y_fit = trace[mask]

    # Estimate initial guess for parameters A and tau
    A_guess = 0.0002  # Initial value
    tau_guess = 0.007  # Rough guess for decay time

    # Perform the curve fit
    popt, _ = curve_fit(exponential_fit, x_fit, y_fit, p0=[A_guess, tau_guess])
    return x_fit, y_fit, popt

# Function to plot the fitting region and the exponential function
def plot_fitting_region_with_custom_exp(time_array, trace, popt):
    plt.figure(figsize=(10, 6))

    # Select the data for fitting between 2.99s and 3.03s
    mask = (time_array >= 3) & (time_array <= 3.03)
    x_fit = time_array[mask]
    y_fit = trace[mask]

    # Plot the original data in the selected range
    plt.plot(x_fit, y_fit, label='Selected Fit Region', color='b')

    # Plot the fitted exponential curve
    plt.plot(x_fit, exponential_fit(x_fit, *popt), 'r--', label='Exponential Fit', linewidth=2)

    # Plot the custom exponential function Y = 0.002 * exp(-t / 0.007)
    # y_custom_exp = custom_exponential(x_fit - 2.99)  # Shift time for comparison
    # plt.plot(x_fit, y_custom_exp, 'g--', label='Custom Exp: 0.002 * exp(-t / 0.007)', linewidth=2)

    plt.title('Exponential Fit and Custom Exponential Function')
    plt.xlabel('Time (s)')
    plt.ylabel('Trace (µV)')
    plt.legend()
    plt.grid()
    plt.show()

# Main function
def main():
    # Load data from a text file
    filename = r"C:\Users\LAB-nanooptomechanic\Desktop\ringdown\time_trace_datab.txt"  # Adjust as needed
    time_array, trace = load_data(filename)

    # Save the data before fitting
    save_data(time_array, trace)

    # Fit the exponential model to the data between 2.99s and 3.03s
    x_fit, y_fit, popt = fit_exponential(time_array, trace)

    # Print the fit results
    print(f"Fitted A: {popt[0]:.6f}, Decay Time (tau): {popt[1]:.6f} seconds")

    # Plot the fitting region along with the custom exponential function
    plot_fitting_region_with_custom_exp(time_array, trace, popt)

if __name__ == "__main__":
    main()
