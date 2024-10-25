import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
from scipy.optimize import curve_fit
import os

# Function to define the exponential decay model
def exponential_fit(x, A, tau):
    return A * np.exp(-x / tau)

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

# Allow user to select two points for fitting
def select_region(time_array, trace):
    plt.figure(figsize=(10, 6))
    plt.plot(time_array, trace, marker='o', linestyle='-', color='b', label='Data')
    plt.title('Select Points for Exponential Fit')
    plt.xlabel('Time (s)')
    plt.ylabel('Trace (µV)')
    plt.grid(True)
    plt.legend()

    # List to store selected points
    selected_points = []

    # Define a function to handle mouse click events
    def onclick(event):
        if event.xdata is not None and event.ydata is not None:
            # Store the selected point
            selected_points.append((event.xdata, event.ydata))
            plt.plot(event.xdata, event.ydata, 'ro')  # Mark selected point
            plt.draw()

            # If two points are selected, draw a line and close the plot
            if len(selected_points) == 2:
                # Draw a line between the selected points
                x_vals = [selected_points[0][0], selected_points[1][0]]
                y_vals = [selected_points[0][1], selected_points[1][1]]
                plt.plot(x_vals, y_vals, 'g--')  # Draw a dashed green line
                plt.text((x_vals[0] + x_vals[1]) / 2, (y_vals[0] + y_vals[1]) / 2,
                         'Selected Region', color='green', fontsize=10, ha='center')
                plt.draw()
                plt.close()  # Close the plot after two points have been selected

    # Connect the click event to the onclick function
    cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    # Check if exactly two points were selected
    if len(selected_points) != 2:
        raise ValueError("Please select exactly two points for the fit.")

    return selected_points

# Fit the exponential model to the selected region
def fit_exponential(time_array, trace, selected_points):
    x_initial, y_initial = selected_points[0]
    x_final, y_final = selected_points[1]

    # Select x and y data for fitting
    mask = (time_array >= min(x_initial, x_final)) & (time_array <= max(x_initial, x_final))
    x_fit = time_array[mask]
    y_fit = trace[mask]

    # Estimate initial guess for parameters A and tau
    A_guess = y_initial  # The initial value at x_initial
    tau_guess = (x_final - x_initial) / 2  # Rough guess for decay time

    # Perform the curve fit
    popt, pcov = curve_fit(exponential_fit, x_fit, y_fit, p0=[A_guess, tau_guess])
    return popt

# Function to plot the fitting region separately
def plot_fitting_region(time_array, trace, selected_points, popt):
    plt.figure(figsize=(10, 6))
    x_initial, y_initial = selected_points[0]
    x_final, y_final = selected_points[1]

    # Select x and y data for fitting
    mask = (time_array >= min(x_initial, x_final)) & (time_array <= max(x_initial, x_final))
    x_fit = time_array[mask]
    y_fit = trace[mask]

    plt.plot(x_fit, y_fit, label='Selected Fit Region', color='b')
    plt.plot(x_fit, exponential_fit(x_fit, *popt), 'r--', label='Exponential Fit', linewidth=2)
    plt.title('Exponential Fit in Selected Region')
    plt.xlabel('Time (s)')
    plt.ylabel('Trace (µV)')
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

    # Select the fitting region
    selected_points = select_region(time_array, trace)

    # Fit the exponential model to the selected region
    A_fit, tau_fit = fit_exponential(time_array, trace, selected_points)

    # Print the fit results
    print(f"Fitted A: {A_fit:.2f}, Decay Time (tau): {tau_fit:.2f} seconds")

    # Plot the original data with the fit
    plt.figure(figsize=(10, 6))
    plt.plot(time_array, trace, label='Original Data', color='b')
    plt.plot(time_array, exponential_fit(time_array, *[A_fit, tau_fit]), 'r--', label='Exponential Fit')
    plt.title('Exponential Fit to Selected Data')
    plt.xlabel('Time (s)')
    plt.ylabel('Trace (µV)')
    plt.legend()
    plt.grid()
    plt.show()

    # Plot the fitting region separately
    plot_fitting_region(time_array, trace, selected_points, [A_fit, tau_fit])

if __name__ == "__main__":
    main()
