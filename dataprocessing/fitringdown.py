import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
from scipy.optimize import curve_fit

# Database location
qc.config["core"]["db_location"] = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part3.db"

# Load dataset
dataset = qc.load_by_id(283)
# Fetch parameter data
data_dict = dataset.get_parameter_data()

# Check if 'v_r' exists
if 'v_r' in data_dict:
    v_r = data_dict['x']['x']

    # Convert to NumPy array
    v_r = np.array(v_r).flatten()

    # Get the number of points in v_r
    num_points = len(v_r)

    # Create a new time array from 2 to 5 seconds with the same number of points
    time_range = np.linspace(2, 5, num_points)

    # Step 1: Show the data between 3 and 3.05 seconds
    start_time = 2.98
    end_time = 3.05

    # Extract the segment of data between 3.00 and 3.05 seconds
    start_index = np.argmin(np.abs(time_range - start_time))
    end_index = np.argmin(np.abs(time_range - end_time))

    fit_time = time_range[start_index:end_index]
    fit_v_r = v_r[start_index:end_index]
    
    # Step 2: Display the plot for the user to visually inspect the data
    plt.figure(figsize=(10, 6))
    plt.plot(fit_time, fit_v_r, marker='o', linestyle='-', color='b', label="Data (2.98s to 3.05s)")
    plt.xlabel("Time (s)")  # Label for X-axis
    plt.ylabel("v_r (V)")   # Label for Y-axis
    plt.title("Data from 2.98s to 3.05s (Please Input Points for Fit)")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Step 3: Let the user manually input the time points for the exponential fit
    print("Please input the time points (x-values) for the start and end of the exponential fit.")
    start_fit_time = float(input("Start time (s): "))
    end_fit_time = float(input("End time (s): "))

    # Step 4: Find the closest indices in the data to the user-input points
    start_index = np.argmin(np.abs(fit_time - start_fit_time))
    end_index = np.argmin(np.abs(fit_time - end_fit_time))

    fit_time_selected = fit_time[start_index:end_index]
    fit_v_r_selected = fit_v_r[start_index:end_index]

    # Step 5: Define the exponential function for fitting
    def exp_fit(t, A, tau):
        return A * np.exp(-(t - fit_time_selected[0]) / tau)  # Shifted time values for stability

    # Step 6: Perform the fit
    # Initial guess: A is the value at the start of the fitting range, and tau is half the time range
    initial_guess = [fit_v_r_selected[0], (end_fit_time - start_fit_time) / 2]

    print(f"Initial guess: A = {initial_guess[0]:.6f}, tau = {initial_guess[1]:.6f}")

    try:
        popt, pcov = curve_fit(exp_fit, fit_time_selected, fit_v_r_selected, p0=initial_guess)

        A_opt = popt[0]  # optimal A value (amplitude)
        tau_opt = popt[1]  # optimal tau value (time constant)

        # Print the fit results
        print(f"Fit result: A = {A_opt:.6f}, tau = {tau_opt:.6f} seconds")

        # Step 7: Plot the selected data and the exponential fit
        plt.figure(figsize=(10, 6))
        plt.plot(fit_time_selected, fit_v_r_selected, marker='o', linestyle='-', color='b', label="Selected Data for Fit")
        plt.plot(fit_time_selected, exp_fit(fit_time_selected, *popt), 'r--', label=f'Exp fit: A={A_opt:.6f}, τ={tau_opt:.6f} s')
        plt.xlabel("Time (s)")  # Label for X-axis
        plt.ylabel("v_r (V)")   # Label for Y-axis
        plt.title("Exponential Fit on Selected Data")
        plt.grid(True)
        plt.legend()
        plt.show()

    except RuntimeError as e:
        print(f"Error during fitting: {e}")

else:
    print("Error: 'v_r' not found. Available keys:", data_dict.keys())

