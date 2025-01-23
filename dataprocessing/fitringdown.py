import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
from scipy.optimize import curve_fit
from matplotlib.widgets import RectangleSelector

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Database location
qc.config["core"]["db_location"] = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_zurichdata.db"

experiments = qc.experiments()

dataset_temp = qc.load_by_id(1134)
df_temp = dataset_temp.to_pandas_dataframe_dict()

# Debug: Inspect df_temp
print("df_temp keys:", df_temp.keys())  # Print the keys to understand the structure

# Extract the appropriate DataFrame from df_temp
# Assuming you want to access the first DataFrame in the dictionary
first_key = list(df_temp.keys())[0]  # Get the first key
v_r_data = df_temp[first_key]  # Extract the DataFrame

# Debug: Check the structure of v_r_data
print("v_r_data structure:\n", v_r_data.head())  # Print the first few rows
print("Columns in v_r_data:", v_r_data.columns.tolist())  # Print all column names

# Check if 'v_r' exists and convert to numeric if it does
if 'v_r' in v_r_data.columns:
    # Ensure v_r is numeric and handle NaNs
    trace = pd.to_numeric(v_r_data["current"], errors='coerce')  # Convert to numeric, coercing errors
else:
    print("'v_r' column not found in the DataFrame. Available columns are:", v_r_data.columns.tolist())
    trace = None  # Set trace to None to avoid further errors

# Proceed only if trace has been set correctly
if trace is not None:
    trace = trace.dropna().to_numpy()  # Drop NaNs and convert to numpy array

    num_points = len(trace)  # Get the number of points in v_r data
    time_array = np.linspace(2, 5, num_points)  # Create time array from 2 to 5 seconds

    # Define an exponential decay function
    def exponential_fit(x, A, tau):
        return A * np.exp(-x / tau)

    # Set up the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(time_array, trace, marker='o', linestyle='-', color='b', label='v_r data')
    ax.set_title('Select Region for Exponential Fit')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('v_r (µV)')
    ax.grid(True)
    ax.legend()

    # Variables to hold selected data
    selected_region = None

    # Define the callback function for rectangle selection
    def onselect(eclick, erelease):
        global selected_region
        # Get coordinates of the rectangle
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        selected_region = (min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))  # (x1, y1, x2, y2)
        print(f"Selected region: {selected_region}")

    # Create a RectangleSelector
    rectangle_selector = RectangleSelector(ax, onselect, drawtype='box', useblit=True,
                                           button=[1],  # Left mouse button
                                           minspanx=5, minspany=5,
                                           spancoords='pixels', rectprops=dict(facecolor='red', alpha=0.5))

    plt.show()

    # After selecting the region, fit the exponential decay
    if selected_region is not None:
        x1, y1, x2, y2 = selected_region
        
        # Select x and y data for fitting
        mask = (time_array >= x1) & (time_array <= x2)
        x_fit = time_array[mask]
        y_fit = trace[mask]

        # Ensure y_fit is a numpy array of floats
        y_fit = np.array(y_fit, dtype=float)

        # Estimate initial guess for parameters A and tau
        A_guess = y_fit[0]  # Initial guess based on first y value in the selection
        tau_guess = (x2 - x1) / 2  # Rough guess for decay time

        # Check the dimensions of x_fit and y_fit before fitting
        if len(x_fit) > 0 and len(y_fit) > 0 and np.all(np.isfinite(y_fit)):
            # Perform the curve fit
            popt, pcov = curve_fit(exponential_fit, x_fit, y_fit, p0=[A_guess, tau_guess])
            A_fit, tau_fit = popt  # Extract fit parameters

            # Print the fit results
            print(f"Fitted A: {A_fit:.2f}, Decay Time (tau): {tau_fit:.2f} seconds")

            # Plot the fit
            plt.figure(figsize=(10, 6))
            plt.plot(time_array, trace, label='Original Data', color='b')
            plt.plot(x_fit, exponential_fit(x_fit, *popt), 'r--', label='Exponential Fit')

            # Display fit parameters on the plot
            plt.title('Exponential Fit to Selected Data')
            plt.xlabel('Time (s)')
            plt.ylabel('v_r (µV)')
            plt.legend()
            plt.grid()
            
            # Add text to display fit results
            plt.text(0.05, 0.95, f'Fitted A: {A_fit:.2f}\nDecay Time (tau): {tau_fit:.2f} s',
                     transform=plt.gca().transAxes, fontsize=12, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            plt.show()
        else:
            print("The fitting data is not valid. Please check the selected region.")
    else:
        print("No region selected for the fit.")
else:
    print("Trace could not be generated due to missing 'v_r' column.")
