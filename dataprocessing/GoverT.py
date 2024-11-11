import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import Boltzmann as kB_eV  # Boltzmann constant in eV/K
from dataprocessing.extract_fkts import extract_1d  # Assuming you have this import available

# Define the thermal_CB_peak function
def thermal_CB_peak(x, peak_V, G_infty, T, Delta_E_V=0.091, alpha=0.15, offset=0):
    Delta_E = alpha * Delta_E_V
    delta = alpha * (x - peak_V)
    term = (Delta_E / (4 * kB_eV * T)) * np.cosh(delta / (2 * kB_eV * T)) ** -2
    G = G_infty * term + offset
    return G

# Define temperature sets and corresponding file numbers
temperature_data = {
    50: [num for num in range(3001, 3005)],
    100: [num for num in range(2982, 2988) if num not in [2983, 2984]],
    200: [num for num in range(2988, 2993) if num not in [2989]],
    300: [num for num in range(3011, 3018) if num not in [3016, 3015, 3017]],
    400: [num for num in range(2996, 3001) if num not in [2999]],
    500: [num for num in range(3022, 3028) if num not in [3023, 3027]],
    600: [num for num in range(3030, 3035) if num not in [3031]]
}

# Data storage for table and plots
results = []
fit_examples = {}

# Process each temperature set
for temperature, file_numbers in temperature_data.items():
    G_max_list = []
    FWHM_list = []
    
    # Process each file in the temperature set
    for file_number in file_numbers:
        # Extract the data for the specified run ID
        gate_sweep, Glist = extract_1d(file_number)
        
        # Extract only the first trace and flatten it
        voltage_data = gate_sweep  # Assuming gate_sweep is the x-axis (voltage)
        conductance_data = Glist[0].flatten()  # Only take the first trace and flatten it

        # Use curve fitting to fit the data to the thermal_CB_peak function
        try:
            # Improved initial parameter guesses
            initial_peak_V = voltage_data[np.argmax(conductance_data)]
            initial_G_infty = np.max(conductance_data)
            initial_T = temperature / 1000

            popt, _ = curve_fit(
                thermal_CB_peak,
                voltage_data,
                conductance_data,
                p0=[initial_peak_V, initial_G_infty, initial_T],
                bounds=([initial_peak_V - 0.01, 0, 0], [initial_peak_V + 0.01, initial_G_infty * 2, initial_T * 2])
            )
            
            peak_V, G_infty, T_fit = popt[:3]
            
            # Calculate G_max and FWHM for this peak
            G_max = thermal_CB_peak(peak_V, *popt)
            half_max = G_max / 2
            FWHM_points = np.where(conductance_data >= half_max)[0]
            FWHM = voltage_data[FWHM_points[-1]] - voltage_data[FWHM_points[0]]
            
            # Store the results for this file
            G_max_list.append(G_max)
            FWHM_list.append(FWHM)
            
            # Save the first fit example for each temperature
            if temperature not in fit_examples:
                fit_examples[temperature] = (voltage_data, conductance_data, popt)
            
        except RuntimeError:
            print(f"Fit did not converge for file meas{file_number}")

    # Compute the averages for this temperature
    G_max_avg = np.mean(G_max_list) if G_max_list else None
    FWHM_avg = np.mean(FWHM_list) if FWHM_list else None

    # Append results to the final data
    results.append({
        "Temperature (mK)": temperature,
        "G_max_avg": G_max_avg,
        "FWHM_avg": FWHM_avg
    })

# Convert results to DataFrame for tabular display
df_results = pd.DataFrame(results)

# Display the table with Temperature, G, and FWHM
display_table = df_results[["Temperature (mK)", "G_max_avg", "FWHM_avg"]]
print(display_table)

# Plot G_max vs Temperature
plt.figure()
plt.plot(df_results["Temperature (mK)"], df_results["G_max_avg"], marker='o', linestyle='-')
plt.xlabel("Temperature (mK)")
plt.ylabel("Average G_max")
plt.title("Average G_max vs Temperature")
plt.grid(True)
plt.show()

# Plot FWHM vs Temperature
plt.figure()
plt.plot(df_results["Temperature (mK)"], df_results["FWHM_avg"], marker='o', linestyle='-')
plt.xlabel("Temperature (mK)")
plt.ylabel("Average FWHM")
plt.title("Average FWHM vs Temperature")
plt.grid(True)
plt.show()

# Plot fit examples for each temperature set
for temperature, (voltage_data, conductance_data, popt) in fit_examples.items():
    plt.figure()
    plt.plot(voltage_data, conductance_data, 'o', label="Data")
    plt.plot(voltage_data, thermal_CB_peak(voltage_data, *popt), '-', label="Fit")
    plt.xlabel("Voltage (V)")
    plt.ylabel("Conductance (G)")
    plt.title(f"Conductance Fit for {temperature} mK")
    plt.legend()
    plt.grid(True)
    plt.show()

