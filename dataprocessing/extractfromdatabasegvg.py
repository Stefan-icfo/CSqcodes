import matplotlib.pyplot as plt 
import qcodes as qc
import numpy as np
import pandas as pd
import time
# Database location
qc.config["core"]["db_location"] = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1.db"
experiments = qc.experiments()

# Load the dataset by its ID
dataset_id = 2052  # Update this to your actual dataset ID
dataset = qc.load_by_id(dataset_id)

# Convert dataset to pandas DataFrame
pdf_2052 = dataset.to_pandas_dataframe_dict()
G_raw = pdf_2052["G"]  # Extract the G values

# Step 1: Create a new voltage array
start_vg = -1.234
stop_vg = -1.288
num_points = len(G_raw)  # Use the same number of points as G_raw

# Create an array of voltages linearly spaced
voltage_array = np.linspace(start_vg, stop_vg, num_points)

# Step 2: Plot G_raw against the new voltage array
plt.figure(figsize=(10, 6))
plt.plot(voltage_array, G_raw, marker='o', linestyle='-', color='b', label='G_raw vs Voltage')
plt.title('Conductance (G_raw) vs Voltage')
plt.xlabel('Voltage (V)')
plt.ylabel('Conductance (S)')
plt.grid()
plt.legend()
plt.show()

# Step 3: Save the data to a CSV file
data_to_save = pd.DataFrame({'Voltage (V)': voltage_array, 'G_raw (S)': G_raw})
csv_file_path = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\G_raw_vs_Voltage.csv"  # Change this path if needed
data_to_save.to_csv(csv_file_path, index=False)

print(f'Data saved to {csv_file_path}')

