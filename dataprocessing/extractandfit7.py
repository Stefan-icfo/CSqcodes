import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
from utils.CS_utils import centered_moving_average
import os

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Database location
qc.config["core"]["db_location"] = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part2.db"

experiments = qc.experiments()

# Function to load and save datasets
def process_and_save_dataset(dataset_id, output_folder):
    dataset_temp = qc.load_by_id(dataset_id)
    df_temp = dataset_temp.to_pandas_dataframe_dict()
    interdeps = dataset_temp.description.interdeps
    param_spec = interdeps.non_dependencies[0]
    data_x = dataset_temp.get_parameter_data(param_spec)

    # Extract trace data
    trace = np.array(df_temp["I_rf"])
    num_points = len(trace)
    time_array = np.linspace(2, 5, num_points)

    # Save data to file
    os.makedirs(output_folder, exist_ok=True)
    output_filename = os.path.join(output_folder, f"dataset_{dataset_id}_trace.txt")
    with open(output_filename, "w") as f:
        f.write("Time (s)\tTrace (µV)\n")
        for t, v in zip(time_array, trace):
            f.write(f"{t:.6f}\t{v:.6f}\n")
    print(f"Data from dataset {dataset_id} saved to {output_filename}")

    # Plot data
    plt.figure(figsize=(10, 6))
    plt.plot(time_array, trace, label=f'Dataset {dataset_id}')
    plt.xlabel('Time (s)')
    plt.ylabel('Trace (µV)')
    plt.title(f'Dataset {dataset_id} - Trace')
    plt.legend()
    plt.grid(True)
    plt.show()

# Main function to process multiple datasets
def main():
    dataset_ids = [1134, 1135, 1136, 1137]  # Add more dataset IDs as needed
    output_folder = r"\\files\\groups\\NanoOptoMechanics\\Users\\Marta\\ringdown2"

    for dataset_id in dataset_ids:
        process_and_save_dataset(dataset_id, output_folder)

if __name__ == "__main__":
    main()


