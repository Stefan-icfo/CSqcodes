

import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import os

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Set the database location
qc.config["core"]["db_location"] = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part2.db"

# Function to process and save datasets
def process_and_save_dataset(experiment_id, output_folder):
    dataset = qc.load_by_id(experiment_id)
    df = dataset.to_pandas_dataframe_dict()

    # Extract data
    trace = np.array(df["v_r"], dtype=float).flatten()
    num_points = len(trace)

    # Generate time axis from 2 to 5 seconds
    time_array = np.linspace(2, 5, num_points)

    # Plot the data
    plt.figure(figsize=(8, 5))
   #  plt.plot(time_array, trace, label=f'Dataset {experiment_id}')
   #  plt.xlabel("Time (s)")
   #  plt.ylabel("G (a.u.)")
   #  plt.title(f'Dataset {experiment_id} - G vs Time')
    # plt.legend()
   #  plt.grid(True)
    # plt.show()

    # Save the data to a text file
    os.makedirs(output_folder, exist_ok=True)
    output_file_path = os.path.join(output_folder, f"{experiment_id}.txt")
    with open(output_file_path, 'w') as f:
        for t, g in zip(time_array, trace):
            f.write(f"{t:.6f}, {float(g):.6f}\n")
    print(f"Data extracted and saved to '{output_file_path}'")

# Main function to process multiple experiment IDs
def main():
    experiment_ids = range(2470, 2515)  # From 314 to 318
    output_folder = r"\\files\\groups\\NanoOptoMechanics\\Users\\Marta\\ringdown2"

    for experiment_id in experiment_ids:
        process_and_save_dataset(experiment_id, output_folder)

if __name__ == "__main__":
    main()
