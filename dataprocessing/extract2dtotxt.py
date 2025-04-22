import math
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
import numpy as np
import os

# Define the database location
dB_location = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part3.db"

# Define the output folder and file name
output_folder = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part3"
output_filename = "CD11_D7_C1_run508.txt"
output_path = os.path.join(output_folder, output_filename)

def extract_2d(run_id, 
               data_2d_name="G",
               setpoints1_name='delta',      # Y-axis: delta (mV)
               setpoints2_name='gateV',      # X-axis: gateV (V)
               plot=True,
               save_txt=True,
               txt_path=output_path,
               dB_location=dB_location):

    # Set database location
    qc.config["core"]["db_location"] = dB_location

    # Load dataset
    dataset = qc.load_by_id(run_id)
    pdf_temp = dataset.to_pandas_dataframe_dict()
    data2d_raw = pdf_temp[data_2d_name]
    data2d_np = np.array(data2d_raw)

    # Get metadata for axes
    interdeps = dataset.description.interdeps
    param_spec = interdeps.non_dependencies[0]  
    param_name = param_spec.name
    data_xy = dataset.get_parameter_data(param_spec)

    # Extract setpoints
    setpoints1_raw = data_xy[param_name][setpoints1_name]
    setpoints2_raw = data_xy[param_name][setpoints2_name]
    setpoints1_np = np.array(setpoints1_raw)
    setpoints2_np = np.array(setpoints2_raw)

    # Get unique values
    setpoints1 = np.unique(setpoints1_np)
    setpoints2 = np.unique(setpoints2_np)

    # Build 2D array
    data_2d = np.zeros((len(setpoints1), len(setpoints2)))
    for i in range(len(data2d_np)): 
        row = np.where(setpoints1 == setpoints1_np[i])[0][0]  
        col = np.where(setpoints2 == setpoints2_np[i])[0][0]  
        data_2d[row, col] = data2d_np[i]

    # Plot heatmap
    if plot:
        plt.figure(figsize=(6, 5))
        plt.pcolor(setpoints2, setpoints1, data_2d, shading='auto', cmap='plasma')
        plt.title(f"Measurement {run_id}")
        plt.colorbar(label=data_2d_name)
        plt.xlabel('Gate voltage (V)')
        plt.ylabel('Delta (mV)')
        plt.tight_layout()
        plt.savefig(f"measurement_{run_id}.png", dpi=300, bbox_inches='tight')
        plt.show()

    # Save as .txt
    if save_txt:
        with open(txt_path, 'w') as f:
            f.write("gateV(V)\tdelta(mV)\tG(μS)\n")
            for i in range(len(setpoints1)):
                for j in range(len(setpoints2)):
                    f.write(f"{setpoints2[j]:.8f}\t{setpoints1[i]:.8f}\t{data_2d[i, j]:.8f}\n")
        print(f"Data saved to: {txt_path}")

    return setpoints1, setpoints2, data_2d

# Run the function
extract_2d(run_id=508)
