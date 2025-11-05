import math
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
import numpy as np
from database import *
from matplotlib.colors import LogNorm
import os
from datetime import datetime

#dB_location = DATABASE_LOCATION

def extract_2d(run_id,
               data_2d_name="signal_shift_Vxn_deriv",
               setpoints1_name='QDAC_ch02_dc_constant_V',  # gate 2
               setpoints2_name='QDAC_ch04_dc_constant_V',  # gate 4
               plot=True, log=False, progress_report=False,
               save_data=True, output_dir="extracted_data"):
    
    dataset = qc.load_by_id(run_id)
    pdf_temp = dataset.to_pandas_dataframe_dict()
    data2d_raw = pdf_temp[data_2d_name]
    data2d_np = np.array(data2d_raw)
    
    if progress_report:
        print("loaded dataset")
    
    interdeps = dataset.description.interdeps
    param_spec = interdeps.non_dependencies[0]  
    param_name = param_spec.name
    data_xy = dataset.get_parameter_data(param_spec)
    setpoints1_raw = data_xy[param_name][setpoints1_name]
    setpoints2_raw = data_xy[param_name][setpoints2_name]
    setpoints1_np = np.array(setpoints1_raw)
    setpoints2_np = np.array(setpoints2_raw)
   
    setpoints1 = np.unique(setpoints1_np)
    setpoints2 = np.unique(setpoints2_np)
    
    if progress_report:
        print("formatted setpoints")
    
    data_2d = np.zeros((len(setpoints2), len(setpoints1)))  # Note: swapped dimensions
    
    for i in range(len(data2d_np)):
        row = np.where(setpoints1 == setpoints1_np[i])[0][0]  
        col = np.where(setpoints2 == setpoints2_np[i])[0][0]  
        data_2d[col, row] = data2d_np[i]  # Note: swapped indices
    
    if progress_report:
        print("formatted 2ddata")
    
    # Save data to files
    if save_data:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_filename = f"{output_dir}/run{run_id}_2d_{timestamp}"
        
        # Save as numpy arrays
        np.save(f"{base_filename}_data.npy", data_2d)
        np.save(f"{base_filename}_setpoints1.npy", setpoints1)
        np.save(f"{base_filename}_setpoints2.npy", setpoints2)
        
        # Save as text files
        np.savetxt(f"{base_filename}_data.txt", data_2d, delimiter='\t',
                   header=f"2D data from run {run_id}\nData: {data_2d_name}\nRows: {setpoints2_name}\nCols: {setpoints1_name}")
        np.savetxt(f"{base_filename}_setpoints1.txt", setpoints1, delimiter='\t',
                   header=f"{setpoints1_name} values")
        np.savetxt(f"{base_filename}_setpoints2.txt", setpoints2, delimiter='\t',
                   header=f"{setpoints2_name} values")
        
        # Save metadata
        with open(f"{base_filename}_metadata.txt", 'w') as f:
            f.write(f"Run ID: {run_id}\n")
            f.write(f"Experiment: {dataset.exp_name}\n")
            f.write(f"Data parameter: {data_2d_name}\n")
            f.write(f"Setpoints 1: {setpoints1_name}\n")
            f.write(f"Setpoints 2: {setpoints2_name}\n")
            f.write(f"Data shape: {data_2d.shape}\n")
            f.write(f"Setpoints 1 range: {setpoints1[0]:.6f} to {setpoints1[-1]:.6f}\n")
            f.write(f"Setpoints 2 range: {setpoints2[0]:.6f} to {setpoints2[-1]:.6f}\n")
            f.write(f"Extraction timestamp: {timestamp}\n")
        
        if progress_report:
            print(f"Data saved to {base_filename}_*")
    
    if plot:
        plt.figure(figsize=(6, 5))
        # pcolormesh(x-values, y-values, 2D-data)
        if log:
            plt.pcolormesh(setpoints1, setpoints2, data_2d, norm=LogNorm())
        else:
            plt.pcolormesh(setpoints1, setpoints2, data_2d)
        
        plt.title(f"Measurement {run_id}")
        plt.colorbar(label=data_2d_name)
        plt.xlabel(setpoints1_name)
        plt.ylabel(setpoints2_name)
        plt.savefig(f"measurement_{run_id}.png", dpi=300, bbox_inches='tight')
        plt.show()
    
    return setpoints1, setpoints2, data_2d


def extract_1d(run_id, data_1d_name="x", setpoint_name='time_param', 
               plot=True, return_exp_name=False,
               save_data=True, output_dir="extracted_data"):

    #qc.config["core"]["db_location"]=dB_location
    experiments = qc.experiments()
    dataset = qc.load_by_id(run_id)

    interdeps = dataset.description.interdeps
    #param_spec = interdeps.non_dependencies[0]  # 
    #param_name = param_spec.name
    data_x = dataset.get_parameter_data(data_1d_name)
    setpoints_raw = data_x[data_1d_name][setpoint_name]
    setpoints_np = np.array(setpoints_raw)

    pdf_temp = dataset.to_pandas_dataframe_dict()
    data1d_raw = pdf_temp[data_1d_name]
    data1d_np = np.array(data1d_raw)
    
    # Flatten the arrays
    setpoints_flat = setpoints_np.flatten()
    data1d_flat = data1d_np.flatten()
    
    # Save data to files
    if save_data:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_filename = f"{output_dir}/run{run_id}_1d_{timestamp}"
        
        # Save as numpy arrays
        np.save(f"{base_filename}_data.npy", data1d_flat)
        np.save(f"{base_filename}_setpoints.npy", setpoints_flat)
        
        # Save as text files (two-column format for easy plotting)
        combined_data = np.column_stack((setpoints_flat, data1d_flat))
        np.savetxt(f"{base_filename}_combined.txt", combined_data, delimiter='\t',
                   header=f"Run {run_id}: {setpoint_name}\t{data_1d_name}")
        
        # Save metadata
        with open(f"{base_filename}_metadata.txt", 'w') as f:
            f.write(f"Run ID: {run_id}\n")
            f.write(f"Experiment: {dataset.exp_name}\n")
            f.write(f"Data parameter: {data_1d_name}\n")
            f.write(f"Setpoint parameter: {setpoint_name}\n")
            f.write(f"Number of points: {len(data1d_flat)}\n")
            f.write(f"Setpoint range: {setpoints_flat[0]:.6f} to {setpoints_flat[-1]:.6f}\n")
            f.write(f"Data range: {np.min(data1d_flat):.6f} to {np.max(data1d_flat):.6f}\n")
            f.write(f"Extraction timestamp: {timestamp}\n")
        
        print(f"Data saved to {base_filename}_*")
    
    if plot:
        plt.plot(setpoints_flat, data1d_flat)
        plt.title(f"measurement {run_id}")
        plt.xlabel(setpoint_name)
        plt.ylabel(data_1d_name)
        plt.show()
    
    if return_exp_name == False:
        return setpoints_flat, data1d_flat
    else:
        return dataset.exp_name, setpoints_flat, data1d_flat