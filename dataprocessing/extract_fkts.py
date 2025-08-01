import math
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc
import numpy as np
#enter here the database location
from database import *

dB_location=DATABASE_LOCATION#"C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD12_B4_F4v2.db'

#qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1_part3.db'

def extract_2d(run_id, 
               data_2d_name="signal_shift_Vxn_deriv",
               setpoints1_name='QDAC_ch02_dc_constant_V',  # gate 2
               setpoints2_name='QDAC_ch04_dc_constant_V',  # gate 4
               plot=True):

    #qc.config["core"]["db_location"]=dB_location
    dataset = qc.load_by_id(run_id)
    pdf_temp = dataset.to_pandas_dataframe_dict()
    data2d_raw = pdf_temp[data_2d_name]
    data2d_np = np.array(data2d_raw)


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


    data_2d = np.zeros((len(setpoints1), len(setpoints2)))
    for i in range(len(data2d_np)): 
        
        row = np.where(setpoints1 == setpoints1_np[i])[0][0]  
        col = np.where(setpoints2 == setpoints2_np[i])[0][0]  
        data_2d[row, col] = data2d_np[i]

   
    if plot:
        plt.figure(figsize=(6, 5))
        # pcolor(x-values, y-values, 2D-data)
        plt.pcolor(setpoints2, setpoints1, data_2d)
        plt.title(f"Measurement {run_id}")
        plt.colorbar(label=data_2d_name)
        plt.xlabel('Gate 4 (V)')
        plt.ylabel('Gate 2 (V)')
        plt.savefig(f"measurement_{run_id}.png", dpi=300, bbox_inches='tight')
        plt.show()

    return setpoints1, setpoints2, data_2d



def extract_1d(run_id, data_1d_name = "x", setpoint_name = 'time_param',  plot = True):


    #qc.config["core"]["db_location"]=dB_location
    experiments=qc.experiments()
    dataset=qc.load_by_id(run_id)

    interdeps = dataset.description.interdeps
    #param_spec = interdeps.non_dependencies[0]  # 
    #param_name = param_spec.name
    data_x = dataset.get_parameter_data(data_1d_name)
    setpoints_raw = data_x[data_1d_name][setpoint_name]
    setpoints_np=np.array(setpoints_raw)

    pdf_temp=dataset.to_pandas_dataframe_dict()
    data1d_raw=pdf_temp[data_1d_name]
    data1d_np=np.array(data1d_raw)
    if plot:
        plt.plot(setpoints_np,data1d_np)
        plt.title(f"measurement {run_id}")
        plt.show()
    return setpoints_np.flatten(), data1d_np.flatten()