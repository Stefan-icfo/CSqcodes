import math
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc
import numpy as np
#enter here the database location
#qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1.db'

#C:\Users\sforstner\Desktop\Triton database
#"data_2d_name = "psd""
#setpoints1_name = 'freq_param'
#setpoints2_name ='time_param'


def extract_2d(run_id, data_2d_name = "I_rf",setpoints1_name = 'delta', setpoints2_name = 'zurich_oscs0_freq' , plot = True):
#database location
    qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"sforstner"+"\\"+"Desktop"+"\\"+"Triton database"+"\\"+'CD11_D7_C1.db'
    experiments=qc.experiments()
    dataset=qc.load_by_id(run_id)


    pdf_temp=dataset.to_pandas_dataframe_dict()
    data2d_raw=pdf_temp[data_2d_name]
    data2d_np=np.array(data2d_raw)
    # ---------------------Geting the data from the database---------------------
    # pprint(dataset.get_parameter_data())
    found_dependency=False
    interdeps = dataset.description.interdeps
    for param_spec in interdeps.non_dependencies:
        if param_spec.name==data_2d_name:
            found_dependency=True
            break
    if found_dependency==False:
        raise Exception("data_2d_name doesnt seem to exist")
    param_name = param_spec.name

    data_xy = dataset.get_parameter_data(param_spec)
    xy = data_xy[param_name][param_name]

    #g1:outer gate
    #g2:inner gate

    setpoints1_raw = data_xy[param_name][setpoints1_name]
    setpoints2_raw = data_xy[param_name][setpoints2_name]
    setpoints1_np=np.array(setpoints1_raw)
    setpoints2_np=np.array(setpoints2_raw)
    setpoints1=np.unique(setpoints1_np)
    setpoints2=np.unique(setpoints2_np)

    data_2d = np.zeros([len(setpoints2), len(setpoints1)])  

    for n in range(len(setpoints2)):  # Loop over the rows (y-axis)
        for m in range(len(setpoints1)):  # Loop over the columns (x-axis)
        
            data_2d[n, m] = data2d_np[n * len(setpoints1) + m]  


        
    if plot:
        plt.pcolor(setpoints1,setpoints2,data_2d)
        plt.title(f"measurement {run_id}")
        plt.show()
    return setpoints1, setpoints2, data_2d


def extract_1d(run_id, data_1d_name = "G", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = True):
    experiments=qc.experiments()
    dataset=qc.load_by_id(run_id)

    interdeps = dataset.description.interdeps
    param_spec = interdeps.non_dependencies[0]  # 
    #param_name = param_spec.name
    data_x = dataset.get_parameter_data(param_spec)
    setpoints_raw = data_x[data_1d_name][setpoint_name]
    setpoints_np=np.array(setpoints_raw)

    pdf_temp=dataset.to_pandas_dataframe_dict()
    data1d_raw=pdf_temp[data_1d_name]
    data1d_np=np.array(data1d_raw)
    if plot:
        plt.plot(setpoints_np,data1d_np)
        plt.title(f"measurement {run_id}")
        plt.show()
    return setpoints_np, data1d_np