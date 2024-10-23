import math
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc
import numpy as np

#database location
qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1.db'
experiments=qc.experiments()

dataset=qc.load_by_id(2975)


pdf=dataset.to_pandas_dataframe_dict()
freq_response_raw=pdf['V_rf']
freq_response_np=np.array(freq_response_raw)
# ---------------------Geting the data from the database---------------------
# pprint(dataset.get_parameter_data())
interdeps = dataset.description.interdeps
param_spec = interdeps.non_dependencies[0]  # hall resistance data
param_name = param_spec.name
data_xy = dataset.get_parameter_data(param_spec)
xy = data_xy[param_name][param_name]

#g1:outer gate
#g2:inner gate

power_raw = data_xy[param_name]['drive_mag']
freq_raw = data_xy[param_name]['zurich_oscs0_freq']
power_np=np.array(power_raw)
freq_np=np.array(freq_raw)
power_np=np.unique(power_np)
freq=np.unique(freq_np)

freq_response=np.zeros([len(power_np), len(freq)])


for m in range(len(power_np)):
    for n in range(len(freq)):
    freq_response[m,n]=freq_response_raw[m*len(freq)+n]

        #integral

integ=np.zeros(len(vg1))

for m in range(len(vg1)):
    integ[m]=sum(GIV[m,:])


#substract data below a certain value
GIV_nonoise=np.zeros([len(vg1), len(vg2)])
integ_nonoise=np.zeros(len(vg1))

for m in range(len(vg1)):
    for n in range(len(vg2)):
        if GIV[m,n]>30e-9:
            GIV_nonoise[m,n]=GIV[m,n]

for m in range(len(vg1)):
    integ_nonoise[m]=sum(GIV_nonoise[m,:])
