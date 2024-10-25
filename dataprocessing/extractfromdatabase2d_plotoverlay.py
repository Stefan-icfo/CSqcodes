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

dataset1=qc.load_by_id(1452)#1462
dataset2=qc.load_by_id(1464)

pdf_temp1=dataset1.to_pandas_dataframe_dict()
pdf_temp2=dataset2.to_pandas_dataframe_dict()

signal_shift_Vx_deriv1_raw=pdf_temp1['signal_shift_Vx_deriv']
signal_shift_Vx_deriv2_raw=pdf_temp2['signal_shift_Vx_deriv']

signal_shift_Vx_deriv1_np=np.array(signal_shift_Vx_deriv1_raw)
signal_shift_Vx_deriv2_np=np.array(signal_shift_Vx_deriv2_raw)
# ---------------------Geting the data from the database---------------------
# pprint(dataset.get_parameter_data())
interdeps = dataset1.description.interdeps
param_spec = interdeps.non_dependencies[0]  # hall resistance data
param_name = param_spec.name

data_xy1 = dataset1.get_parameter_data(param_spec)
xy1 = data_xy1[param_name][param_name]

data_xy2 = dataset2.get_parameter_data(param_spec)
xy2 = data_xy2[param_name][param_name]

#g1:outer gate
#g2:inner gate

x1_raw = data_xy1[param_name]['QDAC_ch04_dc_constant_V']
y1_raw = data_xy1[param_name]['QDAC_ch02_dc_constant_V']
x2_raw = data_xy2[param_name]['QDAC_ch04_dc_constant_V']
y2_raw = data_xy2[param_name]['QDAC_ch02_dc_constant_V']

x1_np=np.array(x1_raw)
y1_np=np.array(y1_raw)
x2_np=np.array(x2_raw)
y2_np=np.array(y2_raw)

x1=np.unique(x1_np)
y1=np.unique(y1_np)
x2=np.unique(x2_np)
y2=np.unique(y2_np)

signal_shift_Vx_deriv1=np.zeros([len(y1), len(x1)])
signal_shift_Vx_deriv2=np.zeros([len(y2), len(x2)])


for m in range(len(y1)):
    for n in range(len(x1)):
        signal_shift_Vx_deriv1[m,n]=signal_shift_Vx_deriv1_np[m*len(x2)+n]

for m in range(len(y2)):
    for n in range(len(x2)):
        signal_shift_Vx_deriv2[m,n]=signal_shift_Vx_deriv2_np[m*len(x2)+n]

signal_shift_Vx_deriv1_truncated=np.clip(signal_shift_Vx_deriv1,-100e-6,100e-6)
signal_shift_Vx_deriv2_truncated=np.clip(signal_shift_Vx_deriv2,-100e-6,100e-6)

plt.figure(1)     
plt.pcolor(x1,y1,signal_shift_Vx_deriv1_truncated)
plt.colorbar()  
plt.show()
plt.close()

plt.figure(2) 
plt.pcolor(x2,y2,signal_shift_Vx_deriv2_truncated)
plt.colorbar()  
plt.show()
plt.close()

plt.figure(3) 
plt.pcolor(x2,y2,signal_shift_Vx_deriv1_truncated+signal_shift_Vx_deriv2_truncated)
plt.colorbar()  
plt.show()
plt.close()
#integ=np.zeros(len(vg1))

#for m in range(len(vg1)):
#    integ[m]=sum(GIV[m,:])


#substract data below a certain value
#GIV_nonoise=np.zeros([len(vg1), len(vg2)])
#integ_nonoise=np.zeros(len(vg1))

#for m in range(len(vg1)):
#    for n in range(len(vg2)):
#        if GIV[m,n]>30e-9:
#            GIV_nonoise[m,n]=GIV[m,n]

#for m in range(len(vg1)):
#    integ_nonoise[m]=sum(GIV_nonoise[m,:])
