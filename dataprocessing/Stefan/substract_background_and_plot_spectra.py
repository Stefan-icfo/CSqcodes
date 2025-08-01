from dataprocessing.extract_fkts import *
import math
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc
import numpy as np
import copy


rbw=1.676
#runid_driven=298
runid_background=294
runids_driven=[290,298,306,314]
p_diffs=copy.copy(runids_driven)#just init

def voltage_to_psd(v_rms, rbw, impedance=50):
  
    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd


freq,v_back=extract_1d(runid_background, data_1d_name = "V_fft_avg_avg", setpoint_name = 'freq_param',  plot = False)

for runid_driven,p_diff in zip(runids_driven,p_diffs):
    freq,v_drive=extract_1d(runid_driven, data_1d_name = "V_fft_avg_avg", setpoint_name = 'freq_param',  plot = False)
    v_diff=v_drive-v_back
    p_diff=voltage_to_psd(v_diff, rbw)
    plt.plot(freq,p_diff)
plt.show()


## second part:driven swept traces

#runids_swept=[291,299,307,315]

#for runid in runids_swept:
#     freq,I_swept_avg=extract_1d(runid_driven, data_1d_name = "I_rf", setpoint_name = 'zurich_osc0_freq',  plot = False)
#     plt.plot(freq,I_swept_avg)
#plt.show()


#for p_diff in p_diffs:
#    plt.plot(freq,p_diff)
#plt.show()