import math
#import math
#import csv
import matplotlib.pyplot as plt
#import os


import pandas as pd
import qcodes as qc
import numpy as np
from dataprocessing.extract_fkts import *
#enter here the database location

def compare_GVgs(run_id1,run_id2):
    Vg1,G1= extract_1d(run_id1, data_1d_name = "G", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = False)
    Vg2,G2= extract_1d(run_id2, data_1d_name = "G", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = False)
    plt.plot(Vg1,G1,label=f"run id {run_id1}")
    plt.plot(Vg2,G2,label=f"run id {run_id2}")
    plt.legend()
    plt.show()

def compare_GVgs_and_sitpos_and_slope(run_id1,run_id2):
    Vg1,G1= extract_1d(run_id1, data_1d_name = "G", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = False)
    Vg2,G2= extract_1d(run_id2, data_1d_name = "G", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = False)
    V_not_used,sitpos1= extract_1d(run_id1, data_1d_name = "sitpos", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = False)
    V_not_used,sitpos2= extract_1d(run_id2, data_1d_name = "sitpos", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = False)
    V_not_used,slope1= extract_1d(run_id1, data_1d_name = "slope", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = False)
    V_not_used,slope2= extract_1d(run_id2, data_1d_name = "slope", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = False)
    
    plt.figure(1)
    plt.plot(Vg1,G1,label=f"run id {run_id1}")
    plt.plot(Vg2,G2,label=f"run id {run_id2}")
    plt.plot(Vg1,sitpos1,"o")
    plt.plot(Vg2,sitpos2,"o")
    plt.plot(Vg1,slope1)
    plt.plot(Vg2,slope2)
    plt.legend()
    plt.show()
    