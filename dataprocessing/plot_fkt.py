#plotting from database

import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *

def compare_GVGs(run_id1,run_id2):
    v1,G1=extract_1d(run_id1,plot=False)
    v2,G2=extract_1d(run_id2,plot=False)
    plt.figure(1)
    plt.title("comparingGVGs")
    plt.ylabel("Conductance [S]")
    plt.xlabel("gate voltage [V]")
    plt.plot(v1,G1, label=f"run_id {run_id1}")
    plt.plot(v2,G2, label=f"run_id {run_id2}")
    plt.legend()
    plt.show()

