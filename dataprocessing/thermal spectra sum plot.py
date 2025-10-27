import qcodes as qc
import numpy as nps
import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *

runids=[1123]
hole_nr=[34,33,32,31,30,29,28,27,26,25,24]

ps_list=[]
for runid in runids:
    freq,psd=extract_1d(data_1d_name='avg_avg_nodrive_avg_substracted',setpoint_name='freq_param',plot=True)
    ps_list.append(np.sum(np.array(psd)))

plt.plot(hole_nr,ps_list)
plt.show()