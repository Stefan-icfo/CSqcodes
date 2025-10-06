#from instruments import *

#measure nonlinearity 310725
import time
import copy
from experiments.zurich_experiments.spectrum_0825 import *#c is 20 and false
from cs_experiment import *
from instruments import *
#Vg=0#init
#zurich.set_mixdown(mech_freq-1e6)
qdac.ch06.dc_constant_V(0.9557)
#run_thermomech_temp_meas(exp_name=f'background_')
time.sleep(100)
while qdac.ch06.dc_constant_V()<0.9585:
    current_V=qdac.ch06.dc_constant_V()
    qdac.ch06.dc_constant_V(current_V+0.1e-3)
    exp.mech_simple_fun_db(start_f=136.65e6,stop_f=136.85e6,step_num_f=200*5,costum_prefix=f"driven_vs_sitpos_{current_V*1e3:6g}mV_g3")
    exp.mech_simple_fun_db(start_f=136.85e6,stop_f=136.65e6,step_num_f=200*5,costum_prefix=f"driven_vs_sitpos_{current_V*1e3:6g}mV_g3")
    #mech_freq=exp.max_thermomech_freq#here setting for rea
    time.sleep(5)
