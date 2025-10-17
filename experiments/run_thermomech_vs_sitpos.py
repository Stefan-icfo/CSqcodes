#from instruments import *

#measure nonlinearity 310725
import time
import copy
from experiments.zurich_experiments.spectrum_0825 import *#c is 20 and false
#Vg=0#init

mech_freq=151.4639e6
#zurich.set_mixdown(mech_freq-1e6)
qdac.ch06.dc_constant_V(0.833)
#run_thermomech_temp_meas(exp_name=f'background_')
time.sleep(100)
while qdac.ch06.dc_constant_V()<0.835:
    print(f"setting zuricch to mech freq")
    zurich.set_mixdown(mech_freq)
    current_V=qdac.ch06.dc_constant_V()
    qdac.ch06.dc_constant_V(current_V+1e-4)
    run_thermomech_temp_meas(exp_name=f'thermal_gcs={current_V*1e3:6g} mV',take_time_resolved_spectrum=True)
    #mech_freq=exp.max_thermomech_freq#here setting for real
    print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
    zurich.set_mixdown(mech_freq)
    print(f"setting ch06  to {current_V:6g} mV")
    time.sleep(5)
