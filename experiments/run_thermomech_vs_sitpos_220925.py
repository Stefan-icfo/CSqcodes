#from instruments import *

#measure nonlinearity 310725
import time
import copy
from experiments.zurich_experiments.spectrum_0825 import *#c is 20 and false
#Vg=0#init
Vg,G,sens=exp.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=True)
peakpos=Vg[np.argmax(G)]
print(f"setting cs to {peakpos*1e3:.5g} mV")
approx_maxpos=peakpos
start_vg=approx_maxpos-1.5e-3
stop_vg=approx_maxpos+1.5e-3
mech_freq=146.82e6
zurich.set_mixdown(mech_freq-1e6)
qdac.ch06.dc_constant_V(start_vg)
run_thermomech_temp_meas(exp_name=f'quick_background_')
time.sleep(100)
while qdac.ch06.dc_constant_V()<stop_vg:
    zurich.set_mixdown(mech_freq)
    current_V=qdac.ch06.dc_constant_V()
    qdac.ch06.dc_constant_V(current_V+0.5e-4)
    run_thermomech_temp_meas(exp_name=f'thermalV_gcs_quick={current_V*1e3:6g} mV')
    #mech_freq=exp.max_thermomech_freq#here setting for real
    print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
    zurich.set_mixdown(mech_freq)
    print(f"setting ch06  to {current_V:6g} mV")
    time.sleep(5)
