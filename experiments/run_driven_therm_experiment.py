from instruments import *

#measure nonlinearity 310725
import time
import copy
from experiments.zurich_experiments.spectrum_0825_in_construction import *
#Vg=0#init
drive_amp=50e-6
mech_freq=137.99e6
while drive_amp<10e-3:
    #sitpos adjust
    zurich.freq2(1.25e6)
    exp.sit_at_max_Isens()
    time.sleep(30)
    zurich.freq2(1.25e6)
    zurich.output1_amp1(drive_amp) 
    zurich.sigout1_amp1_enabled_param.value(0)
    zurich.set_mixdown(mech_freq-1e6)
    #noise
    run_thermomech_temp_meas(exp_name='noise')
    zurich.set_mixdown(mech_freq)
    #thermal
    run_thermomech_temp_meas(exp_name='thermal')
    mech_freq=exp.max_thermomech_freq
    print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
    zurich.set_mixdown(mech_freq)
    zurich.sigout1_amp1_enabled_param.value(1)
    #driven
    run_thermomech_temp_meas(exp_name='driven')
    #driven and swept
    exp.mech_simple_fun_db()
    
    drive_amp=drive_amp*2
    print(f"setting output1 to {drive_amp*1e6:4g} uV")
    time.sleep(5)
