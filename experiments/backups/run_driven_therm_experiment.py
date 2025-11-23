#from instruments import *

#measure nonlinearity 310725
import time
import copy
from experiments.zurich_experiments.spectrum_0825 import *
#Vg=0#init
exp.max_thermomech_freq=137.99e6
drive_amp=0
exp.sit_at_max_Isens()
mech_freq=exp.max_thermomech_freq
zurich.set_mixdown(mech_freq-1e6)
time.sleep(30)
#take noise once
zurich.sigout1_amp1_enabled_param.value(0) 
run_thermomech_temp_meas(exp_name='noise')
drive_amp=50e-6
i=0
while drive_amp<10e-3:
    i+=1
    #sitpos adjust
    #mech_freq=exp.max_thermomech_freq#test
    if i % 10 == 0: 
        zurich.freq2(1.25e6)
        exp.sit_at_max_Isens()
    zurich.set_mixdown(mech_freq)
    zurich.output1_amp1(drive_amp) 
    #noise
    #mech_freq=exp.max_thermomech_freq#test
    #now only once ->noise out of loop   
    #thermal
    time.sleep(30)
    run_thermomech_temp_meas(reps_nodrive=10,exp_name='thermal')
    mech_freq=exp.max_thermomech_freq#here setting for real
    print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
    zurich.set_mixdown(mech_freq)
    zurich.sigout1_amp1_enabled_param.value(1)
    #driven
    time.sleep(30)
    run_thermomech_temp_meas(exp_name=f'driven {drive_amp*1e6:4g} uV')
    #driven and swept
    #exp.mech_simple_fun_db()
    drive_amp=drive_amp*+50e-6
    print(f"setting output1 to {drive_amp*1e6:4g} uV")
    time.sleep(5)
