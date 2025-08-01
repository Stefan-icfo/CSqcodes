from instruments import *
import copy
"""
#frequency sweep
start_f=100e6
stop_f=200e6
while stop_f<700e6:
    exp.do_GVg_and_adjust_sitpos()
    zurich.sigout1_amp1_enabled_param.value(1)
    exp.mech_simple_fun_db(start_f=start_f,stop_f=stop_f)
    zurich.sigout1_amp1_enabled_param.value(0)
    start_f+=100e6
    if stop_f<600e6:
        stop_f+=100e6 
    else:
       stop_f=700e6

       

#extra repeats

exp.do_GVg_and_adjust_sitpos()
zurich.sigout1_amp1_enabled_param.value(1)
exp.mech_simple_fun_db(start_f=600e6,stop_f=700e6)#not measured yet
zurich.sigout1_amp1_enabled_param.value(0)

#for tomorrow
#
qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[0,-3,0.6,0.6,0.6,-2],step_size=0.1e-1,ramp_speed=1e-3)
exp.GVG_fun()
qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[-2,-3,0.6,0.6,0.6,-0.5],step_size=0.1e-1,ramp_speed=1e-3)

qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[-0.5,-3,0.6,0.6,0.6,-0.5],step_size=0.1e-1,ramp_speed=1e-3)
exp.linesweep_parallel(costum_prefix='g1_WIDEEEE')
#peak, mechanics,...





#power sweep

start_p=20e-3
min_power=0.4e-6
division=2
power=copy.copy(start_p)
while power>min_power:
    zurich.output1_amp1(power)
    exp.do_GVg_and_adjust_sitpos()
    zurich.sigout1_amp1_enabled_param.value(1)
    exp.mech_simple_fun_db(start_f=133e6,stop_f=136e6,step_num_f=300*10)
    exp.mech_simple_fun_db(start_f=149.7e6,stop_f=150.2e6,step_num_f=50*50)
    exp.mech_simple_fun_db(start_f=295.7e6,stop_f=296.2e6,step_num_f=50*50)
    exp.mech_simple_fun_db(start_f=314e6,stop_f=316e6,step_num_f=200*10)
    zurich.sigout1_amp1_enabled_param.value(0)
    power=power/2


qdac.ch06.dc_constant_V(0.9655)
start_p=2.5e-3
min_power=0.4e-6
division=2
power=copy.copy(start_p)
while power>min_power:
    zurich.output1_amp1(power)
    #exp.do_GVg_and_adjust_sitpos()
    #zurich.sigout1_amp1_enabled_param.value(1)
    exp.mech_simple_fun_db()
    #exp.mech_simple_fun_db(start_f=149.7e6,stop_f=150.2e6,step_num_f=50*50)
    #exp.mech_simple_fun_db(start_f=295.7e6,stop_f=296.2e6,step_num_f=50*50)
    #exp.mech_simple_fun_db(start_f=314e6,stop_f=316e6,step_num_f=200*10)
    #zurich.sigout1_amp1_enabled_param.value(0)
    power=power/2
"""

#sweeping cs peak
import time
#Vg=0#init
V_sl=0.25e-3
qdac.ch07.dc_constant_V(V_sl)
while V_sl>-0.25e-3:
    exp.GVG_fun_sensitivity(costum_prefix=f"sl_at_{(V_sl*1e6):4g} uV")
    last_V_sl=qdac.ch07.dc_constant_V()
    V_sl=last_V_sl-25e-6
    print(f"setting ch7 to {V_sl} V")
    time.sleep(5)
    qdac.ch07.dc_constant_V(V_sl)
    
import time
import copy
#Vg=0#init
V_sl_rf=0.5e-3
zurich.output0_amp0(V_sl_rf) 
while V_sl_rf<100e-3:
    %run experiments/zurich_experiments/take_simple_spectrum_0725.py
    last_V_sl_rf=copy.copy(V_sl_rf)
    V_sl_rf=2*last_V_sl_rf
    print(f"setting output0 to {V_sl_rf} V")
    time.sleep(5)
    zurich.output0_amp0(V_sl_rf)

#measure nonlinearity 310725
import time
import copy
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
    %run experiments/zurich_experiments/spectrum_0825_in_constructionb.py
    zurich.set_mixdown(mech_freq)
    #thermal
    %run experiments/zurich_experiments/spectrum_0825_in_construction.py
    mech_freq=exp.exp.max_thermomech_freq
    print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
    zurich.sigout1_amp1_enabled_param.value(1)
    #driven
    %run experiments/zurich_experiments/spectrum_0825_in_constructionc.py
    #driven and swept
    exp.mech_simple_fun_db()
    
    drive_amp=drive_amp*2
    print(f"setting output1 to {drive_amp*1e6:4g} uV")
    time.sleep(5)
