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
exp.mech_simple_fun_db(start_f=445e6,stop_f=545e6)#seems hemt batteries were screwed or so
zurich.sigout1_amp1_enabled_param.value(0)
exp.do_GVg_and_adjust_sitpos()
zurich.sigout1_amp1_enabled_param.value(1)
exp.mech_simple_fun_db(start_f=120e6,stop_f=220e6)#possibly demod 1kHz off - repeat to be sure
zurich.sigout1_amp1_enabled_param.value(0)
exp.do_GVg_and_adjust_sitpos()
exp.mech_simple_fun_db(start_f=20e6,stop_f=120e6)#possibly demod 1kHz off - repeat to be sure
zurich.sigout1_amp1_enabled_param.value(0)
exp.do_GVg_and_adjust_sitpos()
exp.mech_simple_fun_db(start_f=600e6,stop_f=700e6)#not measured yet
zurich.sigout1_amp1_enabled_param.value(0)


"""

#power sweep

start_p=40e-3
min_power=0.4e-6
division=2
power=copy.copy(start_p)
while power>min_power:
    zurich.output1_amp1(power)
    exp.do_GVg_and_adjust_sitpos()
    zurich.sigout1_amp1_enabled_param.value(1)
    exp.mech_simple_fun_db(start_f=133e6,stop_f=136e6,step_num_f=3000*10)
    exp.mech_simple_fun_db(start_f=149.7e6,stop_f=150.2e6,step_num_f=500*50)
    exp.mech_simple_fun_db(start_f=295.7e6,stop_f=296.2e6,step_num_f=500*50)
    exp.mech_simple_fun_db(start_f=314e6,stop_f=316e6,step_num_f=2000*10)
    zurich.sigout1_amp1_enabled_param.value(0)
    power=power/2



