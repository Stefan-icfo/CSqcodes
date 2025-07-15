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

       
"""
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



"""

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

"""
qdac.ch06.dc_constant_V(0.9655)
start_p=20e-3
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


