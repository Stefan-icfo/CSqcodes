#SF 210723
#import instruments
from instruments import k2400
import drivers.k2400 as k2


# essential user inputs
ramp_rate=0.001# steps
max_V=2 #max voltage

#optional user inputs
maxV_neg=-max_V #lowest V
maxV_pos=max_V #highest V
rampspeed=1


#name
device_name ='test_device_'
prefix = 'leaktest_'
postfix = '20mK'
exp_name = sample_name(prefix,postfix)

#define experiment
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)

meas.register_parameter(k2400.volt)
meas.register_custom_parameter('Current', 'I', unit='Amp', basis=[], setpoints=[k2400.volt])
meas.register_custom_parameter('Conductance', 'G', unit='Siemens', basis=[], setpoints=[k2400.volt])


with meas.run() as datasaver:
    
    #ramp up
    for i in range(round(maxV_pos/ramp_rate)+1):
        current_V=i*ramp_rate
        k2.ramp_k2400(ramp_param=k2400.volt,final_vg=current_V, step_size = ramp_rate, ramp_speed=rampspeed)
        datasaver.add_result(('Current', k2400.curr()),('Conductance', k2400.curr()/k2400.volt()),(k2400.volt(), current_V))

    #ramp down
    for i in range(round((maxV_pos-maxV_neg)/ramp_rate)+1):
        current_V=maxV_pos-i*ramp_rate
        k2.ramp_k2400(ramp_param=k2400.volt,final_vg=current_V, step_size = ramp_rate, ramp_speed=rampspeed)
        datasaver.add_result(('Current', k2400.curr()),('Conductance', k2400.curr()/k2400.volt()),(k2400.volt(), current_V))

    #ramp to zero
    for i in range(round(-maxV_neg/ramp_rate)+1):
        current_V=maxV_neg+i*ramp_rate
        k2.ramp_k2400(ramp_param=k2400.volt,final_vg=current_V, step_size = ramp_rate, ramp_speed=rampspeed)
        datasaver.add_result(('Current', k2400.curr()),('Conductance', k2400.curr()/k2400.volt()),(k2400.volt(), current_V))
    



