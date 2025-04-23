# Burns leaks on a gate
# Stefan Forstner 

import numpy as np
import os

from instruments import  station, keithley2400
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
import time
from tqdm import tqdm


#------User input----------------

ramp_source = 100e-6 # V/ms
step_v = 10e-6 # source steps
#offset = -10e-6 #voltage offset of k2400
#offset_i=-50e-12

device_name='test'
#device_name = 'CD08_G7_G6_g5'
prefix_name = 'annealing' 
postfix = 'RT'

#####################
start_vs = 0  #
stop_vs =  2 #
step_num = 201      #
#####################

#--------Definition-------------
source = keithley2400 # source 

exp_name = prefix_name+postfix


# ------------------Create a new Experiment-------------------------
burn_sweep = source.volt.sweep(start=start_vs, stop=stop_vs, num = step_num)  # gate parameter and the instrument will go here
measured_parameter = keithley2400.curr  # measured parameters will go here

#------------init--------------------
# applied  voltages at the intrument level before attenuation
#source.mode('VOLT')
#source.compliancev = 21
#source.compliancei = 100e-6



time.sleep(2)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(burn_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
#meas.register_custom_parameter('Resistance', 'R', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[burn_sweep.parameter])
#meas.register_custom_parameter('Conductance', 'G', unit='S', basis=[], setpoints=[vgdc_sweep.parameter])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    for burn_value in tqdm(burn_sweep, leave=False, desc='Voltage Ramp', colour = 'green'):
        #print("startloop")
        #k2.ramp_k2400(ramp_param=source,final_vs=burn_value, step_size = step_v, ramp_speed=ramp_source)
        source.volt(burn_value)
        #print("step0")
        burn_sweep.set(burn_value)
        #print(burn_value)
        #time.sleep(abs(stop_vs-start_vs)/step_num/ramp_source/1e3) # Wait 3 times the time contanst of the lock-in plus the ramping time
        measured_value = measured_parameter()
        #print("step1")
        datasaver.add_result(('Current', measured_value),(burn_sweep.parameter, burn_value))
        #print(gate())
        # vgdc_sweep.reverse()

# Ramp down everything
print("loop done")
source.volt(0)
