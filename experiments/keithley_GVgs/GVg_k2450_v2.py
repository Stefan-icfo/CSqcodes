# GVG with K2450
# Stefan Forstner with template by Parmeshwar Prasad
# same as equivalent k2400 code, not quite working yet as of 28/04/23

import numpy as np
import os

from instruments import bilt, manual, station, k2450
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
#import drivers.k2450 as k2
import time
from tqdm import tqdm


#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
ramp_gate = 100e-6 # V/ms
ramp_source = 100e-6 # V/ms
tc = 10e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd = 10e-3 # source DC voltage in volt
step_v = 1e-3 # source steps
offset = 75e-6 #voltage offset of k2400
# device_name = 'tbg_r_3_1'
# prefix_name = 'Conductance_'
# postfix = 'T_10_K_BG'
# exp_name = 'Test 50 K'

device_name = 'test'
prefix_name = 'Conductance_1' 
postfix = '300K'

#####################
start_vg = -3    #
stop_vg =  3   #
step_num = 1201      #
#####################



#--------Definition-------------
#source = k2450  # source 
gate = bilt.ch01.v  # channel of the back gate
gate.label = 'VgDC' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)


# ------------------Create a new Experiment-------------------------

vgdc_sweep = gate.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here


#------------init--------------------
# applied  voltages at the intrument level before attenuation
#source.mode('VOLT')
#source.compliancev = 21
#source.compliancei = 1e-6

#k2.ramp_k2450(ramp_param=source,final_vg=vsd+offset, step_size = step_v, ramp_speed=ramp_source)
k2450.source.function("voltage")
k2450.source.sweep_setup(0, step_v, vsd)

#k2450.source.voltage(vsd)

k2450.sense.function("current")
measured_parameter = k2450.sense.current  # measured parameters will go here

bilt.channels.output_mode('ramp')
bilt.channels.ramp_slope(ramp_gate)
gate(start_vg)
time.sleep(5)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(vgdc_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas.register_custom_parameter('Resistance', 'R', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Conductance', 'G', unit='S', basis=[], setpoints=[vgdc_sweep.parameter])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
        vgdc_sweep.set(vgdc_value)
        time.sleep(3*tc) # Wait 2 times the time contanst of the lock-in
        measured_value = measured_parameter()
        R = vsd/measured_value
        G=1/R
        datasaver.add_result(('Conductance',G),('Resistance', R),('Current', measured_value),(vgdc_sweep.parameter, vgdc_value))
    
        # vgdc_sweep.reverse()

# Ramp down everything
gate(0)
k2450.source.sweep_setup(vsd, step_v, 0)
