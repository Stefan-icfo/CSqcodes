# The program does frequency sweep at a fixed gate voltage
# Stefan Forstner with template by Parmeshwar Prasad

import numpy as np
import os

from instruments import bilt, station, keithley2400
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
import time
from tqdm import tqdm

k2400=keithley2400
#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
ramp_gate = 100e-6 # V/ms
ramp_source = 10e-6 # V/ms
tc = 10e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd = 20e-3 # source DC voltage in volt
step_v = 10e-6 # source steps
offset = -10e-6 #voltage offset of k2400
offset_i=-50e-12
# device_name = 'tbg_r_3_1'

# prefix_name = 'Conductance_'
# postfix = 'T_10_K_BG'
# exp_name = 'Test 50 K'

#device_name = '100ktest2'
device_name = 'CD08_G7_F6_charge sensor'
prefix_name = 'Conductance_1' 
postfix = '10K'

#####################
start_vg = -2 #
stop_vg =  2#
step_num = 201      #
#####################

#--------Definition-------------
source = k2400 # source 
gate = bilt.ch01.v  # channel of the gate
gate.label = 'VgDC' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)


# ------------------Create a new Experiment-------------------------

vgdc_sweep = gate.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here
measured_parameter = k2400.curr  # measured parameters will go here

#------------init--------------------
# applied  voltages at the intrument level before attenuation
source.mode('VOLT')
source.compliancev = 21
source.compliancei = 1e-6
k2.ramp_k2400(ramp_param=source,final_vg=vsd+offset, step_size = step_v, ramp_speed=ramp_source)

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
        time.sleep(3*tc+abs(stop_vg-start_vg)/step_num/ramp_gate/1e3) # Wait 3 times the time contanst of the lock-in plus the ramping time
        measured_value = measured_parameter()-offset_i
        R = vsd/measured_value
        G=1/R
        datasaver.add_result(('Conductance',G),('Resistance', R),('Current', measured_value),(vgdc_sweep.parameter, vgdc_value))
        #print(gate())
        # vgdc_sweep.reverse()

# Ramp down everything
gate(0)
k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=ramp_source)
