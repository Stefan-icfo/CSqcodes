# Program for the IV measurement using current applifier
import numpy as np
import os
from instruments import station, k2400#,  #triton2
from qcodes.dataset import Measurement, new_experiment
#from experiments_resonance.finisher import end_game
import drivers.k2400 as k2
import time
from tqdm import tqdm

#------User input----------------
wait_time = 0.05   # rise time in current amplifier in seconds
ramp_mode = 'ramp'  #source ramp mode
ramp_speed = 1 # V/ms
exp_name = 'IV_'
# exp_name = 'Test 1 M Ohm '

device_name = 'test_'


#--------Definition-------------
source = k2400  # source 

vsd_range = -1e-3-1e-3 #
step_v = 5e-6
num_of_step = 201

instr_dict = dict(gate=[source])

# ------------------Create a new Experiment-------------------------
i_read = k2400.curr  # measured parameters will go here

vsd_sweep = source.volt.sweep(start=-vsd_range, stop=vsd_range, num = num_of_step)  # source the instrument will go here

#------------init--------------------
#source.mode('VOLT')
#source.compliancev = 2.1
#source.compliancei = 1e-6
k2.ramp_k2400(source,final_vg=vsd_sweep[0], step_size = step_v, ramp_speed=1)

# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(vsd_sweep.parameter)  # register the 1st independent parameter

meas.register_custom_parameter('i', 'I', unit='A', basis=[], setpoints=[vsd_sweep.parameter])  # registering the custom parameters Resistance which is calculated from other parameters such as source drian voltage(Vsd), measured lock-in voltage(measured_param) and Reference Resi
meas.register_custom_parameter('resistance', 'Resistance', unit='Ohm', basis=[], setpoints=[vsd_sweep.parameter])  # registering the custom parameters Resistance which is calculated from other parameters such as source drian voltage(Vsd), measured lock-in voltage(measured_param) and Reference Resi

meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)
time.sleep(5)

with meas.run() as datasaver:
    for vsd_value in tqdm(vsd_sweep, leave=False, desc='Source sweep',colour = 'green'):
        k2.ramp_k2400(source,final_vg=vsd_value, step_size = step_v, ramp_speed=ramp_speed)
        time.sleep(wait_time)
        i_sample = i_read()
        R_sample = vsd_value/i_sample #calculate the sample resistance
        datasaver.add_result(('resistance', R_sample), ('i', i_sample), (vsd_sweep.parameter, vsd_value))

# k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=1)
