# The program does frequency sweep at a fixed gate voltage
# Stefan Forstner with template by Parmeshwar Prasad

import numpy as np
import os

from instruments import qdac, station, keithley2400, Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
from drivers.Qdac_utils import ramp_QDAC_channel

import time
from tqdm import tqdm

k2400=keithley2400
#------User input----------------
slew_rate=1e-2

ramp_source = 10e-6 # V/ms
tc = 20e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd = 5e-3 # source DC voltage in volt
step_v = 2e-3 # source steps
offset = -10e-6 #voltage offset of k2400
offset_i=-44e-12

device_name = 'CD13_E3_C2'
prefix_name = 'g4' 
temp=Triton.T5()
postfix = f'T={temp:4g}'
#Temp=Triton.T5()
#postfix = f"{Temp}K"
vsdkT=temp/11604
if vsdkT<vsd:
    vsd=vsdkT#automatically sets vsd to kT. comment out if wanna do manually
#print(f"vsdkT={vsd}V. ABORT NOW IF FUNKY. u got 10 seconds")
#####################
start_vg = 2 #
stop_vg =  -2
step_num = 400#
#####################

#--------Definition-------------
source = k2400 # source 
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)

gate = qdac.ch04
gate.label = 'VgDC' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)


# ------------------Create a new Experiment-------------------------

vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here

#------------init--------------------
# applied  voltages at the intrument level before attenuation

k2.ramp_k2400(ramp_param=source,final_vg=vsd+offset, step_size = step_v, ramp_speed=ramp_source)


gate.dc_slew_rate_V_per_s(slew_rate)
gate.dc_constant_V(start_vg)
print(f"going to sleep for the time it takes to ramp the gate({abs(start_vg-gate.dc_constant_V())/slew_rate + 30}) plus 30 seconds")
#time.sleep(20)
time.sleep(abs(start_vg-gate.dc_constant_V())/slew_rate + 30)
print("wake up, gates are")

print(gate.dc_constant_V())
#time.sleep(100)


measured_parameter = k2400.curr
# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(vgdc_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas.register_custom_parameter('Conductance', 'G', unit='S', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Resistance', 'R', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[vgdc_sweep.parameter])


# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
        gate.dc_constant_V(vgdc_value)
        time.sleep(1.1*tc+abs(stop_vg-start_vg)/step_num/slew_rate) # Wait 3 times the time contanst of the lock-in plus the ramping time
        measured_value = measured_parameter()-offset_i
        R = vsd/measured_value
        G=1/R
        datasaver.add_result(('Conductance',G),('Resistance', R),('Current', measured_value),(vgdc_sweep.parameter, vgdc_value))
        #print(gate())
        # vgdc_sweep.reverse()

# Ramp down everything
gate.dc_constant_V(0)
k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=ramp_source)
