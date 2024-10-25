# The program does frequency sweep at a fixed gate voltage
# Stefan Forstner with template by Parmeshwar Prasad

import numpy as np
import os

from instruments import  station, keithley2400, zurich, rohde, qdac
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
import time
from tqdm import tqdm
#from utils.CS_utils import save_metadata_var, get_var_name

k2400=keithley2400
#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
#ramp_gate = 100e-6 # V/ms
ramp_source = 10e-6 # V/ms
tc = 10e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd = 350e-6 # source DC voltage in volt
step_v = 10e-6 # source steps
offset = -10e-6 #voltage offset of k2400
offset_i=-50e-12
#freq = zurich.oscs.oscs0.freq
freq = rohde.frequency
# device_name = 'tbg_r_3_1'
rohde.power(-10)

# prefix_name = 'Conductance_'
# postfix = 'T_10_K_BG'
# exp_name = 'Test 50 K'

#device_name = '100ktest2'
device_name = 'CD11_D7_C1_'
prefix_name = 'mechanics_-10db_200microVsource10Khzstepzoomg1-4g2-2.34g3-1g4-4g5-4zoom' 
postfix = '18mK'

#####################
#start_vg = -2 #
#stop_vg =  2#
#step_num = 201      #
#####################
start_f = 100e6 #MHz unit
stop_f =  140e6 #MHz unit
step_num = 40*1000
#--------Definition-------------
source = k2400 # source 
  # channel of the gate
freq.label = 'Frequency (MHz)' # Change the label of the gate chaneel
instr_dict = dict(freq=[freq])
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)


# ------------------Create a new Experiment-------------------------

freq_sweep = freq.sweep(start=start_f, stop=stop_f, num = step_num)  # gate parameter and the instrument will go here
measured_parameter = k2400.curr  # measured parameters will go here

#------------init--------------------
# applied  voltages at the intrument level before attenuation

k2.ramp_k2400(ramp_param=source,final_vg=vsd+offset, step_size = step_v, ramp_speed=ramp_source)



# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(freq_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas.register_custom_parameter('Resistance', 'R', unit='Ohm', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('Conductance', 'G', unit='S', basis=[], setpoints=[freq_sweep.parameter])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V',qdac.ch01.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch02_dc_constant_V',qdac.ch02.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V',qdac.ch03.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch04_dc_constant_V',qdac.ch04.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch05_dc_constant_V',qdac.ch05.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V',qdac.ch06.dc_constant_V())
    # for i in range(2):
    for f_value in tqdm(freq_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
        freq_sweep.set(f_value)
        time.sleep(3*tc) # Wait 3 times the time contanst of the k2400 
        measured_value = measured_parameter()-offset_i
        R = vsd/measured_value
        G=1/R
        datasaver.add_result(('Conductance',G),('Resistance', R),('Current', measured_value),(freq_sweep.parameter, f_value))
        #print(gate())
        # vgdc_sweep.reverse()

# Ramp down everything
#gate(0)
k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=ramp_source)
