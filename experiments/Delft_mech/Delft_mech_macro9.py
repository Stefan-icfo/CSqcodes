# The program does frequency sweep at a fixed gate voltage
# Stefan Forstner with template by Parmeshwar Prasad

import numpy as np


from instruments import station, keithley2400, zurich,qdac
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
import time
from tqdm import tqdm

k2400=keithley2400
#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
#ramp_gate = 100e-6 # V/ms
ramp_source = 10e-6 # V/ms
tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd = -8.9e-3 # source DC voltage in volt
step_v = 10e-6 # source steps
offset = -47e-9 #voltage offset of k2400
offset_i=-50e-12
freq = zurich.oscs.oscs1.freq
# device_name = 'tbg_r_3_1'

# prefix_name = 'Conductance_'
# postfix = 'T_10_K_BG'
# exp_name = 'Test 50 K'

#device_name = '100ktest2'
device_name = 'CD11_D7_C1_QDover5_g2_400mV'#andsource40mV
prefix_name = 'Delft_DC' 
postfix = f'5g={qdac.ch01.dc_constant_V()},T=700mK'

#####################
#start_vg = -2 #
#stop_vg =  2#
#step_num = 201      #
#####################
start_f = 20 #MHz unit
stop_f =  100 #MHz unit
step_num =80*100
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
meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('Resistance', 'R', unit='Ohm', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('Conductance', 'G', unit='S', basis=[], setpoints=[freq_sweep.parameter])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    for f_value in tqdm(freq_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
        freq_sweep.set(1e6*f_value)
        time.sleep(3*tc) # Wait 3 times the time contanst of the k2400 
        measured_value = measured_parameter()-offset_i
        R = vsd/measured_value
        G=1/R
        datasaver.add_result(('Conductance',G),('Resistance', R),('Current', measured_value),(freq_sweep.parameter, f_value))
        #print(gate())
        # vgdc_sweep.reverse()

# Ramp down everything
#gate(0)
