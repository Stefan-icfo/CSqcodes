# The program does frequency sweep at a fixed gate voltage
# Stefan Forstner with template by Parmeshwar Prasad

import numpy as np


from instruments import station, zurich
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name

import time
from tqdm import tqdm
gate_ramp_slope = 1e-2# V/s
tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd_dB = 39 # attenuation at the source in dB
vsdac = 0.1e-3 # source AC voltage in volt, for now set manually
device_name = 'CD11_D7_C1'
prefix_name = 'Conductance_qdac_zurich_'
postfix = '700mK'
# exp_name = 'Test 50 K'

mix_down_f = 1.25e6 # RLC frequency
###################
#----------- defined values------
#####################
gain_RT = 200       #
gain_HEMT = 5.64   #%
Z_tot = 7521        #
###################

#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
#ramp_gate = 100e-6 # V/ms
ramp_source = 10e-6 # V/ms
tc = 10e-3   # in seconds. Doesn't get overwritten by ZI called value.

mix_down_f = 1.25e6 

freq_rlc = zurich.oscs.oscs0.freq
freq_rlc(mix_down_f)

device_name = 'CD11_D7_C1_charge sensor'
prefix_name = 'mechanics_1' 
postfix = '700mK'

start_f = 10 #MHz unit
stop_f =  600 #MHz unit
step_num =560*10
#--------Definition-------------



freq = zurich.oscs.oscs0.freq
freq_sweep = freq.sweep(start=start_f, stop=stop_f, num = step_num)  # gate parameter and the instrument will go here  # measured parameters will go here
measured_parameter = zurich.demods.demods1.sample 




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
starttime=time.time()
with meas.run() as datasaver:
    # for i in range(2):
    for f_value in tqdm(freq_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
        current_time=time.time()-starttime
        freq_sweep.set(1e6*f_value)
        time.sleep(3*tc) # Wait 3 times the time contanst of the k2400 
        measured_value = measured_parameter()
        R = vsdac/measured_value
        G=1/R
        datasaver.add_result(('Conductance',G),('Resistance', R),('Current', measured_value),(freq_sweep.parameter, f_value))
