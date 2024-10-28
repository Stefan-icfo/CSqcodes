# sweeps frequency at source of cs and at gate of 5g section simultaneously with difference of RLC frequency
# demodulates at RLC frequency
# Stefan Forstner 

import numpy as np

from instruments import station, zurich,qdac,Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name



import time
from tqdm import tqdm
import matplotlib.pyplot as plt
from utils.CS_utils import centered_moving_average, zurich_phase_voltage_current_conductance


#------User input----------------
#ramp_speed = 1.2e-3 # V/s
tc = 100e-3   # in seconds
vsd_dB = 46#39+20 # attenuation at the source in dB
vsdac =16e-6 # source AC voltage in volt
device_name = 'CD11_D7_c1'
prefix_name = 'Delft_mech_sideofpeak_g2drive'
postfix = '700mK'

source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value
postfix = f"_{round(gate_amplitude_param()*1000,3)}mV on gate@inst,_{round(source_amplitude_param()*1000,3)}mV on source@inst, g1={round(qdac.ch01.dc_constant_V(),2)},g2={round(qdac.ch02.dc_constant_V(),5)},g3={round(qdac.ch03.dc_constant_V(),2)},g4={round(qdac.ch04.dc_constant_V(),5)},g5={round(qdac.ch05.dc_constant_V(),2)},gcs={round(qdac.ch06.dc_constant_V(),5)}"

# exp_name = 'Test 50 K'

mix_down_f = 1.25e6 # RLC frequency
#####################
start_f = 110e6 #Hz unit
stop_f =  130e6 #Hz unit
step_num_f =20*100#1000Hz
#####################





freq_mech = zurich.oscs.oscs1.freq

freq_rlc = zurich.oscs.oscs0.freq
exp_dict = dict(vsdac = vsdac)
exp_name = sample_name(prefix_name,exp_dict,postfix)
freq_rlc(mix_down_f)
freq_mech(start_f)

time.sleep(10)
#----------- defined values------


# ------------------Create a new Experiment-------------------------
freq_sweep = freq_mech.sweep(start=start_f, stop=stop_f, num = step_num_f)
measured_parameter = zurich.demods.demods0.sample  



# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(freq_sweep.parameter)  # 
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('I_rf', 'current', unit='I', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('I_rf_avg', 'current', unit='I', basis=[], setpoints=[freq_sweep.parameter])



# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V',qdac.ch01.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch02_dc_constant_V',qdac.ch02.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V',qdac.ch03.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch04_dc_constant_V',qdac.ch04.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch05_dc_constant_V',qdac.ch05.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V',qdac.ch06.dc_constant_V())
    # for i in range(2):
    I_list=[]
    for f_value in tqdm(freq_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
        freq_mech(f_value)
        time.sleep(1.1*tc) # Wait 1.1 times the time contanst of the lock-in
        measured_value=measured_parameter()
        theta_calc, v_r_calc, I, G = zurich_phase_voltage_current_conductance(measured_value, vsdac)
                
        #G calculation
        I_list.append(I)
    
        datasaver.add_result(('I_rf', I),
                            ('V_r', v_r_calc),
                            ('Phase', theta_calc),
                            (freq_sweep.parameter,f_value))
        
 