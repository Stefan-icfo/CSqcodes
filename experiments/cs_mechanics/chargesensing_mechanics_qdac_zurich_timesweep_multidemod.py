# sweeps frequency at source of cs and at gate of 5g section simultaneously with difference of RLC frequency
# demodulates at RLC frequency
# Stefan Forstner 

import numpy as np

from instruments import station, zurich,qdac,Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name

from qcodes import Parameter

import time
from tqdm import tqdm
import matplotlib.pyplot as plt
from utils.CS_utils import centered_moving_average, zurich_phase_voltage_current_conductance


#------User input----------------
#ramp_speed = 1.2e-3 # V/s
tc = 1e-3   # in seconds; minimal to get many datapoints!
vsd_dB = 45 # attenuation at the source in dB
vsdac =16e-6 # source AC voltage in volt
device_name = 'CD11_D7_c1'
prefix_name = 'chargesensing_mechanics_'
postfix = '22mK'
Temp=Triton.MC()
source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value
mech_f=121.9099e6 

postfix = f"_mech_f={mech_f},{round(gate_amplitude_param()*1000,3)}mV on gate@inst,_{round(source_amplitude_param()*1000,3)}mV on source@inst, g1={round(qdac.ch01.dc_constant_V(),2)},g2={round(qdac.ch02.dc_constant_V(),5)},g3={round(qdac.ch03.dc_constant_V(),2)},g4={round(qdac.ch04.dc_constant_V(),5)},g5={round(qdac.ch05.dc_constant_V(),2)},gcs={round(qdac.ch06.dc_constant_V(),5)}"

# exp_name = 'Test 50 K'

mix_down_f = 1.25e6 # RLC frequency
#####################
#start_f = 200e6 #Hz unit
#stop_f =  400e6 #Hz unit
#step_num_f =200*1000#1kHz

#####################
timeout=1000#s





#qdac.ch01.dc_constant_V(+0.6)
#qdac.ch02.dc_constant_V(-1.96570)#y
#qdac.ch03.dc_constant_V(+1.1)
#qdac.ch04.dc_constant_V(-1.962365)#x
#qdac.ch05.dc_constant_V(0.6)
#qdac.ch06.dc_constant_V(-0.33963)#side of cb peak




freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq
exp_dict = dict(vsdac = vsdac)
exp_name = sample_name(prefix_name,exp_dict,postfix)
freq_rlc(mix_down_f)
freq_mech(mech_f)
freq_rf(mech_f-mix_down_f)
time.sleep(10)
#----------- defined values------
#####################
gain_RT = 200       #
gain_HEMT = 5.64   #
Z_tot = 7521        #
###################


# ------------------Create a new Experiment-------------------------
#freq_sweep = freq_rf.sweep(start=start_f, stop=stop_f, num = step_num_f)
measured_parameter = zurich.demods.demods2.sample  

start_time = time.time()
time_param = Parameter('time', label='Time', unit='s',
                       get_cmd=lambda: time.time() - start_time)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(time_param)  # 
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[time_param])
meas.register_custom_parameter('I_rf', 'current', unit='I', basis=[], setpoints=[time_param])
#meas.register_custom_parameter('I_rf_avg', 'current', unit='I', basis=[], setpoints=[time_param])



# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V',qdac.ch01.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch02_dc_constant_V',qdac.ch02.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V',qdac.ch03.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch04_dc_constant_V',qdac.ch04.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch05_dc_constant_V',qdac.ch05.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V',qdac.ch06.dc_constant_V())

    datasaver.dataset.add_metadata('f_mech',mech_f)
    datasaver.dataset.add_metadata('source_amp_instr',source_amplitude_param())
    datasaver.dataset.add_metadata('gate_amp_instr',gate_amplitude_param())
    # for i in range(2):
    #I_list=[]
    while time_param()<timeout:
        #freq_rf(f_value-freq_rlc())
        #freq_mech(f_value)
        time.sleep(1.1*tc) # Wait 1.1 times the time contanst of the lock-in
        measured_value=measured_parameter()
        theta_calc, v_r_calc, I, G = zurich_phase_voltage_current_conductance(measured_value, vsdac)
                
        #G calculation
        #I_list.append(I)
    
        datasaver.add_result(('I_rf', I),
                            ('V_r', v_r_calc),
                            ('Phase', theta_calc),
                            (time_param,time_param()))
        
    
    #datasaver.dataset.add_metadata('rohde.power()',rohde.power())
# Ramp down everything
#print(gate())
#gate(0)

#AvgI=centered_moving_average(I_list,120)
#plt.plot(list(freq_sweep),AvgI)
#plt.title(f'meas{run_id}')
#plt.show()    

#for i in range(8):
#        param1 = getattr(zurich.sigouts.sigouts0.enables, f'enables{i}')
#        param2 = getattr(zurich.sigouts.sigouts1.enables, f'enables{i}')
#        param1.value(0)
#        param2.value(0)


#rohde.power(-50)
