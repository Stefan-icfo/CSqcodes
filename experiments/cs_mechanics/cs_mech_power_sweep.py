import numpy as np
import scipy as scp
import os
import copy
import math
import time
from tqdm import tqdm
import matplotlib.pyplot as plt

from instruments import   station, qdac,  Triton, zurich
from qcodes.dataset import Measurement, new_experiment
from qcodes import Parameter

from utils.sample_name import sample_name
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk

from utils.CS_utils import *
from experiments.cs_mechanics.cs_mechanics_simple_setpoint_adjust_fun import *
import experiment_parameters


#------User input----------------
#costum name
device_name = experiment_parameters.device_name#'CD11_D7_c1'
prefix_name = 'cs_mech_powersweep_'

postfix = '34mK'
postfix = f"_g1={round(qdac.ch01.dc_constant_V(),4)},g2={round(qdac.ch02.dc_constant_V(),4)},g3={round(qdac.ch03.dc_constant_V(),4)},g4={round(qdac.ch04.dc_constant_V(),4)},g5={round(qdac.ch05.dc_constant_V(),4)}"
exp_name = prefix_name+device_name+postfix
#adjustable hardware params

tc = 30e-3#experiment_parameters.tc   # in seconds. Doesn't get overwritten by ZI called value.
att_source_dB = experiment_parameters.attn_dB_source # attenuation at the source in dB# 
att_gate_dB =experiment_parameters.attn_dB_gate
mix_down_f = experiment_parameters.mix_down_f # RLC frequency
#source_amplitude_instrumentlevel_GVg = 20e-3

#power_sweep
start_value=9e-3
length=12
instr_power_sweep=[start_value / (1.5 ** i) for i in range(length)]
#instr_power_sweep=10*[1e-6]

#gate sweep params
start_vg = 0.962#-0.788
stop_vg = 0.967#-0.776
step_num = 5*20#40uV
step_vgi=np.absolute((start_vg-stop_vg)/step_num)

#frequency sweep params
start_f = 300.1e6 #Hz unit
stop_f =  300.6e6#Hz unit
step_num_f = 5*500+1 #

#source_amp
source_amplitude_instrumentlevel_GVg = experiment_parameters.source_amplitude_instrumentlevel_GVg
source_amplitude_instrumentlevel_mech = 20e-3

#other function params

fit_type='data'
data_avg_num=7
sitfraction=0.85
#"l_max_slope"
freq_sweep_avg_nr=21
#return_GVgs=False
return_all_fit_data=False

vars_to_save = [tc, att_source_dB, att_gate_dB, mix_down_f,  source_amplitude_instrumentlevel_GVg, start_value, length, instr_power_sweep, start_vg, stop_vg, step_num, step_vgi, stop_f, start_f, step_num_f, source_amplitude_instrumentlevel_GVg, source_amplitude_instrumentlevel_mech]


source_amplitude_CNT_GVg=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_GVg)-att_source_dB)
print(f"source amp at CNT for GVg:{source_amplitude_CNT_GVg*1e6} uV")

source_amplitude_CNT_mech=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_mech)-att_source_dB)
print(f"source amp at CNT for mech:{source_amplitude_CNT_mech*1e6} uV")

vars_to_save.extend([source_amplitude_CNT_GVg,source_amplitude_CNT_mech])



#channel assignments
#cs_gate=qdac.ch06  # swept gate voltage
source_amplitude_param = zurich.output0_amp0
gate_amplitude_param = zurich.output1_amp1
gate_rf_enabled_param = zurich.sigout1_amp1_enabled_param
freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq
measured_parameter = zurich.demod2 
measured_aux_parameter = zurich.demod0


drive_mag_param = Parameter('drive_mag', label='drive_mag', unit='Vrms')

freq_param = Parameter('freq', label='freq', unit='Hz')

gateV_param = Parameter('gateV', label='gateV', unit='V')



# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(drive_mag_param)
meas.register_parameter(freq_param)
#meas.register_parameter(inner_gate_sweep.parameter)   # 
meas.register_custom_parameter('I_rf', 'current', unit='I', basis=[], setpoints=[drive_mag_param,freq_param])
meas.register_custom_parameter('I_rf_avg', 'current_avg', unit='I', basis=[], setpoints=[drive_mag_param,freq_param])
meas.register_custom_parameter('V_rf', 'Amplitude', unit='V', basis=[], setpoints=[drive_mag_param,freq_param])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[drive_mag_param,freq_param])

#meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])

experiment_aux = new_experiment(name=exp_name+"aux", sample_name=device_name)
meas_aux = Measurement(exp=experiment_aux)
meas_aux.register_parameter(drive_mag_param)  # 
meas_aux.register_parameter(gateV_param)
meas_aux.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[drive_mag_param,gateV_param])
#meas_aux.register_custom_parameter('V_aux', 'Amplitude_aux', unit='V', basis=[], setpoints=[drive_mag_param,freq_param])
#meas_aux.register_custom_parameter('Phase_aux', 'Phase_aux', unit='rad', basis=[], setpoints=[drive_mag_param,freq_param])

# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:
    #saving metadata parameters
    qdac.add_dc_voltages_to_metadata(datasaver)
    zurich.save_config_to_metadata(datasaver)
    
    #saving metadata variables
    varnames=[]
  #  for i in range(len(vars_to_save)):
  #      varnames.append(get_var_name(vars_to_save[i]))
  #  save_metadata_var(datasaver.dataset,varnames,vars_to_save)

    with meas_aux.run() as datasaver_aux:
        zurich.sigout1_amp1_enabled_param.value(1)
        last_gate_amplitude_CNT=0#init
        i=0
        for instr_magVrms in tqdm(instr_power_sweep, leave=False, desc='outer Gate Sweep', colour = 'green'): 
            i+=1
            gate_amplitude_instrumentlevel=instr_magVrms
            gate_amplitude_param(instr_magVrms)
            gate_amplitude_CNT=d2v(v2d(np.sqrt(1/2)*gate_amplitude_param())-att_gate_dB)
            print(f"__wanna be {instr_magVrms*1e6} uV and is {round(gate_amplitude_param()*1e6,4)} uV")
            if gate_amplitude_CNT==last_gate_amplitude_CNT:
                gate_amplitude_CNT+=i*1e-12#add to make sure values are distinct
            print(f"step_num_f={step_num_f}")
            single_sweep_results=cs_mechanics_simple_setpoint(start_f=start_f, stop_f=stop_f, step_num_f=step_num_f, 
                                                              start_vg=start_vg, stop_vg=stop_vg, step_num=step_num, 
                                                              fit_type=fit_type, data_avg_num=data_avg_num, sitfraction=sitfraction,
                                                                freq_sweep_avg_nr=freq_sweep_avg_nr, tc=tc,check_at_end=False, 
                                                                return_GVgs=True, return_all_fit_data=return_all_fit_data,switch_off_gate_drive_for_GVg=True)
            
            datasaver.add_result(('I_rf', single_sweep_results["I"]),
                                 ('I_rf_avg', single_sweep_results["I_avg"]),
                                ('V_rf', single_sweep_results["V"]),
                                ('Phase', single_sweep_results["Phase"]),
                                (drive_mag_param,gate_amplitude_CNT),
                                (freq_param,single_sweep_results["freq"]))
            
            datasaver_aux.add_result(('G', single_sweep_results["G_vals_before"]),
                                     (drive_mag_param,gate_amplitude_CNT),
                                (gateV_param,single_sweep_results["Vg_before"]))
            
            last_gate_amplitude_CNT=copy.copy(gate_amplitude_CNT)             