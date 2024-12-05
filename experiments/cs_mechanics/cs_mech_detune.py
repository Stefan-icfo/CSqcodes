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


#------User input----------------
#costum name
device_name = 'CD11_D7_c1'
prefix_name = 'cs_mech_'
exp_name = "delta_detune_"
postfix = '30mK'

#adjustable hardware params
manual_attenuation_gate=20
tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
att_source_dB = 39 # attenuation at the source in dB# 
att_gate_dB =46+manual_attenuation_gate
mix_down_f = 1.25e6 # RLC frequency


#define delta sweep
idt_point1_x=-1.6747
idt_point1_y=-1.645
idt_point2_x=-1.67108
idt_point2_y=-1.6407
delta=2500e-6

step_vgo_num =6+1 #
xi=0#move along ict (take traces not through centerbut closer to  triple pt)
epsilon_0 =-400e-6#move prependicular to ict (compensate for drift)
start_vgo2,start_vgo1,stop_vgo2,stop_vgo1=make_detuning_axis_noncenterM(idt_point1_x,idt_point1_y,idt_point2_x,idt_point2_y,delta,xi,epsilon_0) 

step_vgo1=np.absolute((start_vgo1-stop_vgo1)/step_vgo_num)
step_vgo2=np.absolute((start_vgo2-stop_vgo2)/step_vgo_num)

vars_to_save=[tc,att_source_dB,att_gate_dB,mix_down_f,idt_point1_x,idt_point1_y,idt_point2_x,idt_point2_y,delta,step_vgo_num]



#inner gate sweep params
start_vgi = -2.2325#-0.788
stop_vgi = -2.2305#-0.776
step_num = 2*50+1#40uV


#frequency sweep params
start_f = 401.8e6 #Hz unit
stop_f =  402.8e6 #Hz unit
step_num_f = 3000+1 #

#source_amp
#source_amplitude_instrumentlevel_GVg = 20e-3 NOT IN USE NOW
source_amplitude_instrumentlevel = 20e-3
gate_amplitude_instrumentlevel = 9.5e-3

#other function params

fit_type='data'
data_avg_num=3
sitfraction="l_max_slope"
freq_sweep_avg_nr=5

return_GVgs=False
return_all_fit_data=False

vars_to_save = [tc, att_source_dB, att_gate_dB, mix_down_f, manual_attenuation_gate, stop_f, start_f, step_num_f,  source_amplitude_instrumentlevel]






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
outer_gate1=qdac.ch02
outer_gate2=qdac.ch04

# ------------------define sweep axes-------------------------

outer_gate1_sweep=outer_gate1.dc_constant_V.sweep(start=start_vgo1, stop=stop_vgo1, num = step_vgo_num)
outer_gate2_sweep=outer_gate2.dc_constant_V.sweep(start=start_vgo2, stop=stop_vgo2, num = step_vgo_num)
outer_gate1_list=list(outer_gate1_sweep)
outer_gate2_list=list(outer_gate2_sweep)

g1_array=np.array(outer_gate1_list)
g2_array=np.array(outer_gate2_list)
delta_array=np.sqrt(abs((g1_array-start_vgo1)**2+(g2_array-start_vgo2)**2))
delta_array-=delta

#amplitudes
#source_amplitude_CNT_GVg=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_GVg)-att_source_dB)
#print(f"source amp at CNT for GVg:{source_amplitude_CNT_GVg*1e6} uV") #NOT IN USE FOR NOW

source_amplitude_CNT=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel)-att_source_dB)
print(f"source amp at CNT:{source_amplitude_CNT*1e6} uV")

gate_amplitude_CNT=d2v(v2d(np.sqrt(1/2)*gate_amplitude_instrumentlevel)-att_gate_dB)
print(f"gate amp at CNT for mech:{gate_amplitude_CNT*1e6} uV")

vars_to_save.extend([gate_amplitude_CNT,source_amplitude_CNT])


#INIT
source_amplitude_param(source_amplitude_instrumentlevel)
gate_amplitude_param(gate_amplitude_instrumentlevel)
outer_gate1.ramp_ch(start_vgo1)
outer_gate2.ramp_ch(start_vgo2)
print("wake up, outer gates are")
print(outer_gate1.dc_constant_V())
#print(outer_auxgate1())
print(outer_gate2.dc_constant_V())




#drive_mag_param = Parameter('drive_mag', label='drive_mag', unit='Vrms',
#                       get_cmd=lambda: drive_now)

freq_param = Parameter('freq', label='freq', unit='Hz',
                       get_cmd=lambda: freq_now)

gateV_param = Parameter('gateV', label='gateV', unit='V',
                       get_cmd=lambda: gateV_now)

#delta_current=0#just to define param
delta_param = Parameter('delta', label='delta', unit='V',
                       get_cmd=lambda: delta_now)



# ----------------Create a measurement-------------------------
experiment_freqdata = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment_freqdata)
meas.register_parameter(delta_param)
meas.register_parameter(freq_param)

#meas.register_parameter(inner_gate_sweep.parameter)   # 
meas.register_custom_parameter('V_rf', 'Amplitude', unit='V', basis=[], setpoints=[delta_param,freq_param])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[delta_param,freq_param])
meas.register_custom_parameter('I_rf', 'current', unit='I', basis=[], setpoints=[delta_param,freq_param])
meas.register_custom_parameter('I_rf_avg', 'current_avg', unit='I', basis=[], setpoints=[delta_param,freq_param])
meas.register_custom_parameter('I_rf/slope', 'current_normalized', unit='a.u.', basis=[], setpoints=[delta_param,freq_param])
#meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])

experiment_G_data = new_experiment(name=exp_name+"_G_data", sample_name=device_name)
meas_aux = Measurement(exp=experiment_G_data)
meas_aux.register_parameter(delta_param)  # 
meas_aux.register_parameter(gateV_param)
meas_aux.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[delta_param,gateV_param])
meas_aux.register_custom_parameter('G_with_sitpos', 'G_with_sitpos', unit='S', basis=[], setpoints=[delta_param,gateV_param])
#meas_aux.register_custom_parameter('V_aux', 'Amplitude_aux', unit='V', basis=[], setpoints=[drive_mag_param,freq_param])
#meas_aux.register_custom_parameter('Phase_aux', 'Phase_aux', unit='rad', basis=[], setpoints=[drive_mag_param,freq_param])

experiment_1D_data = new_experiment(name=exp_name+"_1D_data", sample_name=device_name)
meas_aux_aux = Measurement(exp=experiment_1D_data)
meas_aux_aux.register_parameter(delta_param)  # 
meas_aux_aux.register_custom_parameter('slope', 'slope', unit='S/V', basis=[], setpoints=[delta_param])
meas_aux_aux.register_custom_parameter('sitpos', 'sitpos', unit='V', basis=[], setpoints=[delta_param])
meas_aux_aux.register_custom_parameter('G_at_sitpos', 'G_at_sitpos', unit='S', basis=[], setpoints=[delta_param])
meas_aux_aux.register_custom_parameter('peakpos(max)', 'peakpos(max)', unit='V', basis=[], setpoints=[delta_param])

# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:
    #saving metadata parameters
    qdac.add_dc_voltages_to_metadata(datasaver)
    zurich.save_config_to_metadata(datasaver)
    
    #saving metadata variables
    varnames=[]
    for i in range(len(vars_to_save)):
        varnames.append(get_var_name(vars_to_save[i]))
    save_metadata_var(datasaver.dataset,varnames,vars_to_save)

    with meas_aux.run() as datasaver_aux:
        with meas_aux_aux.run() as datasaver_aux_aux:
            
            for outer_gate1_value, outer_gate2_value,delta_value in tqdm(zip(outer_gate1_sweep, outer_gate2_sweep,delta_array), 
                                                    total=len(outer_gate1_sweep),
                                                    leave=False, 
                                                    desc='Outer Gate Sweep', 
                                                    colour='green'):
                qdac.ramp_multi_ch_fast([outer_gate1, outer_gate2], [outer_gate1_value, outer_gate2_value])
                
                
                single_sweep_results=cs_mechanics_simple_setpoint(start_f=start_f, stop_f=stop_f, step_num_f=step_num_f, 
                                                                start_vg=start_vgi, stop_vg=stop_vgi, step_num=step_num, 
                                                                fit_type=fit_type, data_avg_num=data_avg_num, sitfraction=sitfraction,
                                                                    freq_sweep_avg_nr=freq_sweep_avg_nr, check_at_end=False, 
                                                                    return_GVgs=True, return_all_fit_data=return_all_fit_data)
                Vg=single_sweep_results["Vg_before"]
                G_vals=single_sweep_results["G_vals_before"]
                sitpos=single_sweep_results["sitpos_before"]
                slope=single_sweep_results["slope_before"]

                approx_sitpos_index = np.argmin(np.abs(Vg - sitpos))

                if approx_sitpos_index in {0, len(Vg)-1}:
                    raise ValueError("sitpos is at beginning or end of sweep")
                # Define the approx_sitpos_array
                approx_sitpos_array = copy.copy(G_vals)
                approx_sitpos_array[approx_sitpos_index] = 2*G_vals[approx_sitpos_index]

                datasaver.add_result(('I_rf', single_sweep_results["I"]),
                                    ('I_rf/slope', single_sweep_results["I"]/slope),
                                    ('I_rf_avg', single_sweep_results["I_avg"]),
                                    ('V_rf', single_sweep_results["V"]),
                                    ('Phase', single_sweep_results["Phase"]),
                                    (delta_param,delta_value),
                                    (freq_param,single_sweep_results["freq"]))
                
                datasaver_aux.add_result(('G', G_vals),
                                        ('G_with_sitpos', approx_sitpos_array),
                                        (delta_param,delta_value),
                                    (gateV_param,single_sweep_results["Vg_before"]))
                
                datasaver_aux_aux.add_result(('slope', slope),
                                        (delta_param,delta_value))
                                   
                
            