# weeps gates of both sections for charge sensing
# Stefan Forstner



import numpy as np


from instruments import   station, qdac,  Triton, zurich
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
#from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk

import time
from tqdm import tqdm
import scipy as scp
import matplotlib.pyplot as plt
from utils.CS_utils import breit_wigner_fkt, breit_wigner_detuning, zurich_phase_voltage_current_conductance, zurich_phase_voltage_current_conductance_compensate, idt_perpendicular_angle, make_detuning_axis, save_metadata_var, get_var_name, zurich_working
from experiment_functions.CS_functions import *
import os
from qcodes import Parameter
import copy
import math

script_path = __file__
print("Full path of the script:", script_path)
debug=False
#------User input----------------
slew_rate=1e-2
tc = 30e-3   # in seconds.
att_source_dB = 39 # attenuation at the source in dB
att_gate_dB =46+20
#vsdac = 200e-6 # source AC voltage in volt
device_name = 'CD11_D7_C1'

prefix_name = 'cs_mechg2_adjustg2andg4toremoverightdot'#
postfix = 'coolingfrom700mK'


x_avg=+4.38e-6
y_avg=-4.41e-6
mix_down_f = 1.25e6 # RLC frequency

step_vgo_num = 45#
start_vgo=-1.4
stop_vgo=0
#others_offset=+1
crosscap24=-0.1

source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value
postfix = f"_{round(gate_amplitude_param()*1000,3)}mV on gate@inst,_{round(source_amplitude_param()*1000,3)}mV on source@inst, g1={round(qdac.ch01.dc_constant_V(),2)},g3={round(qdac.ch03.dc_constant_V(),2)},g5={round(qdac.ch05.dc_constant_V(),2)},gcs={round(qdac.ch06.dc_constant_V(),2)}"

#outer gates (slow axis)
step_vgo=np.absolute((start_vgo-stop_vgo)/step_vgo_num)
vars_to_save=[slew_rate,tc,att_source_dB,att_gate_dB,x_avg,y_avg,mix_down_f,step_vgo_num]

#inner gate voltage range (fast axis, CS)
#####################
start_vgi = -1.26#-0.788
stop_vgi = -1.225#-0.776
step_vgi_num = 35*10#40uV
step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)

#initial_guess = [-2.45, 1e-4, 3e-7]#initial guess for peakV, Gamma,height for first GVg
sitfraction="l_max_slope"#"l_max_slope"#0.6#where to sit on Coulomb peak. For now on left side

vars_to_save.extend([start_vgi,stop_vgi,step_vgi_num])
#####################
start_f =122e6#Hz unit
stop_f = 126e6 #Hz unit
step_num_f = 4*1000#1kHz

vars_to_save.extend([start_f,stop_f,step_num_f])

source_amplitude_instrumentlevel_GVg = 20e-3
source_amplitude_CNT_GVg=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_GVg)-att_source_dB)
print(f"source amp at CNT for GVg:{source_amplitude_CNT_GVg*1e6} uV")
source_amplitude_instrumentlevel_mech = 20e-3
source_amplitude_CNT_mech=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_mech)-att_source_dB)
print(f"source amp at CNT for mech:{source_amplitude_CNT_mech*1e6} uV")
gate_amplitude_instrumentlevel = 5e-3
gate_amplitude_CNT=d2v(v2d(np.sqrt(1/2)*gate_amplitude_instrumentlevel)-att_gate_dB)
print(f"gate amp at CNT for mech:{gate_amplitude_CNT*1e6} uV")

vars_to_save.extend([source_amplitude_instrumentlevel_GVg, source_amplitude_CNT_GVg, source_amplitude_instrumentlevel_mech,source_amplitude_CNT_mech,gate_amplitude_instrumentlevel,gate_amplitude_CNT])

#--------Definitions-------------

#swept contacts
cs_gate=qdac.ch06.dc_constant_V  # swept gate voltage
qdac.ch06.dc_slew_rate_V_per_s(slew_rate)


outer_gate=qdac.ch04.dc_constant_V
#outer_other_gates=[qdac.ch04.dc_constant_V,qdac.ch05.dc_constant_V]
currentvgcs=cs_gate()
cs_gate((start_vgi+currentvgcs)/2)
time.sleep(10)
cs_gate(start_vgi)
sleeptime=10
print(f"manual sleeptime={sleeptime} s")
time.sleep(sleeptime)
print("wake up, gates are")

print(cs_gate())


cs_gate.label = 'CS_gate(inner)' # Change the label of the source chaneel
exp_dict = dict(uVrfgate = gate_amplitude_CNT*1e6)
exp_name = sample_name(prefix_name,exp_dict,postfix)

# ------------------define DC sweep axes-------------------------

outer_gate_sweep=outer_gate.sweep(start=start_vgo, stop=stop_vgo, num = step_vgo_num)
outer_gate1_list=list(outer_gate_sweep)
g1_array=np.array(outer_gate1_list)
cs_sweep=cs_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)

#-------------------define RF frequencies-----------------------

freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq
freq_rlc(mix_down_f)
freq_mech(start_f)
freq_rf(start_f - mix_down_f)

#-------------------define RF sweep axis-----------------------
freq_sweep = freq_rf.sweep(start=start_f, stop=stop_f, num = step_num_f)
freq_sweep_list=list(freq_sweep)
measured_parameter = zurich.demods.demods2.sample  
measured_aux_parameter = zurich.demods.demods0.sample

source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value
gate_rf_enabled_param = getattr(zurich.sigouts.sigouts1.enables, f'enables{1}')


#delta_current=0#just to define param

# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(outer_gate_sweep.parameter)
meas.register_parameter(freq_sweep.parameter)
#meas.register_parameter(inner_gate_sweep.parameter)   # 
meas.register_custom_parameter('V_rf', 'Amplitude', unit='V', basis=[], setpoints=[outer_gate_sweep.parameter,freq_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[outer_gate_sweep.parameter,freq_sweep.parameter])
meas.register_custom_parameter('I_rf', 'current', unit='I', basis=[], setpoints=[outer_gate_sweep.parameter,freq_sweep.parameter])
#meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])

experiment_aux = new_experiment(name=exp_name+"aux", sample_name=device_name)
meas_aux = Measurement(exp=experiment_aux)
meas_aux.register_parameter(outer_gate_sweep.parameter)  # 
meas_aux.register_parameter(cs_sweep.parameter)
meas_aux.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[outer_gate_sweep.parameter,cs_sweep.parameter])
meas_aux.register_custom_parameter('V_aux', 'Amplitude_aux', unit='V', basis=[], setpoints=[outer_gate_sweep.parameter,cs_sweep.parameter])
meas_aux.register_custom_parameter('Phase_aux', 'Phase_aux', unit='rad', basis=[], setpoints=[outer_gate_sweep.parameter,cs_sweep.parameter])
meas_aux.register_custom_parameter('approx_sitpos', 'approx_sitpos', unit='V', basis=[], setpoints=[outer_gate_sweep.parameter,cs_sweep.parameter])
# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:
    #saving metadata parameters
    datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V',qdac.ch01.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch02_dc_constant_V',qdac.ch02.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V',qdac.ch03.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch04_dc_constant_V',qdac.ch04.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch05_dc_constant_V',qdac.ch05.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V',qdac.ch06.dc_constant_V())
    datasaver.dataset.add_metadata('script_file',script_path)
    
    #saving metadata variables
    varnames=[]
    for i in range(len(vars_to_save)):
        varnames.append(get_var_name(vars_to_save[i]))
    save_metadata_var(datasaver.dataset,varnames,vars_to_save)

    with meas_aux.run() as datasaver_aux:

        cs_sweep_list = list(cs_sweep) #
        
        sitposlist=[]
        current_csvg=copy.copy(start_vgi)
        for outer_gate_value in tqdm(outer_gate_sweep, leave=False, desc='outer Gate Sweep', colour = 'green'): #slow axis loop (gate)
            #print('temperature')
            #Triton.MC()
            auxgatecurrentvalue=qdac.ch02.dc_constant_V()
            time.sleep(0.5)
            outer_gate(outer_gate_value)
            time.sleep(0.5)
            qdac.ch02.dc_constant_V(auxgatecurrentvalue+crosscap24*step_vgo)
            print(f"settingg2 to {auxgatecurrentvalue+crosscap24}")
            
            time.sleep(abs(step_vgo/slew_rate)) # Wait  the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
           
            #lists for frequency sweeps
            I_list=[]
            V_list=[]
            Phase_list=[]

            #run GVg to find sitpos
            freq_rf(mix_down_f)#set to RLC frequency
            #########################
            gate_rf_enabled_param.value(0)#switch off gate output
            source_amplitude_param(source_amplitude_instrumentlevel_GVg)#set source aplitude
            cs_gate((start_vgi+currentvgcs)/2)
            time.sleep(10)
            cs_gate(start_vgi)
            time.sleep(abs(start_vgi-current_csvg)/slew_rate)

            if debug:
                print("starting GVg")
            
            Glist,V_aux_list,I_aux_list,Phase_aux_list=GVG_simple(gate_sweep=cs_sweep,
                                                   measured_parameter=measured_aux_parameter,
                                                   step_sleep_time=1.1*tc+abs(step_vgi)/slew_rate,
                                                   vsdac=source_amplitude_CNT_GVg,
                                                   x_avg=x_avg,
                                                   y_avg=y_avg)
            
            
            fit_calculated_sitpos=fit_and_find_sitpos_singlepeak(gate_sweep=cs_sweep,
                                                  Glist=Glist,
                                                  initial_guess=None, 
                                                  sitfraction=sitfraction,
                                                  return_full_fit_data=False)
            #plt.plot(cs_sweep,Glist)
            #plt.plot(fit_calculated_sitpos,0.75*max(Glist),"o")
            #plt.show()
            #initial_guess=[peak_fit,1e-4,10e-9]
            
            sitposlist.append(fit_calculated_sitpos)     
            cs_gate(fit_calculated_sitpos)
            time.sleep(abs(fit_calculated_sitpos-stop_vgi)/slew_rate+1)
            current_csvg=copy.copy(fit_calculated_sitpos)
            
            closest_element_index = np.argmin(np.abs(cs_sweep - fit_calculated_sitpos))#find closest element in cs_sweep to sitpos
            sitpos_approximation_for_show=np.zeros(len(cs_sweep))
            sitpos_approximation_for_show[closest_element_index]=1.1*Glist[closest_element_index]

            datasaver_aux.add_result(('G', Glist),
                            ('V_aux', V_aux_list),
                            ('Phase_aux', Phase_aux_list),
                            ('approx_sitpos',sitpos_approximation_for_show),
                            (outer_gate_sweep.parameter,outer_gate_value),
                            (cs_sweep.parameter,cs_sweep_list))
            freq_rf(start_f-mix_down_f)#set to start point of frequency measurement
            gate_rf_enabled_param.value(1)#switch on gate rf
            for f_value in tqdm(freq_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
                freq_rf(f_value-freq_rlc())
                freq_mech(f_value)
                time.sleep(1.1*tc) # Wait 1.1 times the time contanst of the lock-in
                measured_value=measured_parameter()
                theta, v_r, I, G = zurich_phase_voltage_current_conductance(measured_value, source_amplitude_CNT_mech)
                I_list.append(I)
                V_list.append(v_r)
                Phase_list.append(theta)

            datasaver.add_result(('I_rf', I_list),
                                ('V_rf', V_list),
                                ('Phase', Phase_list),
                                (outer_gate_sweep.parameter,outer_gate_value),
                                (freq_sweep.parameter,freq_sweep_list))
            
            
        
#time.sleep(abs(stop_vg)/ramp_speed/1000 + 10)
gate_rf_enabled_param.value(0)
print("wake up, gates are")
print(cs_gate())


run_id = datasaver.run_id
foldername=f'C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_\\meas{run_id}'
if not os.path.exists(foldername):
    os.makedirs(foldername) 

#filename=f'meas{run_id}_sitpos_V.npy'
filename='sitpos_V.npy'
path = os.path.join(foldername, filename)
np.save(path, np.array(sitposlist))

#figfilename=f'meas{run_id}_GVgs_and_sitpoints.png'
figfilename='GVgs_and_sitpoints.png'
path = os.path.join(foldername, figfilename)
plt.savefig(path)  # Saves the plot as a PNG file
#plt.show()
#plt.close()  


