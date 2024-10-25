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
import os
from qcodes import Parameter
import copy
import math

script_path = __file__
print("Full path of the script:", script_path)

debug=False
#------User input----------------
slew_rate=1e-2

tc = 100e-3   # in seconds.
att_source_dB = 39 # attenuation at the source in dB
att_gate_dB =46
#vsdac = 200e-6 # source AC voltage in volt
device_name = 'CD11_D7_C1'
#device_name =  'CD05_G6_E3_'# 
prefix_name = '_cs_mechanics_detune_'#

postfix = '20mK'

#additional_diagonal_detuning=100e-3

#Temp=Triton.MC()
#postfix = f"{Temp}K"
#vsdkT=Temp/11604
#vsd=vsdkT
#compensation values
x_avg=+4.38e-6
y_avg=-4.41e-6

mix_down_f = 1.25e6 # RLC frequency
#outer gate voltage range (slow axis, 5gate)
#####################
#idt_point1_x=-1.96894
#idt_point1_y=-1.97167
#idt_point2_x=-1.96544
#idt_point2_y=-1.96808
delta_simple=5e-3#45deg
middle_point_g2=-1.2674
middle_point_g4=-1.1818
#angle=3.14/4#45deg

step_vgo_num = 5#
#start_vgo2,start_vgo1,stop_vgo2,stop_vgo1=make_detuning_axis(idt_point1_x,idt_point1_y,idt_point2_x,idt_point2_y,delta) 
start_vgo1=middle_point_g2-delta_simple
stop_vgo1=middle_point_g2+delta_simple
start_vgo2=middle_point_g4+delta_simple
stop_vgo2=middle_point_g4-delta_simple
#reverse
#stop_vgo1=middle_point_g2-delta_simple
#start_vgo1=middle_point_g2+delta_simple
#stop_vgo2=middle_point_g4+delta_simple
#start_vgo2=middle_point_g4-delta_simple
source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value
postfix = f"_{round(gate_amplitude_param()*1000,3)}mV on gate@inst,_{round(source_amplitude_param()*1000,3)}mV on source@inst, g1={round(qdac.ch01.dc_constant_V(),2)},g2mp={middle_point_g2},g3={round(qdac.ch03.dc_constant_V(),2)},g4mp={middle_point_g2},g5={round(qdac.ch05.dc_constant_V(),2)},gcs={round(qdac.ch06.dc_constant_V(),2)}"




print(start_vgo2,start_vgo1,stop_vgo2,stop_vgo1)

step_vgo1=np.absolute((start_vgo1-stop_vgo1)/step_vgo_num)
step_vgo2=np.absolute((start_vgo2-stop_vgo2)/step_vgo_num)

vars_to_save=[slew_rate,tc,att_source_dB,att_gate_dB,x_avg,y_avg,mix_down_f,middle_point_g2,middle_point_g4,delta_simple,step_vgo_num]

#inner gate voltage range (fast axis, CS)
#####################
start_vgi = -1.457#-0.788
stop_vgi = -1.453#-0.776
step_vgi_num = 4*20#40uV
#step_vgi_num = round((stop_vgi-start_vgi)/vsd*upper_bound_lever_arm)
#print(f"step i num={step_vgi_num}")
step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)

initial_guess = [-1.455, 1e-4, 3e-7]#initial guess for peakV, Gamma,height for first GVg
sitfraction=0.6#where to sit on Coulomb peak. For now on left side

vars_to_save.extend([start_vgi,stop_vgi,step_vgi_num])
#####################
start_f =100e6#Hz unit
stop_f = 150e6 #Hz unit
step_num_f = 50*1000#2kHz

vars_to_save.extend([start_f,stop_f,step_num_f])

source_amplitude_instrumentlevel_GVg = 20e-3
source_amplitude_CNT_GVg=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_GVg)-att_source_dB)
print(f"source amp at CNT for GVg:{source_amplitude_CNT_GVg*1e6} uV")
source_amplitude_instrumentlevel_mech = 20e-3
source_amplitude_CNT_mech=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_mech)-att_source_dB)
print(f"source amp at CNT for mech:{source_amplitude_CNT_mech*1e6} uV")
gate_amplitude_instrumentlevel = 20e-6
gate_amplitude_CNT=d2v(v2d(np.sqrt(1/2)*gate_amplitude_instrumentlevel)-att_gate_dB)
print(f"gate amp at CNT for mech:{gate_amplitude_CNT*1e6} uV")

vars_to_save.extend([source_amplitude_instrumentlevel_GVg, source_amplitude_CNT_GVg, source_amplitude_instrumentlevel_mech,source_amplitude_CNT_mech,gate_amplitude_instrumentlevel,gate_amplitude_CNT])

#--------Definitions-------------

#swept contacts
cs_gate=qdac.ch06.dc_constant_V  # swept gate voltage
qdac.ch06.dc_slew_rate_V_per_s(slew_rate)

outer_gate1=qdac.ch02.dc_constant_V
qdac.ch02.dc_slew_rate_V_per_s(slew_rate)

outer_gate2=qdac.ch04.dc_constant_V
qdac.ch04.dc_slew_rate_V_per_s(slew_rate)


#constant gate voltages, labelled by the channels they are connected to; 
#gate_V_ch3=+1
#gate_V_ch1=-3
#gate_V_ch5=-3

#initialize constant gates, comment out for single-gate device

#qdac.ch03.dc_slew_rate_V_per_s(slew_rate)
#qdac.ch03.dc_constant_V(gate_V_ch3)
#qdac.ch05.dc_slew_rate_V_per_s(slew_rate)
#qdac.ch05.dc_constant_V(gate_V_ch5)
#qdac.ch01.dc_slew_rate_V_per_s(slew_rate)
#qdac.ch01.dc_constant_V(gate_V_ch1)



outer_gate1(start_vgo1)
outer_gate2(start_vgo2)
cs_gate(start_vgi)
print('wait time')
#time.sleep(10)
sleeptime=max(abs(start_vgo1-outer_gate1()),abs(start_vgo2-outer_gate2()),abs(start_vgi-cs_gate()))/slew_rate+10
print(sleeptime)
time.sleep(sleeptime)
print("wake up, gates are")
print(outer_gate1())
#print(outer_auxgate1())
print(outer_gate2())
print(cs_gate())






#outer_gate1.label = '5g(outer)' # Change the label of the gate chanel
cs_gate.label = 'CS_gate(inner)' # Change the label of the source chaneel
#instr_dict = dict(gate=[outer_gate1])
exp_dict = dict(uVrfgate = gate_amplitude_CNT*1e6)
exp_name = sample_name(prefix_name,exp_dict,postfix)
#----------- defined values------
#----------- defined values------
#####################
#gain_RT = 200       #
#gain_HEMT = 5.64   #
#Z_tot = 7521        #
###################

# ------------------define sweep axes-------------------------

outer_gate1_sweep=outer_gate1.sweep(start=start_vgo1, stop=stop_vgo1, num = step_vgo_num)
outer_gate2_sweep=outer_gate2.sweep(start=start_vgo2, stop=stop_vgo2, num = step_vgo_num)
outer_gate1_list=list(outer_gate1_sweep)
outer_gate2_list=list(outer_gate2_sweep)

g1_array=np.array(outer_gate1_list)
g2_array=np.array(outer_gate2_list)
delta_array=np.sqrt(abs((g1_array-start_vgo1)**2+(g2_array-start_vgo2)**2))
delta_array-=delta_simple*math.sqrt(2)

cs_sweep=cs_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)

freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq
freq_rlc(mix_down_f)
freq_mech(start_f)
freq_rf(start_f-mix_down_f)
time.sleep(tc)

freq_sweep = freq_rf.sweep(start=start_f, stop=stop_f, num = step_num_f)
freq_sweep_list=list(freq_sweep)
measured_parameter = zurich.demods.demods2.sample  
measured_aux_parameter = zurich.demods.demods0.sample

source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value
gate_rf_enabled_param = getattr(zurich.sigouts.sigouts1.enables, f'enables{1}')

gate_amplitude_param(gate_amplitude_instrumentlevel)

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT
#vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
# applied  voltages at the intrument level before attenuation
#vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #SF:this sets up the zurich LIA

delta_current=0#just to define param
delta_param = Parameter('delta', label='delta', unit='V',
                       get_cmd=lambda: delta_now)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(delta_param)
meas.register_parameter(freq_sweep.parameter)
#meas.register_parameter(inner_gate_sweep.parameter)   # 
meas.register_custom_parameter('V_rf', 'Amplitude', unit='V', basis=[], setpoints=[delta_param,freq_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[delta_param,freq_sweep.parameter])
meas.register_custom_parameter('I_rf', 'current', unit='I', basis=[], setpoints=[delta_param,freq_sweep.parameter])
#meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])

experiment_aux = new_experiment(name=exp_name+"aux", sample_name=device_name)
meas_aux = Measurement(exp=experiment_aux)
meas_aux.register_parameter(delta_param)  # 
meas_aux.register_parameter(cs_sweep.parameter)
meas_aux.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[delta_param,cs_sweep.parameter])
meas_aux.register_custom_parameter('V_aux', 'Amplitude_aux', unit='V', basis=[], setpoints=[delta_param,cs_sweep.parameter])
meas_aux.register_custom_parameter('Phase_aux', 'Phase_aux', unit='rad', basis=[], setpoints=[delta_param,cs_sweep.parameter])
meas_aux.register_custom_parameter('outer_gate1', 'outer_gate1', unit='V', basis=[], setpoints=[delta_param,cs_sweep.parameter])
meas_aux.register_custom_parameter('outer_gate2', 'outer_gate2', unit='V', basis=[], setpoints=[delta_param,cs_sweep.parameter])

# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:
    #saving metadata parameters
    datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V',qdac.ch01.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V',qdac.ch03.dc_constant_V())
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
        i=0 #outer sweep counter
        First_run=True
        sitposlist=[]
        current_csvg=copy.copy(start_vgi)
        if debug:
            print("starting loops")
        for outer_gate_value in tqdm(outer_gate1_sweep, leave=False, desc='outer Gate Sweep', colour = 'green'): #slow axis loop (gate)
            i=i+1#outergatesweepcounter
            #print('temperature')
            #Triton.MC()
            outer_gate1_sweep.set(outer_gate_value)
            outer_gate2_sweep.set(outer_gate2_sweep[i-1])
            time.sleep(max(abs(step_vgo1/slew_rate),abs(step_vgo2/slew_rate))) # Wait  the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
            #lists for GVgs
            G_list=[]
            V_aux_list=[]
            Phase_aux_list=[]

            #lists for frequency sweeps
            I_list=[]
            V_list=[]
            Phase_list=[]

            #run GVg to find sitpos
            freq_rf(mix_down_f)#set to RLC frequency
            #########################
           # gate_rf_enabled_param.value(0)#switch off gate output
            source_amplitude_param(source_amplitude_instrumentlevel_GVg)#set source aplitude
            #########################
            cs_gate(start_vgi)
            time.sleep(abs(start_vgi-current_csvg)/slew_rate)

            if debug:
                print("starting GVg")
            
            #while zurich_working(measured_aux_parameter)==False:
            #    print("zurich not working, waiting 10 min")
            #    time.sleep(600)
            for gatecs_value in tqdm(cs_sweep, leave=False, desc='cs Sweep', colour = 'cyan'):
                cs_sweep.set(gatecs_value)
                time.sleep(1.1*tc+abs(step_vgi)/slew_rate)
                measured_value = measured_aux_parameter()
                theta_aux, v_aux, I_aux, G = zurich_phase_voltage_current_conductance_compensate(measured_value, source_amplitude_CNT_GVg, x_avg,y_avg)
                G_list.append(G)
                V_aux_list.append(v_aux)
                Phase_aux_list.append(theta_aux)

            #initial_guess=[peak_fit,1e-4,10e-9]
            popt, pcov = scp.optimize.curve_fit(breit_wigner_fkt, cs_sweep_list, G_list, p0=initial_guess)
            peak_fit, hgamma_fit, peak_G_fit=popt
            fit_detuning=-breit_wigner_detuning(peak_G_fit*sitfraction,peak_G_fit,hgamma_fit)
            fit_calculated_sitpos=peak_fit+fit_detuning
            sitposlist.append(fit_calculated_sitpos)

            #set cs gate to sitpos
            if (start_vgi-0.1) < fit_calculated_sitpos < (start_vgi+0.1):
                cs_gate(fit_calculated_sitpos)
            time.sleep(abs(fit_calculated_sitpos-stop_vgi)/slew_rate+1)
            current_csvg=copy.copy(fit_calculated_sitpos)

            if debug:
                plt.figure(1)
                plt.title("GVg, fit, and sitpos")
                plt.plot(cs_sweep_list,G_list)
                plt.plot(cs_sweep_list,breit_wigner_fkt(cs_sweep_list,peak_fit,hgamma_fit,peak_G_fit))
                plt.plot(current_csvg,breit_wigner_fkt(current_csvg,peak_fit,hgamma_fit,peak_G_fit),'go')
                plt.show()
            else:
                plt.figure(1)
                plt.title("GVg, fit, and sitpos")
                plt.plot(cs_sweep_list,G_list)
                plt.plot(cs_sweep_list,breit_wigner_fkt(cs_sweep_list,peak_fit,hgamma_fit,peak_G_fit))
                plt.plot(current_csvg,breit_wigner_fkt(current_csvg,peak_fit,hgamma_fit,peak_G_fit),'go')


            datasaver_aux.add_result(('G', G_list),
                            ('V_aux', V_aux_list),
                            ('Phase_aux', Phase_aux_list),
                            ('outer_gate1', [outer_gate1_sweep[i-1]]*len(cs_sweep_list)),
                            ('outer_gate2', [outer_gate2_sweep[i-1]]*len(cs_sweep_list)),
                            (delta_param,delta_array[i-1]),
                            (cs_sweep.parameter,cs_sweep_list))

            freq_rf(start_f-mix_down_f)#set to start point of frequency measurement
            gate_rf_enabled_param.value(1)#switch on gate rf
            source_amplitude_param(source_amplitude_instrumentlevel_mech)#set source aplitude
            #gate_amplitude_param(gate_amplitude_instrumentlevel/2**(i-1))
            #if debug:
            #print(f"starting frequency sweep with gaterf={gate_amplitude_param()}")
            #while
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
                                (delta_param,delta_array[i-1]),
                                (freq_sweep.parameter,freq_sweep_list))
            
            if First_run==False: #ie if it's not the first run and popt has been measured, then redefine the initial guess by using the last fitted values
                initial_guess = [popt[0], popt[1], popt[2]]
            First_run=False
        
#time.sleep(abs(stop_vg)/ramp_speed/1000 + 10)
gate_rf_enabled_param.value(0)
print("wake up, gates are")
print(outer_gate1())
print(outer_gate2())
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


