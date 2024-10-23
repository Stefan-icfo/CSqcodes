
#IN CONSTRUCTION

import numpy as np


from instruments import station, zurich, qdac, keithley2400, Triton
from qcodes.dataset import Measurement, new_experiment
from qcodes import Parameter
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
from utils.CS_utils import zurich_phase_voltage_current_conductance_compensate, save_metadata_var, get_var_name
from drivers.triton_utils import relevant_T
from drivers.Qdac_utils import ramp_time,read_channels

import drivers.k2400 as k2
import time
from tqdm import tqdm
import matplotlib.pyplot as plt
import scipy as scp
import copy
k2400=keithley2400

#------User input----------------
#gate basics
gate_ramp_slope = 1e-2# V/s

#rf source
tc_zurich = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd_dB = 0# zero for test 39 # attenuation at the source in dB
upperbound_lever_arm_rf_gate=0.2
#source_amplitude_instrumentlevel_GVg = 50e-3
#vsdac=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_GVg)-vsd_dB)
#compensation values
x_avg=8.9e-6
y_avg=-10.6e-6

#dc source
#vsd_dc = 200e-6
tc_keithley = 100e-3 
step_vdc_max = 200e-6 # source steps; for microvolts, one step is ok
offset_v = 27e-6 #voltage offset, to find zero point ie. if voltage offset has to be chosen 10uV to reach zero current then wrtie here +10uV
offset_i=-46e-12 #current measurement offset
upperbound_lever_arm_dc_gate=0.5

StartTemp=relevant_T()

#print(f"source amp at CNT for GVg:{vsdac*1e6} uV")
device_name = 'CD11_D7_C1'
prefix_name = 'Conductance_cooldown'
postfix = '__'
# exp_name = 'Test 50 K'

vars_to_save=[gate_ramp_slope,tc_zurich,vsd_dB,upperbound_lever_arm_rf_gate,tc_keithley,x_avg,y_avg,step_vdc_max,offset_v,offset_i,upperbound_lever_arm_dc_gate,StartTemp]

mix_down_f = 1.25e6 # RLC frequency
#####################
start_vgrf = 0#
stop_vgrf = 0.1#
#step_numrf =1001 #50uV
#step_vgrf=np.absolute((start_vgrf-stop_vgrf)/step_numrf)
#####################
start_vgdc = 0#
stop_vgdc = 0.1#

#temp_setpoints=list(np.linspace(300,50,num=6))
#temp_setpoints.extend(list(np.linspace(40,10,num=4)))
#temp_setpoints.extend(list(np.linspace(10,2,num=9)))
temp_setpoints=[300,200,150]
continuous_meas_at_low_T=True
#below that continuous measurement
irrelevant_difference=5e-3#change value below which no more extra measurement is needed


gate_rf=qdac.ch06#gate corresponding to RF source, eg charge sensor

gates_dc=[qdac.ch01,qdac.ch02,qdac.ch03,qdac.ch04,qdac.ch05] #gate corresponding to DC source, eg 5-gate section


freq = zurich.oscs.oscs0.freq
freq(mix_down_f)
#####################
gain_RT = 200       #
gain_HEMT = 5.64   #%
Z_tot = 7521        #
###################


#source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
#source_amplitude_param(source_amplitude_instrumentlevel_GVg)#set source aplitude
gate_rf.label = 'cs_gate' # Change the label of the rf gate chaneel
#instr_dict = dict(gate=[gate])
#exp_dict = dict(vsdac = vsdac)
#exp_name = sample_name(prefix_name,exp_dict,postfix)
exp_name_rf='cooldown_both_sides_rf'
exp_name_dc='cooldown_both_sides_dc'
#----------- defined values------



# ------------------Create a new Experiment-------------------------

#vgdc_sweep = gate.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here
#rf_gate_sweep = gate_rf.dc_constant_V.sweep(start=start_vgrf, stop=stop_vgrf, num = step_numrf)
measured_parameter1 = zurich.demods.demods0.sample
source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
measured_parameter2 = k2400.curr
source = k2400
#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT
# applied  voltages at the intrument level before attenuation
#vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #what is this?
#bilt.channels.output_mode('ramp')
#bilt.channels.ramp_slope(ramp_speed)


gate_rf.dc_slew_rate_V_per_s(gate_ramp_slope)
gate_rf.dc_constant_V(start_vgrf)

for gate_dc in gates_dc:
    gate_dc.dc_slew_rate_V_per_s(gate_ramp_slope)
    gate_dc.dc_constant_V(start_vgdc)
all_gates=copy.copy(gates_dc)
all_gates.append(gate_rf)
setpoints_dc=[start_vgdc,start_vgdc,start_vgdc,start_vgdc,start_vgdc,0]#for convencience
setpoints_rf=[0,0,0,0,0,start_vgrf]
setpoints=[start_vgdc,start_vgdc,start_vgdc,start_vgdc,start_vgdc,start_vgrf]

#initial ramp
print(f"sleep for {ramp_time(all_gates,setpoints)+3} seconds to ramp the gates")
time.sleep(ramp_time(all_gates,setpoints)+3)
print("wake up, gates are")
read_channels()
#T=relevant_T()
Temp_param=Parameter('T', label='T', unit='K',
                       get_cmd=lambda: Temp_now)

Vgrf_param = Parameter('Vgrf', label='Vgrf', unit='V',
                       get_cmd=lambda: Vgrf_now)



# ----------------Create a measurement-------------------------
experiment_rf = new_experiment(name=exp_name_rf, sample_name=device_name)
meas1 = Measurement(exp=experiment_rf)
meas1.register_parameter(Vgrf_param)  # register1the 1st1indepe n 10e-3ent parameter
meas1.register_parameter(Temp_param)
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas1.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[Vgrf_param,Temp_param])
meas1.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[Vgrf_param,Temp_param])
meas1.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[Vgrf_param,Temp_param])
meas1.register_custom_parameter('R', 'R', unit='Ohm', basis=[], setpoints=[Vgrf_param,Temp_param])



Vgdc_param = Parameter('Vgdc', label='Vgdc', unit='V',
                       get_cmd=lambda: Vgdc_now)


experiment_dc = new_experiment(name=exp_name_dc, sample_name=device_name)
meas2 = Measurement(exp=experiment_dc)
meas2.register_parameter(Temp_param)
meas2.register_parameter(Vgdc_param)  # register1the 1st1indepe n 10e-3ent parameter
meas2.register_custom_parameter('Conductance_IV', 'G_IV', unit='Siemens', basis=[], setpoints=[Vgdc_param,Temp_param])
meas2.register_custom_parameter('Current_top', 'I_top', unit='Amps', basis=[], setpoints=[Vgdc_param,Temp_param])
meas2.register_custom_parameter('Current_zero', 'I_zero', unit='Amps', basis=[], setpoints=[Vgdc_param,Temp_param])
meas2.register_custom_parameter('Current_bottom', 'I_bottom', unit='Siemens', basis=[], setpoints=[Vgdc_param,Temp_param])
meas2.register_custom_parameter('Resistance_IV', 'R_IV', unit='Ohm', basis=[], setpoints=[Vgdc_param,Temp_param])


# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas1.run() as datasaver1:
    datasaver1.dataset.add_metadata('qdac_ch01_dc_constant_V_start',qdac.ch01.dc_constant_V())
    datasaver1.dataset.add_metadata('qdac_ch02_dc_constant_V_start',qdac.ch02.dc_constant_V())
    datasaver1.dataset.add_metadata('qdac_ch03_dc_constant_V_start',qdac.ch03.dc_constant_V())
    datasaver1.dataset.add_metadata('qdac_ch04_dc_constant_V_start',qdac.ch04.dc_constant_V())
    datasaver1.dataset.add_metadata('qdac_ch05_dc_constant_V_start',qdac.ch05.dc_constant_V())
    datasaver1.dataset.add_metadata('qdac_ch06_dc_constant_V_start',qdac.ch06.dc_constant_V())
    varnames=[]
    for i in range(len(vars_to_save)):
        varnames.append(get_var_name(vars_to_save[i]))
    save_metadata_var(datasaver1.dataset,varnames,vars_to_save)
    # for i in range(2):
    with meas2.run() as datasaver2:

        #current_temp=relevant_T()
        next_temp_setpoint=temp_setpoints[0]

        for current_T_setpoint in temp_setpoints:
           #start by ramping gates to start value, others to zero 
            gate_rf.dc_constant_V(start_vgrf)
            for gate_dc in gates_dc:
                gate_dc.dc_constant_V(0)
            print(f"sleep for {ramp_time(all_gates,setpoints_rf)+1} seconds to ramp the dc gates")
            time.sleep(ramp_time(all_gates,setpoints_rf)+1)
 

            while relevant_T()>current_T_setpoint:
                if waiting==False:
                    print(f"waiting for temp to decrease to {current_T_setpoint} T")
                waiting=True
                time.sleep(10)
            waiting=False

            ####################################
            current_temp=temp_setpoints[0]#setpoint for test relevant_T()
            vsdkTrf=current_temp/11604
            vsdkTrf_instrument=d2v(v2d(np.sqrt(2)*vsdkTrf)+vsd_dB) 
            
            source_amplitude_param(vsdkTrf_instrument)##############FIX THIS FIX THIS FIX THIS FIX THIS!!!!
            step_numrf = round((stop_vgrf-start_vgrf)/vsdkTrf*upperbound_lever_arm_rf_gate)+1  
            step_rf=abs(stop_vgrf-start_vgrf)/(step_numrf-1)

            ####FOR INITIAL TEST HERE OVERRIDE STEP NUMBER(WITH CONSTANT VALUE) AND CURRENT TEMP (WITH SETPOINT)
            vgrf_sweep = gate_rf.dc_constant_V.sweep(start=start_vgrf, stop=stop_vgrf, num = step_numrf)
            ###############################
            print(f"step_numrf={step_numrf}")
            for vgrf_value in tqdm(vgrf_sweep, leave=False, desc='csgate voltage Sweep', colour = 'green'):
                
                gate_rf.dc_constant_V(vgrf_value)
                time.sleep(1.1*tc_zurich+step_rf/gate_ramp_slope) # Wait 3 times the time contanst of the lock-in plus gate ramp speed
                measured_value1 = measured_parameter1()
                theta_calc, v_r_calc, I,  G = zurich_phase_voltage_current_conductance_compensate(measured_value1, vsdkTrf,x_avg, y_avg)
                R = 1/G
                #NOT SURE THIS CAN HANDLE THE DIFFERENT SIZED OF THE ARRAYS...MIGHT HAVE TO DEFINE SOME FINEGRAINED ARRAY FOR ALL AND THEN FILL IN
                ###################for test replaced current T with setpoint
                datasaver1.add_result(('R', R),
                                    ('G', G),
                                    ('V_r', v_r_calc),
                                    ('Phase', theta_calc),
                                    (Vgrf_param ,vgrf_value),(Temp_param,current_T_setpoint))
                

             #now re-ramp gates in new configuration and perform second measurement loop. Maybe...test first!

            source_amplitude_param(0)#switch off RF source for other side'd DC GVg
             
            gate_rf.dc_constant_V(0)
            for gate_dc in gates_dc:
                gate_dc.dc_constant_V(start_vgdc)
            print(f"sleep for {ramp_time(all_gates,setpoints_dc)+1} seconds to ramp the gates")
            time.sleep(ramp_time(all_gates,setpoints_dc)+1)
            ###############################
            current_temp=current_T_setpoint#setpoint for test relevant_T()
            vsdkTdc=current_temp/11604
            if vsdkTdc>step_vdc_max:
                divider=vsdkTdc/step_vdc_max
                step_v=vsdkTdc/round(divider)#better: round down!
            else:
                step_v=vsdkTdc
            step_numdc = round((stop_vgdc-start_vgdc)/vsdkTdc*upperbound_lever_arm_dc_gate)+1  
            print(f"step_numdc={step_numdc}")
            step_dc=abs(stop_vgdc-start_vgdc)/(step_numdc-1)
            gate1=gates_dc[0]

            ####FOR INITIAL TEST HERE OVERRIDE STEP NUMBER(WITH CONSTANT VALUE) AND CURRENT TEMP (WITH SETPOINT)
            vgdc_sweep = gate1.dc_constant_V.sweep(start=start_vgdc, stop=stop_vgdc, num = step_numdc)


            for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='Gate Sweep', colour = 'green'):
        
                for gate in gates_dc:
                    gate.dc_constant_V(vgdc_value)
                time.sleep(step_dc/gate_ramp_slope)
                #measure lower vsd value (~kT)
                k2.ramp_k2400(ramp_param=source,final_vg=vsdkTdc+offset_v, step_size = step_v, ramp_speed=1)
                #time.sleep(tc_keithley) #
                measured_value2 = measured_parameter2()
                I_bottom=measured_value2-offset_i+1e-15
                R_bottom = -vsdkTdc/I_bottom

                #measure effective zero value 
                k2.ramp_k2400(ramp_param=source,final_vg=offset_v, step_size = step_v, ramp_speed=1)
                #time.sleep(2*ts) # Wait 2 times the time constant for the source to settle
                measured_value = measured_parameter2()
                I_zero=measured_value2-offset_i+1e-15

            #measure upper vsd value (~ + kT) 
                k2.ramp_k2400(ramp_param=source,final_vg=vsdkTdc+offset_v, step_size = step_v, ramp_speed=1)
                #time.sleep(2*ts) # Wait 2 times the time constant for the source to settle
                measured_value = measured_parameter2()
                I_top=measured_value2-offset_i+1e-15
                R_top = vsdkTdc/I_top

            #linear regression fit
                x=[-vsdkTdc,0,vsdkTdc]
                y=[I_bottom,I_zero,I_top]
                result=scp.stats.linregress(x,y)
            
            
                G_IV=result.slope
                R_IV=1/G_IV


                #pt.plot(x,y,'bo')#just to see individual IV fits for test of code
                #pt.show()
                ###################for test replaced current T with setpoint
                datasaver2.add_result(('Conductance_IV', G_IV),('Resistance_IV', R_IV),('Current_top', I_top),('Current_zero', I_zero),('Current_bottom', I_bottom),(Vgdc_param, vgdc_value),(Temp_param,current_T_setpoint))
            k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=1)
        # vgdc_sweep.reverse()
        # Ramp down everything
   ####################################################
   ##ONCE THIS WORKS:
   # -1)CONVERT ZURICH AMP TO AMP AT INSTRUMENT LEVEL; AFTER THAT COULD BE READY TO TEST --- should be done!
   # 0) make sure data arrays fit and save properly; 
   # 1)Take care of continuous continuation of loop below cutoff temp (modify while statement at beginning)
   # 2)plot both streams as lines on plot
   # 3)small things to fix: switch off output of zurich while keithley's on, ramp down the gates

