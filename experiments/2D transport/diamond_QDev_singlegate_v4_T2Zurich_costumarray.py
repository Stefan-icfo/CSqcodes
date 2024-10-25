# Sweeps gate and source voltages while reading out with Zurich LIA 
# Sweeps only single gate while keeping other gates constant
# Stefan Forstner using template by Parameshwar Prasad
# snake is working now. watch out - there might still be issues when setting step_ramp_speed to slow and the step_size too large!
# uses triton 1 Zurich

import numpy as np
import os

from instruments import  manual, station, qdac, zurich, Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk

import time
from tqdm import tqdm
import copy

from qcodes import Parameter



#------User input----------------
ramp_speed_gate = 0.01 # V/s for large ramps
ramp_speed_source = 0.01 # V/s for large ramps
step_ramp_speed=0.01 # between steps, V/s
tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd_dB = 39 # attenuation at the source in dB
vsdac = 158e-6/10 # source AC voltage in volt
device_name = 'CD11_D7_C1_cs'
prefix_name = 'Diamond_costumsweep'
postfix = '___'
# exp_name = 'Test 50 K'
Temp=Triton.MC()
postfix = f"_{Temp}K,g1={round(qdac.ch01.dc_constant_V(),2)},g2={round(qdac.ch02.dc_constant_V(),2)},g3={round(qdac.ch03.dc_constant_V(),2)},g4={round(qdac.ch04.dc_constant_V(),2)},g5={round(qdac.ch05.dc_constant_V(),2)}"
mix_down_f = 1.25e6 #

#gate voltage range (slow axis)
#####################

start_vg = -2.5    #
stop_vg = -0   #
step_vg_num = 2500    #
step_vg=np.absolute((start_vg-stop_vg)/step_vg_num)


#source voltage range (fast axis)
#####################
source_sweep=list(np.linspace(-30,-4,26)/1e3)
source_sweep.extend(list(np.linspace(-3,3,6*2+1)/1e3))
source_sweep.extend(list(np.linspace(4,30,26)/1e3))
start_vs=source_sweep[0]
    #
#step_vs=np.absolute((start_vs-stop_vs)/step_vs_num)

#constant gate voltages, labelled by the channels they are connected to; 


#initialize constant gates, comment out for single-gate device

# qdac.ch03.dc_slew_rate_V_per_s(ramp_speed)
# qdac.ch03.dc_constant_V(gate_V_ch3)
# qdac.ch04.dc_slew_rate_V_per_s(ramp_speed)
# qdac.ch04.dc_constant_V(gate_V_ch4)
# qdac.ch05.dc_slew_rate_V_per_s(ramp_speed)
# qdac.ch05.dc_constant_V(gate_V_ch5)
# qdac.ch06.dc_slew_rate_V_per_s(ramp_speed)
# qdac.ch06.dc_constant_V(gate_V_ch6)
# qdac.ch07.dc_slew_rate_V_per_s(ramp_speed)
# qdac.ch07.dc_constant_V(gate_V_ch7)


#--------Definitions-------------

#swept contacts
gate=qdac.ch06  # swept gate voltage
source=qdac.ch07 #swept source voltage


freq = zurich.oscs.oscs0.freq
gate.label = 'gate_label' # Change the label of the gate chanel
source.label = 'source_label' # Change the label of the source chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = vsdac)
exp_name = sample_name(prefix_name,exp_dict,postfix)

#----------- defined values------
#####################
gain_RT = 200       #
gain_HEMT = 5.64    #
Z_tot = 7521        #
#####################


# ------------------define sweep axes-------------------------

gate_sweep=gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_vg_num)
#source_sweep=source.dc_constant_V.sweep(start=start_vs, stop=stop_vs, num = step_vs_num)

#source_logsweep=[]




measured_parameter = zurich.demods.demods0.sample   # lock-in amplitude measured # SF:CHANGED FROM 'demod_complex'

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT

# applied  voltages at the intrument level before attenuation
vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #SF:this sets up the zurich LIA



#initialize swept contacts

#slow ramp and intial voltage
gate.dc_slew_rate_V_per_s(ramp_speed_gate)
gate.dc_constant_V(start_vg)

source.dc_slew_rate_V_per_s(ramp_speed_source)
source.dc_constant_V(start_vs)



time.sleep(max([abs((start_vg-gate.dc_constant_V())/ramp_speed_gate),abs((source_sweep[0]-source.dc_constant_V())/ramp_speed_source)])+1)  #wait for the time it takes to do both ramps plus one second

last_source_value=source.dc_constant_V()

source_sweep_param = Parameter('source_sweep_param',
                            label='Source Sweep',
                            unit='V',  # Or the appropriate unit
                            set_cmd=None,  # If there is no instrument to communicate with directly
                            get_cmd=None)  # Define get_cmd if you need to read a value

source_sweep_param.set(last_source_value)
#set fast ramp speeds
#gate.dc_slew_rate_V_per_s(step_ramp_speed)
#source.dc_slew_rate_V_per_s(step_ramp_speed)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(gate_sweep.parameter)  # 
meas.register_parameter(source_sweep_param)  # 
meas.register_custom_parameter('R', 'R', unit='Ohm', basis=[], setpoints=[gate_sweep.parameter,source_sweep_param])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[gate_sweep.parameter,source_sweep_param])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[gate_sweep.parameter,source_sweep_param])
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[gate_sweep.parameter,source_sweep_param])


# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:

    fast_axis_unreversible_list = list(source_sweep) #(to deal with snake)
    reversed_sweep=False

    for gate_value in tqdm(gate_sweep, leave=False, desc='Gate Sweep', colour = 'green'): #slow axis loop (gate)
        gate_sweep.set(gate_value)
        time.sleep(1.1*tc+step_vg/step_ramp_speed) # Wait 3 times the time contanst of the lock-in plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
        #init lists (to deal with snake)
        Ilist=[]
        Rlist=[]
        Glist=[]
        VRlist=[]
        PHASElist=[]
        # temp_fast_axis_list=[] #I think this line is unnecessary
        #following is necessary only if no snake #COMMENT OUT ONCE SNAKE IS WORKING
        #source.dc_slew_rate_V_per_s(ramp_speed)
        #source.dc_constant_V(start_vs)
        #time.sleep(abs((start_vs-stop_vs)/ramp_speed)) 
        #source.dc_slew_rate_V_per_s(step_ramp_speed)

        for source_value in source_sweep:
            source.dc_constant_V(source_value)
            source_sweep_param.set(source_value)
            step_vs=abs(source_value-last_source_value)
            last_source_value=copy.copy(source_value)

            time.sleep(1.1*tc+step_vs/step_ramp_speed) # Wait 3 times the time constant of the lock-in, plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
            measured_value = measured_parameter()
            x = measured_value['x'][0] #SF: COMMENTED OUT 
            y = measured_value['y'][0]#SF: COMMENTED OUT
            xy_complex =complex(x,y)#measured_value
            v_r_calc = np.absolute(xy_complex)
            theta_calc = np.angle(xy_complex)
                    
            #G calculation
            I = v_r_calc/(gain_RT*gain_HEMT*Z_tot)
            G = 1/((vsdac/I)-Z_tot)-12e-9
            R = 1/G

            #add to lists (to deal with snake)
            Ilist=Ilist+[I]
            Rlist=Rlist+[R]
            Glist=Glist+[G]
            VRlist=VRlist+[v_r_calc]
            PHASElist=PHASElist+[theta_calc]
            
        #Rlist.reverse
        
        #temp_fast_axis_list.reverse()
        if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Rlist.reverse()
            Glist.reverse()
            VRlist.reverse()
            PHASElist.reverse()
        datasaver.add_result(('R',Rlist),
                            ('G',Glist),
                            ('V_r',VRlist),
                            ('Phase',PHASElist),
                            (gate_sweep.parameter,gate_value),
                            (source_sweep_param,fast_axis_unreversible_list))
        source_sweep.reverse() 
        reversed_sweep= not reversed_sweep

# Ramp down everything
gate.dc_slew_rate_V_per_s(ramp_speed_gate)
source.dc_slew_rate_V_per_s(ramp_speed_source)


#gate.dc_constant_V(0)
source.dc_constant_V(0)


