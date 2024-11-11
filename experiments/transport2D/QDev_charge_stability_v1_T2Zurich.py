# Sweeps gate voltages while reading out with Zurich LIA 
# Sweeps two gates while keeping others constant
# Stefan Forstner using template by Parameshwar Prasad

# snake is working now. watch out - there might still be issues when setting step_ramp_speed to slow and the step_size too large!

import numpy as np
import os

from instruments import  manual, station, qdac, zurich
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk

import time
from tqdm import tqdm


#------User input----------------
ramp_speed = 1 # V/s for large ramps
step_ramp_speed=2 # between steps, V/s
tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd_dB = 60 # attenuation at the source in dB
vsdac = 100e-6 # source AC voltage in volt
device_name = 'device_name_xxX'
prefix_name = 'Charge_stability__QDev'
postfix = 'postfix'
# exp_name = 'Test 50 K'

mix_down_f = 1.25e6 #

#outer voltage range (slow axis)
#####################
start_vg1 = -2     #
stop_vg1 = 2       #
step_vg1_num = 11     #
step_vg1=np.absolute((start_vg1-stop_vg1)/step_vg1_num)


#inner voltage range (fast axis)
#####################
start_vg2 = -3     #
stop_vg2 = 3       #
step_vg2_num = 11    #
step_vg2=np.absolute((start_vg2-stop_vg2)/step_vg2_num)

#constant gate voltages, labelled by the channels they are connected to; 
gate_V_ch3=0
gate_V_ch4=0
gate_V_ch5=0
gate_V_ch6=0
gate_V_ch7=0

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
gate1=qdac.ch01  # swept outer gate voltage
gate2=qdac.ch02 #swept inner gate voltage


freq = zurich.oscs.oscs1.freq
gate1.label = 'gate1_label' # Change the label of the gate1 chanel
gate2.label = 'gate2_label' # Change the label of the gate2 chaneel
instr_dict = dict(gate1=[gate1])
exp_dict = dict(vsdac = vsdac)
exp_name = sample_name(prefix_name,exp_dict,postfix)

#----------- defined values------
#####################
gain_RT = 200       #
gain_HEMT = 5.64    #
Z_tot = 7521        #
#####################


# ------------------define sweep axes-------------------------

gate1_sweep=gate1.dc_constant_V.sweep(start=start_vg1, stop=stop_vg1, num = step_vg2_num)
gate2_sweep=gate2.dc_constant_V.sweep(start=start_vg2, stop=stop_vg2, num = step_vg2_num)
measured_parameter = zurich.source.voltage   # lock-in amplitude measured # SF:CHANGED FROM 'zurich.source.demod_complex'

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT

# applied  voltages at the intrument level before attenuation
vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #SF:this sets up the zurich LIA



#initialize swept contacts

#slow ramp and intial voltage
gate1.dc_slew_rate_V_per_s(ramp_speed)
gate1.dc_constant_V(start_vg1)

gate2.dc_slew_rate_V_per_s(ramp_speed)
gate2.dc_constant_V(start_vg2)

time.sleep(max([abs(start_vg1/ramp_speed),abs(start_vg2/ramp_speed)])+1)  #wait for the time it takes to do both ramps plus one second

#set fast ramp speeds
gate1.dc_slew_rate_V_per_s(step_ramp_speed)
gate2.dc_slew_rate_V_per_s(step_ramp_speed)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(gate1_sweep.parameter)  # 
meas.register_parameter(gate2_sweep.parameter)  # 
meas.register_custom_parameter('R', 'R', unit='Ohm', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])


# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:

    fast_axis_unreversible_list = list(gate2_sweep) #(to deal with snake)
    reversed_sweep=False
    
    for gate1_value in tqdm(gate1_sweep, leave=False, desc='outer gate sweep', colour = 'green'): #slow axis loop (gate)
        gate1_sweep.set(gate1_value)
        time.sleep(3*tc+step_vg1/step_ramp_speed) # Wait 3 times the time contanst of the lock-in plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
        #init lists (to deal with snake)
        Ilist=[]
        Rlist=[]
        Glist=[]
        VRlist=[]
        PHASElist=[]
        # temp_fast_axis_list=[] #I think this line is unnecessary
        #following is necessary only if no snake #COMMENT OUT ONCE SNAKE IS WORKING
        #gate2.dc_slew_rate_V_per_s(ramp_speed)
        #gate2.dc_constant_V(start_vg2)
        #time.sleep(abs((start_vg2-stop_vg2)/ramp_speed)) 
        #gate2.dc_slew_rate_V_per_s(step_ramp_speed)

        for gate2_value in tqdm(gate2_sweep, leave=False, desc='inner gate sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
            gate2_sweep.set(gate2_value)
            time.sleep(3*tc+step_vg2/step_ramp_speed) # Wait 3 times the time constant of the lock-in, plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
            measured_value = measured_parameter()
            x = measured_value['x'][0] #SF: COMMENTED OUT 
            y = measured_value['y'][0]#SF: COMMENTED OUT
            xy_complex = np.complex(x,y) #measured_value #
            v_r_calc = np.absolute(xy_complex)
            theta_calc = np.angle(xy_complex)
                    
            #G calculation
            I = v_r_calc/(gain_RT*gain_HEMT*Z_tot)
            G = 1/((vsdac/I)-Z_tot)
            R = 1/G

            #add to lists (to deal with snake)
            Ilist=Ilist+[I]
            Rlist=Rlist+[R]
            Glist=Glist+[G]
            VRlist=VRlist+[v_r_calc]
            PHASElist=PHASElist+[theta_calc]
            
        #Rlist.reverse
        if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Rlist.reverse()
            Glist.reverse()
            VRlist.reverse()
            PHASElist.reverse()
        #temp_fast_axis_list = list(gate2_sweep) #(to deal with snake)
        #temp_fast_axis_list.reverse
        datasaver.add_result(('R',Rlist),
                            ('G',Glist),
                            ('V_r',VRlist),
                            ('Phase',PHASElist),
                            (gate1_sweep.parameter,gate1_value),
                            (gate2_sweep.parameter,fast_axis_unreversible_list))
        gate2_sweep.reverse()
        reversed_sweep= not reversed_sweep 

# Ramp down everything
gate1.dc_slew_rate_V_per_s(ramp_speed)
gate2.dc_slew_rate_V_per_s(ramp_speed)

gate1.dc_constant_V(0)
gate2.dc_constant_V(0)

qdac.ch03.dc_constant_V(0)
qdac.ch04.dc_constant_V(0)
qdac.ch05.dc_constant_V(0)
qdac.ch06.dc_constant_V(0)
qdac.ch07.dc_constant_V(0)

for i in range(8):
        param1 = getattr(zurich.sigouts.sigouts0.enables, f'enables{i}')
        param2 = getattr(zurich.sigouts.sigouts1.enables, f'enables{i}')
        param1.value(0)
        param2.value(0)
