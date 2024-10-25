# Sweeps gate voltages while reading out with Zurich LIA 
# Sweeps two gates while keeping others constant
# Stefan Forstner using template by Parameshwar Prasad

# snake is working now. watch out - there might still be issues when setting step_ramp_speed to slow and the step_size too large!

import numpy as np
import os

from instruments import  manual, station, qdac, keithley2400
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name

import drivers.k2400 as k2

import time
from tqdm import tqdm

k2400=keithley2400
#------User input----------------
ramp_source = 1e-6 # V/ms
step_v = 1e-6 # source steps
ramp_speed = 0.01 # V/s for large ramps
step_ramp_speed=0.1 # between steps, V/s
tc = 10e-3   # in seconds. 
vsd = 1000e-6 # source DC voltage in volt
device_name = 'CD11_D7_C1'
prefix_name = 'Charge_stability_k2400_QDev'
postfix = '20mK_constantgates135at0.6+1.1+0_postrampingtest'
offset = -10e-6 #voltage offset of k2400
offset_i=-44e-12



#outer voltage range (slow axis)
#####################
start_vg1 = -1.05   #
stop_vg1 = -0.95     #
step_vg1_num = 400+1     #
step_vg1=np.absolute((start_vg1-stop_vg1)/step_vg1_num)


#inner voltage range (fast axis)
#####################
start_vg2 = -1.05     #
stop_vg2 = -0.95       #
step_vg2_num = 400+1    #
step_vg2=np.absolute((start_vg2-stop_vg2)/step_vg2_num)

#constant gate voltages, labelled by the channels they are connected to; 
gate_V_ch3=+1.1
gate_V_ch1=0.6
gate_V_ch5=0
gate_V_ch6=-0.04

#initialize constant gates, comment out for single-gate device

qdac.ch03.dc_slew_rate_V_per_s(ramp_speed)
qdac.ch03.dc_constant_V(gate_V_ch3)
qdac.ch05.dc_slew_rate_V_per_s(ramp_speed)
qdac.ch05.dc_constant_V(gate_V_ch5)
qdac.ch01.dc_slew_rate_V_per_s(ramp_speed)
qdac.ch01.dc_constant_V(gate_V_ch1)
qdac.ch06.dc_slew_rate_V_per_s(ramp_speed)
qdac.ch06.dc_constant_V(gate_V_ch1)
time.sleep(100)#for CS ramp
# qdac.ch06.dc_slew_rate_V_per_s(ramp_speed)
# qdac.ch06.dc_constant_V(gate_V_ch6)
# qdac.ch07.dc_slew_rate_V_per_s(ramp_speed)
# qdac.ch07.dc_constant_V(gate_V_ch7)

#--------Definitions-------------

#swept contacts
gate1=qdac.ch02
  # swept outer gate voltage
gate2=qdac.ch04 #swept inner gate voltage
source = k2400 # source 


gate1.label = 'gate2' # Change the label of the gate1 chanel
gate2.label = 'gate4' # Change the label of the gate2 chaneel
instr_dict = dict(gate1=[gate1])
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)

#----------- defined values------



# ------------------define sweep axes-------------------------

gate1_sweep=gate1.dc_constant_V.sweep(start=start_vg1, stop=stop_vg1, num = step_vg1_num)
gate2_sweep=gate2.dc_constant_V.sweep(start=start_vg2, stop=stop_vg2, num = step_vg2_num)
  # lock-in amplitude measured # SF:CHANGED FROM 'zurich.source.demod_complex'

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT

# applied  voltages at the intrument level before attenuation




#initialize swept contacts

#slow ramp and intial voltage
gate1.dc_slew_rate_V_per_s(ramp_speed)
gate1.dc_constant_V(start_vg1)

gate2.dc_slew_rate_V_per_s(ramp_speed)
gate2.dc_constant_V(start_vg2)
k2.ramp_k2400(ramp_param=source,final_vg=vsd+offset, step_size = step_v, ramp_speed=ramp_source)
#time.sleep(max([abs(start_vg1/ramp_speed),abs(start_vg2/ramp_speed)])+1)  #wait for the time it takes to do both ramps plus one second
measured_parameter = k2400.curr
time.sleep(30)
print("gate channels")
print(qdac.ch01.dc_constant_V())
print(qdac.ch02.dc_constant_V())
print(qdac.ch03.dc_constant_V())
print(qdac.ch04.dc_constant_V())
print(qdac.ch05.dc_constant_V())
print(qdac.ch06.dc_constant_V())
print(qdac.ch07.dc_constant_V())
#set fast ramp speeds
gate1.dc_slew_rate_V_per_s(step_ramp_speed)
gate2.dc_slew_rate_V_per_s(step_ramp_speed)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(gate1_sweep.parameter)  # 
meas.register_parameter(gate2_sweep.parameter)  # 
meas.register_custom_parameter('Conductance', 'G', unit='S', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('Resistance', 'R', unit='Ohm', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])




# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:

    fast_axis_unreversible_list = list(gate2_sweep) #(to deal with snake)
    reversed_sweep=False
    
    for gate1_value in tqdm(gate1_sweep, leave=False, desc='outer gate sweep', colour = 'green'): #slow axis loop (gate)
        gate1_sweep.set(gate1_value)
        time.sleep(1*tc+step_vg1/step_ramp_speed) # Wait 3 times the time contanst of the lock-in plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
        #init lists (to deal with snake)
        Ilist=[]
        Rlist=[]
        Glist=[]
 
        # temp_fast_axis_list=[] #I think this line is unnecessary
        #following is necessary only if no snake #COMMENT OUT ONCE SNAKE IS WORKING
        #gate2.dc_slew_rate_V_per_s(ramp_speed)
        #gate2.dc_constant_V(start_vg2)
        #time.sleep(abs((start_vg2-stop_vg2)/ramp_speed)) 
        #gate2.dc_slew_rate_V_per_s(step_ramp_speed)

        for gate2_value in tqdm(gate2_sweep, leave=False, desc='inner gate sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
            gate2_sweep.set(gate2_value)
            time.sleep(1.1*tc+step_vg2/step_ramp_speed) # Wait 3 times the time constant of the lock-in, plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
            measured_value = measured_parameter()-offset_i
            I=measured_value
            R = vsd/measured_value
            G=1/R

            #add to lists (to deal with snake)
            Ilist=Ilist+[I]
            Rlist=Rlist+[R]
            Glist=Glist+[G]
            
            
        #Rlist.reverse
        if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Rlist.reverse()
            Glist.reverse()
            Ilist.reverse()
          
        #temp_fast_axis_list = list(gate2_sweep) #(to deal with snake)
        #temp_fast_axis_list.reverse
        datasaver.add_result(('Conductance',Glist),('Resistance', Rlist),('Current',Ilist),
                            (gate1_sweep.parameter,gate1_value),
                            (gate2_sweep.parameter,fast_axis_unreversible_list))
        gate2_sweep.reverse()
        reversed_sweep= not reversed_sweep 

# Ramp down everything
gate1.dc_slew_rate_V_per_s(ramp_speed)
gate2.dc_slew_rate_V_per_s(ramp_speed)


k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=ramp_source)
