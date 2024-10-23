# The program does freq uency sweep at a fixed gate voltage
# Stefan Forstner with template by Parmeshwar Prasad


import pickle

with open('my_list.pkl', 'rb') as f:
    resultlist_minus950 = pickle.load(f)

import numpy as np
import os

from instruments import bilt, station, keithley2400
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
import time
from tqdm import tqdm
import scipy as scp

import matplotlib.pyplot as pt

k2400=keithley2400
#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
ramp_speed = 1.2e-6 # V/ms
ts = 20e-3   # in seconds; settling + measurement time for source (keithley)
tg = 5e-3   # in seconds.settling time for gate (bilt)
vsd = 2e-3 # source DC voltage in volt
step_v = 100e-6 # source steps; for microvolts, one step is ok
offset = 13e-6 #voltage offset, to find zero point
offset_i=-50e-12 #current measurement offset
# device_name = 'tbg_r_3_1'
# prefix_name = 'Conductance_'
# postfix = 'T_10_K_BG'
# exp_name = 'Test 50 K'

#device_name = 'test'
device_name = 'CD08_G7_F6_CSw5gat-950mV'
prefix_name = 'Conductance_IV__1' 
postfix = '1k_2mV_1'

#####################
start_vg = 0.06  #
stop_vg =  0.08  #
step_num = 201     #10uV

avg_num=10
vglist=[]#list of gate voltages
resultlist_now=[]#change notation, change before every meas
resultlist_before=resultlist_minus950#use this for previously saved data
resultlist=resultlist_now
#####################


step_size=abs(stop_vg-start_vg)/(step_num-1)
step_time=step_size/ramp_speed/1000
print('step time')
print(step_time)

#--------Definition-------------
source = k2400  # source 
gate = bilt.ch01.v  # channel of the back gate
#auxgate1=bilt.ch02.v
#auxgate2=bilt.ch03.v
#auxgate1=bilt.ch02.v
#auxgate2=bilt.ch03.v
gate.label = 'VgDC' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)


# ------------------Create a new Experiment-------------------------

vgdc_sweep = gate.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here
measured_parameter = k2400.curr  # measured parameters will go here

#------------init--------------------
# applied  voltages at the intrument level before attenuation
source.mode('VOLT')
source.compliancev = 21
source.compliancei = 1e-6


bilt.channels.output_mode('ramp')
bilt.channels.ramp_slope(ramp_speed)
gate(start_vg)
#auxgate1(start_vg)
#auxgate2(start_vg)
print('wait time')
print(abs(start_vg-gate())/ramp_speed/1000+30)
time.sleep(abs(start_vg-gate())/ramp_speed/1000+30)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(vgdc_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas.register_custom_parameter('Current_top', 'I_top', unit='Amps', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_zero', 'I_zero', unit='Amps', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_bottom', 'I_bottom', unit='Siemens', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Resistance_top', 'R_top', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Resistance_bottom', 'R_bottom', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Resistance_IV', 'R_IV', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Conductance_IV', 'G_IV', unit='Siemens', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Saved_conductance', 'G_IV_saved', unit='Siemens', basis=[], setpoints=[vgdc_sweep.parameter])
#meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[vgdc_sweep.parameter])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    j=0
    for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='Voltage Sweep', colour = 'green'):
         vgdc_sweep.set(vgdc_value)
         #auxgate1(vgdc_value)
         #auxgate2(vgdc_value)
         time.sleep(2*ts+tg+step_time)
         avglist=[]
         for i in range(1,avg_num):
        #auxgate1(vgdc_value)
        #auxgate2(vgdc_value)
        #measure lower vsd value (~kT)
            k2.ramp_k2400(ramp_param=source,final_vg=-vsd+offset, step_size = step_v, ramp_speed=1)
            time.sleep(2*ts) # Wait 2 times the time constant for the source and gate to settle
            measured_value = measured_parameter()
            I_bottom=measured_value-offset_i
            R_bottom = -vsd/I_bottom

        #measure effective zero value 
            k2.ramp_k2400(ramp_param=source,final_vg=+offset, step_size = step_v, ramp_speed=1)
            time.sleep(2*ts) # Wait 2 times the time constant for the source to settle
            measured_value = measured_parameter()
            I_zero=measured_value-offset_i

         #measure lower vsd value (~ - kT) 
            k2.ramp_k2400(ramp_param=source,final_vg=vsd+offset, step_size = step_v, ramp_speed=1)
            time.sleep(2*ts) # Wait 2 times the time constant for the source to settle
            measured_value = measured_parameter()
            I_top=measured_value-offset_i
            R_top = vsd/I_top

        #linear regression fit
            x=[-vsd,0,vsd]
            y=[I_bottom,I_zero,I_top]
            result=scp.stats.linregress(x,y)
        
            avglist=avglist+[result.slope]
            avgslope=sum(avglist)/len(avglist)
         G_IV=avgslope
         R_IV=1/G_IV
         G_IV_saved=resultlist_before[j]
         j=j+1

        #pt.plot(x,y,'bo')#just to see individual IV fits for test of code
        #pt.show()
         #('Saved conductance', G_IV_saved),
         datasaver.add_result(('Saved_conductance', G_IV_saved),('Conductance_IV', G_IV),('Resistance_IV', R_IV),('Resistance_top', R_top),('Resistance_bottom', R_bottom),('Current_top', I_top),('Current_zero', I_zero),('Current_bottom', I_bottom),(vgdc_sweep.parameter, vgdc_value))
        # vgdc_sweep.reverse()
         vglist=vglist+[vgdc_value]
         #resultlist_allghighZ=resultlist_allghighZ+[G_IV]#ADAPT NAME
         resultlist=resultlist+[G_IV]#ADAPT NAME

#pt.plot(vglist,resultlist_1)

# Ramp down everything

#gate(0)
#auxgate1(0)
#auxgate2(0)
resultlist_backup=resultlist
k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=1)
