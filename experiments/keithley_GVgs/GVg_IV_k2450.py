# The program does freq uency sweep at a fixed gate voltage
# Stefan Forstner with template by Parmeshwar Prasad


import numpy as np
import os

from instruments import bilt, station, k2450
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
#import drivers.k2400 as k2
import time
from tqdm import tqdm
import scipy as scp
#import matplotlib.pyplot as pt

#k2400=keithley2400
#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
ramp_speed = 2e-6 # V/ms
ts = 20e-3   # in seconds; settling + measurement time for source (keithley)
tg = 5e-3   # in seconds.settling time for gate (bilt)
vsd = 100e-6 # source DC voltage in volt
step_v = vsd # source steps; for microvolts, one step is ok
offset = 0 #voltage offset, to find zero point
offset_i=+9.5e-12 #current measurement offset
# device_name = 'tbg_r_3_1'
# prefix_name = 'Conductance_'
# postfix = 'T_10_K_BG'
# exp_name = 'Test 50 K'

#device_name = 'test'
device_name = 'CD08_G7_F6_5gates_csat0_200mK'
prefix_name = 'Conductance_IV_k2450_' 
postfix = '1k_200uV_'

#####################
start_vg = -0.7 #
stop_vg =  0  #
step_num = 3501   #0.2mV
#####################
step_size=abs(stop_vg-start_vg)/(step_num-1)
step_time=step_size/ramp_speed/1000
print('step time')
print(step_time)

#--------Definition-------------
source = k2450  # source 
gate = bilt.ch04.v  # channel of the back gate
auxgate1=bilt.ch02.v
auxgate2=bilt.ch03.v
gate.label = 'VgDC' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)


# ------------------Create a new Experiment-------------------------

vgdc_sweep = gate.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here
measured_parameter = k2450.sense.current  # measured parameters will go here

#------------init--------------------
# applied  voltages at the intrument level before attenuation
#source.mode('VOLT')
#source.compliancev = 21
#source.compliancei = 1e-6
#k2450.reset()
k2450.sense.function('current')
k2450.sense.range(1E-8)
k2450.source.function("voltage")





bilt.channels.output_mode('ramp')
bilt.channels.ramp_slope(ramp_speed)
gate(start_vg)
auxgate1(start_vg)
auxgate2(start_vg)
print('wait time')
print(abs(start_vg-gate())/ramp_speed/1000+10)
time.sleep(abs(start_vg-gate())/ramp_speed/1000+10)
#time.sleep(100)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(vgdc_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
meas.register_custom_parameter('Conductance_IV', 'G_IV', unit='Siemens', basis=[], setpoints=[vgdc_sweep.parameter])
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas.register_custom_parameter('Current_top', 'I_top', unit='Amps', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_zero', 'I_zero', unit='Amps', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_bottom', 'I_bottom', unit='Siemens', basis=[], setpoints=[vgdc_sweep.parameter])
#meas.register_custom_parameter('Resistance_top', 'R_top', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
#meas.register_custom_parameter('Resistance_bottom', 'R_bottom', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Resistance_IV', 'R_IV', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[vgdc_sweep.parameter])
#meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[vgdc_sweep.parameter])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='Gate Sweep', colour = 'green'):
        
        vgdc_sweep.set(vgdc_value)
        auxgate1(vgdc_value)
        auxgate2(vgdc_value)
        #measure lower vsd value (~kT)
        #k2.ramp_k2400(ramp_param=source,final_vg=-vsd+offset, step_size = step_v, ramp_speed=1)
        k2450.source.voltage(-vsd+offset)
        time.sleep(2*ts+tg+step_time) # Wait 2 times the time constant for the source and gate to settle
        measured_value = measured_parameter()
        I_bottom=measured_value-offset_i
        R_bottom = -vsd/I_bottom

        #measure effective zero value 
        k2450.source.voltage(+offset)
        time.sleep(2*ts) # Wait 2 times the time constant for the source to settle
        measured_value = measured_parameter()
        I_zero=measured_value-offset_i

         #measure lower vsd value (~ - kT) 
        k2450.source.voltage(+vsd+offset)
        time.sleep(2*ts) # Wait 2 times the time constant for the source to settle
        measured_value = measured_parameter()
        I_top=measured_value-offset_i
        R_top = vsd/I_top

        #linear regression fit
        x=[-vsd,0,vsd]
        y=[I_bottom,I_zero,I_top]
        result=scp.stats.linregress(x,y)
        
        
        G_IV=result.slope
        R_IV=1/G_IV


        #pt.plot(x,y,'bo')#just to see individual IV fits for test of code
        #pt.show()
        datasaver.add_result(('Conductance_IV', G_IV),('Resistance_IV', R_IV),('Current_top', I_top),('Current_zero', I_zero),('Current_bottom', I_bottom),(vgdc_sweep.parameter, vgdc_value))
        # vgdc_sweep.reverse()

# Ramp down everything
#gate(0)
#auxgate1(0)
#auxgate2(0)

#k2450.source.voltage(0)
