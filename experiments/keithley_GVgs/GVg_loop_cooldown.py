# GVG loop
# Stefan Forstner with template by Parmeshwar Prasad
# IN CONSTRUCTION


import numpy as np
import os

from instruments import qdac, manual, station, k2400, triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
import time
from tqdm import tqdm
import scipy as scp
#import matplotlib.pyplot as pt


#------User input----------------
gate_ramp_slope = 1e-2 # V/s
ramp_source = 10e-6 # V/ms
ts = 20e-3   # in seconds; settling + measurement time for source (keithley)
tg = 5e-3   # in seconds.settling time for gate (bilt)
vsd = 0.1e-3 # source DC voltage in volt
step_v = abs(vsd) # source steps; for microvolts, one step is ok
offset = 20e-6 #voltage offset, to find zero point
offset_i=-50e-12 #current measurement offset
# device_name = 'tbg_r_3_1'
# prefix_name = 'Conductance_'
# postfix = 'T_10_K_BG'
# exp_name = 'Test 50 K'
timeout=1000 #in seconds

device_name = 'CD06_D3_B3_CS'
prefix_name = '_testnewcode_' 
postfix = '20mK_100uV_gateCS'

#####################
start_vg = 0  #
stop_vg =  0.01   #
step_num = 11    # 200uV steps  
#####################

#--------Definition-------------
source = k2400  # source 
gate = qdac.ch02  # channel of the back gate #ch02 is connected to Gc and ch01 to all other gates
gate.label = 'VgDC' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)


# ------------------Create a new Experiment-------------------------

gate_sweep=gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here
measured_parameter = k2400.curr  # measured parameters will go here


#------------init--------------------
# applied  voltages at the intrument level before attenuation
source.mode('VOLT')
#source.compliancev = 21
source.compliancei = 1e-6


gate.dc_slew_rate_V_per_s(gate_ramp_slope)
gate.dc_constant_V(start_vg)
print("going to sleep for the time it takes to ramp the gate plus 10 seconds")
#time.sleep(10)
time.sleep(abs(start_vg)/gate_ramp_slope + 10)
print("wake up, gate is")
print(gate.dc_constant_V())


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(gate_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas.register_custom_parameter('Current_top', 'I_top', unit='Amps', basis=[], setpoints=[gate_sweep.parameter])
meas.register_custom_parameter('Current_zero', 'I_zero', unit='Amps', basis=[], setpoints=[gate_sweep.parameter])
meas.register_custom_parameter('Current_bottom', 'I_bottom', unit='Amps', basis=[], setpoints=[gate_sweep.parameter])
meas.register_custom_parameter('Resistance_IV', 'R_IV', unit='Ohm', basis=[], setpoints=[gate_sweep.parameter])
meas.register_custom_parameter('Conductance_IV', 'G_IV', unit='Siemens', basis=[], setpoints=[gate_sweep.parameter])
#meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[vgdc_sweep.parameter])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    timezero=time.time()

    for vgdc_value in tqdm(gate_sweep, leave=False, desc='Gate Sweep', colour = 'green'):
        
        gate_sweep.set(vgdc_value)
        #measure lower vsd value (~kT)
        k2.ramp_k2400(ramp_param=source,final_vg=-vsd+offset, step_size = step_v, ramp_speed=1)
        time.sleep(2*ts+2*tg) # Wait 2 times the time constant for the source and gate to settle
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
        
        
        G_IV=result.slope
        R_IV=1/G_IV

        #pt.plot(x,y,'bo')#just to see individual IV fits for test of code
        #pt.show()
        datasaver.add_result(('Conductance_IV', G_IV),('Resistance_IV', R_IV),('Current_top', I_top),('Current_zero', I_zero),('Current_bottom', I_bottom),(gate_sweep.parameter, vgdc_value))
        # vgdc_sweep.reverse()

# Ramp down everything
gate.dc_constant_V(0.0)
k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=1)


print("going to sleep for the time it takes to ramp the gate")
time.sleep(abs(stop_vg)/gate_ramp_slope + 15)
print("wake up, gate is")
print(gate.dc_constant_V())
