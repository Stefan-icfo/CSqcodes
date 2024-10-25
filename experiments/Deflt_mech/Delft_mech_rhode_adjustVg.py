# The program does frequency sweep at a fixed gate voltage
# Stefan Forstner with template by Parmeshwar Prasad

import numpy as np
import os

from instruments import bilt, station, keithley2400, zurich, rohde
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
import time
from tqdm import tqdm

k2400=keithley2400
#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
ramp_gate = 1.2e-6 # V/ms
ramp_source = 10e-6 # V/ms
tc = 10e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd = 3*86e-6 # source DC voltage in volt
step_v = 10e-6 # source steps
offset = 27e-6 #voltage offset of k2400
offset_i=-46e-12
#freq = zurich.oscs.oscs0.freq
freq = rohde.frequency
# device_name = 'tbg_r_3_1'

# prefix_name = 'Conductance_'
# postfix = 'T_10_K_BG'
# exp_name = 'Test 50 K'


#device_name = '100ktest2'
device_name = 'CD08_G7_F6_'
prefix_name = 'mechanics_CB-0.77V39MHz' 
postfix = '1K'


#####################
start_f = 37e6 #Hz unit
stop_f =  41e6 #Hz unit
step_num_f = 2000 #2kHz
#####################
start_vg = -0.775#
stop_vg = -0.765  #
step_num_vg = 101  #0.1mV
step_vg=abs(stop_vg-start_vg)/(step_num_vg-1)
#########################
adjustment_period=60*60*10#never s

#--------Definition-------------
source = k2400 # source 
gate=bilt.ch01.v
auxgate1=bilt.ch02.v

auxgate2=bilt.ch04.v
freq.label = 'Frequency (MHz)' # Change the label of the gate chaneel
instr_dict = dict(freq=[freq])
exp_dict = dict(power = rohde.power())
exp_name = sample_name(prefix_name,exp_dict,postfix)


# ------------------Create a new Experiment-------------------------

freq_sweep = freq.sweep(start=start_f, stop=stop_f, num = step_num_f)  # gate parameter and the instrument will go here
vgdc_sweep = gate.sweep(start=start_vg, stop=stop_vg, num = step_num_vg)
measured_parameter = k2400.curr  # measured parameters will go here

#------------init--------------------


bilt.channels.output_mode('ramp')
bilt.channels.ramp_slope(ramp_gate)
gate(start_vg)
auxgate1(start_vg)
auxgate2(start_vg)
print('wait time')
#time.sleep(10)
sleeptime=max(abs(start_vg-gate()),abs(start_vg-auxgate1()),abs(start_vg-auxgate2()))/ramp_gate/1000+20
print(sleeptime)
time.sleep(sleeptime)
print("wake up, gates are")
print(gate())
#print(outer_auxgate1())
print(auxgate1())
print(auxgate2())


# applied  voltages at the intrument level before attenuation

k2.ramp_k2400(ramp_param=source,final_vg=vsd+offset, step_size = step_v, ramp_speed=ramp_source)



# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(freq_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas.register_custom_parameter('Resistance', 'R', unit='Ohm', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('Current', 'I', unit='A', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('Conductance', 'G', unit='S', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('set_vg', 'V', unit='V', basis=[], setpoints=[freq_sweep.parameter])


# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    GVglist=[]

    next_meas_time=time.monotonic()
    
    for f_value in tqdm(freq_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
    #first adjust Vg if not done so in a while
        if time.monotonic()>=next_meas_time:
            gate(start_vg)
            auxgate1(start_vg)
            auxgate2(start_vg)
            sleeptime=max(abs(start_vg-gate()),abs(start_vg-auxgate1()),abs(start_vg-auxgate2()))/ramp_gate/1000+2
            time.sleep(sleeptime)
            IList=[]#for temporary GVg
            for vgdc_value in vgdc_sweep:
                vgdc_sweep.set(vgdc_value)
                auxgate1(vgdc_value)
                auxgate2(vgdc_value)
                time.sleep(step_vg/ramp_gate/1000)
                measured_value = measured_parameter()-offset_i
                I = measured_value
                IList=IList+[I]
            maxid=IList.index(max(IList))
            set_voltage=vgdc_sweep[maxid]
            gate(set_voltage)
            auxgate1(set_voltage)
            auxgate2(set_voltage)
            sleeptime=max(abs(set_voltage-gate()),abs(set_voltage-auxgate1()),abs(set_voltage-auxgate2()))/ramp_gate/1000+2
            time.sleep(sleeptime)
            next_meas_time=next_meas_time+adjustment_period
            GVglist =  GVglist + [IList]
        
        freq_sweep.set(f_value)
        time.sleep(3*tc) # Wait 3 times the time contanst of the k2400 
        measured_value = measured_parameter()-offset_i
        R = vsd/measured_value
        G=1/R
        datasaver.add_result(('Conductance',G),('Resistance', R),('Current', measured_value),('set_vg',set_voltage),(freq_sweep.parameter, f_value))
        #print(gate())
        # vgdc_sweep.reverse()

# Ramp down everything
#gate(0)
k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=ramp_source)
