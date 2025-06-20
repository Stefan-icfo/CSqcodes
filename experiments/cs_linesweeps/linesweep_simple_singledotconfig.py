# weeps gates of both sections for charge sensing
# Stefan Forstner



import numpy as np

from qcodes import Parameter, SweepValues
import numpy as np
from instruments import   station, qdac,  Triton, zurich
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk

import time
from tqdm import tqdm
import scipy as scp
#from utils.CS_utils import zurich_phase_voltage_current_conductance_compensate
from experiment_functions.CS_functions import *


#k2450 = Keithley2450
#------User input----------------
slew_rate=1e-2


tc = 0.1   # in seconds.
tg = 5e-3 
tc = 50e-3   # in seconds.
vsd_dB = 45 # attenuation at the source in dB
vsdac = 16e-6 # source AC voltage in volt at device
device_name = 'CD13_E3_C2'
#device_name =  'CD05_G6_E3_'# 
prefix_name = '_zurich_chargesensing_2d'#

postfix = '32mK'

#Temp=Triton.MC()
#postfix = f"{Temp}K"
#vsdkT=Temp/11604
#vsd=vsdkT
x_avg=+1.7e-5  #+1.51e-5@75#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-1.6e-5
mix_down_f = 1.25e6 # RLC frequency
#outer gate voltage range (slow axis, 5gate)
#####################
start_vg5 =0#y#-1.96 gate2
stop_vg5 =2# #-1.94

step_vg5_num =2000*2 #10uV

step_vg5=np.absolute((start_vg5-stop_vg5)/step_vg5_num)

sleeptime=10

#inner gate voltage range (fast axis, CS)
#####################
start_vgi = 1.4
stop_vgi =  1.8 #-0.776
step_vgi_num = 400*5
#step_vgi_num = round((stop_vgi-start_vgi)/vsd*upper_bound_lever_arm)
#print(f"step i num={step_vgi_num}")
step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)

#--------Definitions-------------

#swept contacts
inner_gate=qdac.ch06.dc_constant_V  # swept gate voltage

outer_gates=[qdac.ch01.dc_constant_V,qdac.ch02.dc_constant_V,qdac.ch03.dc_constant_V,qdac.ch04.dc_constant_V,qdac.ch05.dc_constant_V]

# Define a function to set the same value to all outer gates
def set_outer_gates(value):
    for gate in outer_gates:
        gate(value)  # Set each gate to the same value

# Create a QCoDes parameter to handle the sweeping
sweep_param = Parameter('outer_gates', set_cmd=set_outer_gates)


for outer_gate in outer_gates:
    outer_gate(start_vg5)

inner_gate(start_vgi)
print('manual wait time')
#time.sleep(10)

print(sleeptime)
time.sleep(sleeptime)
for outer_gate in outer_gates:
    print(f"outer gate at {outer_gate()}V ")
print("inner gate")
print(inner_gate())





#freq = zurich.oscs.oscs1.freq


exp_dict = dict(mV = vsdac*1000)
exp_name = sample_name(prefix_name,exp_dict,postfix)
#----------- defined values------
#----------- defined values------


# ------------------define sweep axes-------------------------
sweep_values = np.linspace(start_vg5, stop_vg5, step_vg5_num)

# Use SweepValues to perform the sweep
#outer_gates_sweep = SweepValues(sweep_param, sweep_values)
#outer_gates_sweep=outer_gates.sweep(start=start_vg5, stop=stop_vg5, num = step_vg5_num)
inner_gate_sweep=inner_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)
measured_parameter = zurich.demods.demods0.sample

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT
vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
# applied  voltages at the intrument level before attenuation
#vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #SF:this sets up the zurich LIA



# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(outer_gates[0])  # Register one gate as a representative
#meas.add_before_run(outer_gates[0], {}) 
meas.register_parameter(inner_gate_sweep.parameter)  # 
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[outer_gates[0],inner_gate_sweep.parameter])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[outer_gates[0],inner_gate_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[outer_gates[0],inner_gate_sweep.parameter])
#meas.register_custom_parameter('R', 'R', unit='Ohm', basis=[], setpoints=[outer_gates_sweep.parameter,inner_gate_sweep.parameter])
#meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])


# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:

    fast_axis_unreversible_list = list(inner_gate_sweep) #(to deal with snake)
    reversed_sweep=False
   
    for outer_gates_value in tqdm(sweep_values, leave=False, desc='outer Gates Sweep', colour = 'green'): #slow axis loop (gate)
        
        #print('temperature')
        #Triton.MC()
        set_outer_gates(outer_gates_value)
        
        time.sleep(abs(step_vg5/slew_rate)) # Wait  the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
       
        
        Glist,Vlist,Ilist,Phaselist=GVG_simple(gate_sweep=inner_gate_sweep,
                                                measured_parameter=measured_parameter,
                                                step_sleep_time=1.1*tc+step_vgi/slew_rate,
                                                vsdac=vsdac,
                                                x_avg=x_avg,
                                                y_avg=y_avg,
                                                reverse=reversed_sweep)
        datasaver.add_result(
                            ('G', Glist),
                            ('V_r', Vlist),
                            ('Phase', Phaselist),
                            (outer_gates[0],outer_gates_value),
                            (inner_gate_sweep.parameter,fast_axis_unreversible_list))
        
        
        inner_gate_sweep.reverse() 
        reversed_sweep= not reversed_sweep
  



#print("going to sleep for the time it takes to ramp the gate plus 10 seconds")
#time.sleep(10)

#time.sleep(abs(stop_vg)/ramp_speed/1000 + 10)
for outer_gate in outer_gates:
    print(f"outer gate at {outer_gate()}V ")
print("inner gate")
print(inner_gate())
#print("and source is")
#print(k2400.volt())