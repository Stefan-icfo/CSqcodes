# The program does 


import numpy as np


from instruments import qdac, station, k2450,Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2450 as k5
from drivers.Qdac_utils import ramp_QDAC_channel
import time
from tqdm import tqdm
import scipy as scp
#import matplotlib.pyplot as pt

#k2400=keithley2400
#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
gate_ramp_slope = 1e-2 # V/s
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
device_name = 'CD11_d7_C1_'
prefix_name = 'Conductance_IV_cs'
upper_bound_lever_arm=0.5#
Temp=0.02 #Triton.T5()
#if Temp<3:
#    print("Temp<3K")
#    Temp=Triton.MC()
postfix = f"{Temp}K"

vsdkT=Temp/11604
#vsd=vsdkT#automatically sets vsd to kT. comment out if wanna do manually
print(f"vsdkT={vsd}V. ABORT NOW IF FUNKY. u got 10 seconds")
#time.sleep(10)

#####################
#####################
start_vg = -0.57#
stop_vg = -0.54  #
step_num = 401  #
step_num = round((stop_vg-start_vg)/vsd*upper_bound_lever_arm)  #500uV
#####################
step_size=abs(stop_vg-start_vg)/(step_num-1)
step_time=step_size/gate_ramp_slope
print('step time')
print(step_time)
print(f"step num={step_num}")

gates=[qdac.ch01,qdac.ch02,qdac.ch03,qdac.ch04,qdac.ch05]

#--------Definition-------------
source = k2450  # source 
gate1 = qdac.ch01


gate1.label = 'VgDC' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate1])
exp_dict = dict(vsdac = vsd)
exp_name = sample_name(prefix_name,exp_dict,postfix)


# ------------------Create a new Experiment-------------------------

vgdc_sweep = gate1.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here
measured_parameter = k2450.sense.current  # measured parameters will go here

#------------init--------------------
# applied  voltages at the intrument level before attenuation
#source.mode('VOLT')
#source.compliancev = 21
#source.compliancei = 1e-6
#k2450.reset()
k2450.sense.function('current')
k2450.sense.range(1E-5)
k2450.source.function("voltage")





for gate in gates:
    gate.dc_slew_rate_V_per_s(gate_ramp_slope)
    #ramp_QDAC_channel(gate, slew_rate = 1e-2,final_vg = start_vg, ramp_speed = gate_ramp_slope)
    gate.dc_constant_V(start_vg)
print(f"going to sleep for the time it takes to ramp the gate({abs(start_vg-gate.dc_constant_V())/gate_ramp_slope + 30}) plus 30 seconds")
#time.sleep(20)
time.sleep(abs(start_vg-gate.dc_constant_V())/gate_ramp_slope + 30)
print("wake up, gates are")
for gate in gates:
    print(gate.dc_constant_V())
#time.sleep(100)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(vgdc_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
meas.register_custom_parameter('Conductance_IV', 'G_IV', unit='Siemens', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_top', 'I_top', unit='Amps', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_zero', 'I_zero', unit='Amps', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_bottom', 'I_bottom', unit='Amps', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Resistance_IV', 'R_IV', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[vgdc_sweep.parameter])


# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='Gate Sweep', colour = 'green'):
        
        for gate in gates:
            gate.dc_constant_V(vgdc_value)

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
for gate in gates:
            gate.dc_constant_V(0)
time.sleep(abs(stop_vg)/gate_ramp_slope + 30)
print("gaterampsleeping")
print("wake up, gates are")
for gate in gates:
    print(gate.dc_constant_V())

k2450.source.voltage(1e-5)#set to low value but not switching off
#auxgate1(0)
#auxgate2(0)

#k2450.source.voltage(0)
