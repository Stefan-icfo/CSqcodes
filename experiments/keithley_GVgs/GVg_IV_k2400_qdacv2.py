# The program does gate voltage sweep with 3-point IV curves
# Stefan Forstner 


from instruments import manual, station, keithley2400, qdac
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
import time

from tqdm import tqdm
import scipy as scp



k2400=keithley2400

#------User input----------------
gate_ramp_slope = 1e-2 # V/s
ts = 20e-3   # in seconds; settling + measurement time for source (keithley)
tg = 5e-3   # in seconds.settling time for gate (bilt)
vsd = 16e-6 # source DC voltage in volt
bias=50e-6#7.5e-3
step_v = vsd # source steps; for microvolts, one step is ok
offset = 27e-6 #voltage offset, to find zero point ie. if voltage offset has to be chosen 10uV to reach zero current then wrtie here +10uV
offset_i=-47e-12 #current measurement offset


#device_name = 'test'
device_name = 'CD11_D7_C1'
#device_name = 'test'
prefix_name = '5gateall'
#upper_bound_lever_arm=0.5#
#Temp=Triton.T5()
#if Temp<3:
#    print("Temp<3K")
#    Temp=Triton.MC()
# Temp=20e-3
# postfix = f"{Temp}K"
postfix="30mk"
#vsdkT=Temp/11604
#vsd=vsdkT#automatically sets vsd to kT. comment out if wanna do manually
#print(f"vsdkT={vsd}V. ABORT NOW IF FUNKY. u got 10 seconds")
#time.sleep(10)

upper_bound_lever_arm=0.5
#####################
start_vg = 0
stop_vg =2#
step_num =2000*5#0.2mV
#step_num = round((stop_vg-start_vg)/vsd*upper_bound_lever_arm)+1  #500uV
#####################
step_size=abs(stop_vg-start_vg)/(step_num-1)
step_time=step_size/gate_ramp_slope
print('step time')
print(step_time)
print(f"step num={step_num}")

#--------Definition-------------
source = k2400  # source 3
#gate=bilt.ch03.
gate=qdac.ch04
#gate=qdac.ch06
auxgate1=qdac.ch02
auxgate2=qdac.ch03
#auxgate3=qdac.ch04
#auxgate4=qdac.ch05
#qdac.ch01.dc_constant_V(-3)
#qdac.ch03.dc_constant_V(-3)
#qdac.ch04.dc_constant_V(-3)
#qdac.ch05.dc_constant_V(-3)
gate.label = 'VgDC' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(mV = vsd*1000)
exp_name = sample_name(prefix_name,exp_dict,postfix)



# ------------------Create a new Experiment-------------------------

vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here
measured_parameter = k2400.curr  # measured parameters will go here

#------------init--------------------
# applied  voltages at the intrument level before attenuation



gate.dc_constant_V(start_vg)
#auxgate1.dc_constant_V(start_vg)
#auxgate2.dc_constant_V(start_vg)
#auxgate3.dc_constant_V(start_vg)
#auxgate4.dc_constant_V(start_vg)
print(f"going to sleep for the time it takes to ramp the gate({abs(start_vg-gate.dc_constant_V())/gate_ramp_slope}) plus 30 seconds")
#time.sleep(20)
time.sleep(abs(start_vg-gate.dc_constant_V())/gate_ramp_slope + 30)
print("wake up, gate is")
print(gate.dc_constant_V())


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(vgdc_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
meas.register_custom_parameter('Conductance_IV', 'G_IV', unit='Siemens', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_top', 'I_top', unit='Amps', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_zero', 'I_zero', unit='Amps', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Current_bottom', 'I_bottom', unit='Siemens', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Resistance_IV', 'R_IV', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[vgdc_sweep.parameter])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------
k2.ramp_k2400(ramp_param=source,final_vg=bias, step_size = step_v, ramp_speed=1)
with meas.run() as datasaver:
    # for i in range(2):
    for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='Gate Sweep', colour = 'green'):
        
        gate.dc_constant_V(vgdc_value)
        #auxgate1.dc_constant_V(vgdc_value)
        #auxgate2.dc_constant_V(vgdc_value)
        #auxgate3.dc_constant_V(vgdc_value)
        #auxgate4.dc_constant_V(vgdc_value)
        
        #measure lower vsd value (~kT)
        k2.ramp_k2400(ramp_param=source,final_vg=bias-vsd+offset, step_size = step_v, ramp_speed=1)
        time.sleep(2*ts+tg+step_time) # Wait 2 times the time constant for the source and gate to settle
        measured_value = measured_parameter()
        I_bottom=measured_value-offset_i+1e-15
        R_bottom = -vsd/I_bottom

        #measure effective zero value 
        k2.ramp_k2400(ramp_param=source,final_vg=bias+offset, step_size = step_v, ramp_speed=1)
        time.sleep(2*ts) # Wait 2 times the time constant for the source to settle
        measured_value = measured_parameter()
        I_zero=measured_value-offset_i+1e-15

         #measure upper vsd value (~ + kT) 
        k2.ramp_k2400(ramp_param=source,final_vg=bias+vsd+offset, step_size = step_v, ramp_speed=1)
        time.sleep(2*ts) # Wait 2 times the time constant for the source to settle
        measured_value = measured_parameter()
        I_top=measured_value-offset_i+1e-15
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
#gate.dc_constant_V(0)
#print("gaterampsleep")
#time.sleep(abs(stop_vg)/gate_ramp_slope + 30)
#print("wake up, gate is")
#print(gate.dc_constant_V())
#auxgate1(0)
#auxgate2(0)
#Temp=Triton.T5()
#print(f"Temp after meas={Temp}K")

k2.ramp_k2400(source,final_vg=0, step_size = step_v, ramp_speed=1)


