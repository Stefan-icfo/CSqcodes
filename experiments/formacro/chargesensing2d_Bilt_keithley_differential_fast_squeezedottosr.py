# Sweeps gate and source voltages while reading out with keithly
# Sweeps only single gate while keeping other gates constant
# Stefan Forstner



import numpy as np
import os

from instruments import   station, bilt, keithley2400, Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
#from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
#from utils.d2v import d2v
#from utils.v2d import v2d
#from utils.rms2pk import rms2pk
import drivers.k2400 as k2
import time
from tqdm import tqdm
import scipy as scp

k2400 = keithley2400
#k2450 = Keithley2450
#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
gate_ramp_speed = 5e-6 # V/ms
source_ramp_speed=0.2e-3 #between steps, !
tc = 0.02   # in seconds.
tg = 5e-3 
vsd=0.1e-3
vbias=0#+8.5e-3#0#-13e-3
step_source = vsd
step_v = 100e-6 # source steps; for microvolts, one step is ok
offset = 40e-6 #voltage offset, to find zero point
offset_i=0 #current measurement offset
device_name = 'CD08_chipG7_devF6'
#device_name =  'CD05_G6_E3_'# 
prefix_name = '_k2400_chargesensing_postanneal_sqeezedotright_useg45andg1+1Vg23-3V_srfloating_zoom'
#prefix_name = 'test'
postfix = 'no_postfix'
offset = 27e-6 #voltage offset of k2400
offset_i=-15e-12
upper_bound_lever_arm=0.2
# exp_name = 'Test 50 K'

#mix_down_f = 1.25e6 #

Temp=Triton.MC()
postfix = f"{Temp}K"
vsdkT=Temp/11604
vsd=vsdkT
print(f"vsdkT={vsd}V. ABORT NOW IF FUNKY. u got 10 seconds")
time.sleep(10)

#outer gate voltage range (slow axis, 5gate)
#####################
start_vgo = -3
stop_vgo = -2.8
step_vgo_num =101 #2mV
step_vgo=np.absolute((start_vgo-stop_vgo)/step_vgo_num)


#inner gate voltage range (fast axis, CS)
#####################
start_vgi = -0.025
stop_vgi = 0.025
step_vgi_num = 51#1mV
#step_vgi_num = round((stop_vgi-start_vgi)/vsd*upper_bound_lever_arm)
#print(f"step i num={step_vgi_num}")
step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)

#--------Definitions-------------

#swept contacts
outer_gate=bilt.ch02.v  # swept gate voltage
outer_auxgate1=bilt.ch04.v
outer_auxgate2=bilt.ch01.v
inner_gate=bilt.ch03.v
source = k2400 #swept source voltage


#constant gate voltages, labelled by the channels they are connected to; 
#initialize constant gates




bilt.channels.output_mode('ramp')
bilt.channels.ramp_slope(gate_ramp_speed)
outer_gate(start_vgo)
outer_auxgate1(+1)
outer_auxgate2(-3)
inner_gate(start_vgi)
print('wait time')
#time.sleep(10)
sleeptime=max(abs(start_vgo-outer_gate()),abs(1-outer_auxgate1()),abs(-3-outer_auxgate2()),abs(start_vgi-inner_gate()))/gate_ramp_speed/1000+100
print(sleeptime)
time.sleep(sleeptime)
print("wake up, gates are")
print(outer_gate())
print(outer_auxgate1())
print(outer_auxgate2())
print(inner_gate())





#freq = zurich.oscs.oscs1.freq
outer_gate.label = '5g(outer)' # Change the label of the gate chanel
inner_gate.label = 'CS(inner)' # Change the label of the source chaneel
instr_dict = dict(gate=[outer_gate])
exp_dict = dict(mV = vsd*1000)
exp_name = sample_name(prefix_name,exp_dict,postfix)
#----------- defined values------


# ------------------define sweep axes-------------------------

outer_gate_sweep=outer_gate.sweep(start=start_vgo, stop=stop_vgo, num = step_vgo_num)
inner_gate_sweep=inner_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)
measured_parameter = source.curr   # keithley 2400 current

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT

# applied  voltages at the intrument level before attenuation
#vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #SF:this sets up the zurich LIA



#initialize swept contacts

#slow ramp and intial voltage
#gate.dc_slew_rate_V_per_s(gate_ramp_speed)
#gate.dc_constant_V(start_vg)

#source.dc_slew_rate_V_per_s(ramp_speed)
#source.dc_constant_V(start_vs)
#k2.ramp_k2400(ramp_param=source,final_vg=vsd+offset, step_size = step_vs, ramp_speed=source_ramp_speed)

#print("sleep while ramping gate")
#time.sleep(max([abs(start_vg/gate_ramp_speed),abs(start_vs/source_ramp_speed)])+1)  #wait for the time it takes to do both ramps plus one second

#set fast ramp speeds
#gate.dc_slew_rate_V_per_s(gate_ramp_speed)
#source.dc_slew_rate_V_per_s(step_ramp_speed)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(outer_gate_sweep.parameter)  # 
meas.register_parameter(inner_gate_sweep.parameter)  # 
meas.register_custom_parameter('Current_top', 'I_top', unit='Amps', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('Current_zero', 'I_zero', unit='Amps', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('Current_bottom', 'I_bottom', unit='Siemens', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('Resistance_IV', 'R_IV', unit='Ohm', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('Conductance_IV', 'G_IV', unit='Siemens', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])

k2.ramp_k2400(ramp_param=source,final_vg=vbias-offset, step_size = step_v, ramp_speed=source_ramp_speed)
#time.sleep(source_ramp_speed/vbias)

# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:

    fast_axis_unreversible_list = list(inner_gate_sweep) #(to deal with snake)
    reversed_sweep=False

    for outer_gate_value in tqdm(outer_gate_sweep, leave=False, desc='outer Gate Sweep', colour = 'green'): #slow axis loop (gate)
        #print('temperature')
        #Triton.MC()
        outer_gate_sweep.set(outer_gate_value)
        #outer_auxgate1(outer_gate_value)
        #outer_auxgate2(outer_gate_value)
        time.sleep(step_vgo/gate_ramp_speed/1000) # Wait  the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
        #init lists (to deal with snake)
        Itoplist=[]
        Izerolist=[]
        Ibottomlist=[]
        RIVlist=[]
        GIVlist=[]
       
        #VRlist=[]
        #PHASElist=[]
        # temp_fast_axis_list=[] #I think this line is unnecessary
        #following is necessary only if no snake #COMMENT OUT ONCE SNAKE IS WORKING
        #source.dc_slew_rate_V_per_s(ramp_speed)
        #source.dc_constant_V(start_vs)
        #time.sleep(abs((start_vs-stop_vs)/ramp_speed)) 
        #source.dc_slew_rate_V_per_s(step_ramp_speed)
       
        
        for inner_gate_value in tqdm(inner_gate_sweep, leave=False, desc='inner gate Sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
            inner_gate_sweep.set(inner_gate_value)
            k2.ramp_k2400(ramp_param=source,final_vg=vbias-vsd+offset, step_size = step_v, ramp_speed=source_ramp_speed)
           
            time.sleep(2*tc+step_vgi/gate_ramp_speed/1000) 
            measured_value = measured_parameter()
            I_bottom=measured_value-offset_i

             #measure effective zero value 
            
             # Wait 2 times the time constant for the source to settle
            
            k2.ramp_k2400(ramp_param=source,final_vg=vbias+offset, step_size = step_v, ramp_speed=source_ramp_speed)
            time.sleep(2*tc)
            measured_value = measured_parameter()
            I_zero=measured_value-offset_i

            #measure lower vsd value (~ - kT) 
            k2.ramp_k2400(ramp_param=source,final_vg=vbias+vsd+offset, step_size = step_v, ramp_speed=source_ramp_speed)
            time.sleep(2*tc) # Wait 2 times the time constant for the source to settle
            measured_value = measured_parameter()
            I_top=measured_value-offset_i
            

        #linear regression fit
            x=[-vsd,0,vsd]
            y=[I_bottom,I_zero,I_top]
            result=scp.stats.linregress(x,y)
        
        
            G_IV=result.slope
            R_IV=1/G_IV
            #add to lists (to deal with snake)
            Itoplist=Itoplist+[I_top]
            Izerolist=Izerolist+[I_zero]
            Ibottomlist=Ibottomlist+[I_bottom]
            RIVlist=RIVlist+[R_IV]
            GIVlist=GIVlist+[G_IV]
            
           
            #VRlist=VRlist+[v_r_calc]
            #PHASElist=PHASElist+[theta_calc]
            
        #Rlist.reverse
        
        #temp_fast_axis_list.reverse()
        if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Itoplist.reverse()
            Izerolist.reverse()
            Ibottomlist.reverse()
            RIVlist.reverse()
            GIVlist.reverse()
            #VRlist.reverse()
            #PHASElist.reverse()
        datasaver.add_result(('Conductance_IV', GIVlist),
                            ('Resistance_IV', RIVlist),
                            ('Current_top', Itoplist),
                            ('Current_zero', Izerolist),
                            ('Current_bottom', Ibottomlist),
                            (outer_gate_sweep.parameter,outer_gate_value),
                            (inner_gate_sweep.parameter,fast_axis_unreversible_list))
       
        inner_gate_sweep.reverse() 
        reversed_sweep= not reversed_sweep

# Ramp down everything
#gate.dc_slew_rate_V_per_s(ramp_speed)
#source.dc_slew_rate_V_per_s(ramp_speed)

#gate(0)
#auxgate1(0)
k2.ramp_k2400(source,final_vg=0, step_size = step_source , ramp_speed=source_ramp_speed)
time.sleep(vsd/source_ramp_speed)
#source.dc_constant_V(0)

#qdac.ch03.dc_constant_V(0)
#qdac.ch04.dc_constant_V(0)
#qdac.ch05.dc_constant_V(0)r
#qdac.ch06.dc_constant_V(0)
#qdac.ch07.dc_constant_V(0)


#print("going to sleep for the time it takes to ramp the gate plus 10 seconds")
#time.sleep(10)

#time.sleep(abs(stop_vg)/ramp_speed/1000 + 10)
print("wake up, gates are")
print(outer_gate())
print(outer_auxgate1())
print(outer_auxgate2())
print(inner_gate())
print("and source is")
print(k2400.volt())