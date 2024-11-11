# Sweeps gate and source voltages while reading out with keithly
# Stefan Forstner



import numpy as np


from instruments import   station, qdac, keithley2400
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
#from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
#from utils.d2v import d2v
#from utils.v2d import v2d
#from utils.rms2pk import rms2pk
import drivers.k2400 as k2
import time
from tqdm import tqdm

k2400 = keithley2400
#------User input----------------

gate_ramp_slope = 1e-2 # V/s
source_ramp_speed=100e-6 #between steps, V/s
tc = 0.02   # in seconds. 
step_source = 0.1e-3
#vsd_dB = 60 # attenuation at the source in dB
#vsdac = 100e-6 # source AC voltage in volt
device_name = 'CD11_D7_C1_all5g'
#device_name =  'CD05_G6_E3_'# 
prefix_name = '_k2400_'

#prefix_name = 'test'
postfix = '30mK'#'1K5gtrycrossbandgapcs0andg1at1Vconstant'
offset = 30e-6 #voltage offset of k2400
offset_i=-15e-12
# exp_name = 'Test 50 K'

#mix_down_f = 1.25e6 #

#gate voltage range (slow axis)
#####################
start_vg = -2.5 #
stop_vg = 0  #
step_vg_num = 2500 #0.5mV
step_vg=np.absolute((start_vg-stop_vg)/step_vg_num)


#source voltage range (fast axis)
####################
start_vs = -10e-3     #
stop_vs = 10e-3       #
step_vs_num = 101 #  #1mV
step_vs=np.absolute((start_vs-stop_vs)/step_vs_num)

#--------Definitions-------------

#swept contacts0
gates=[qdac.ch01,qdac.ch02,qdac.ch03,qdac.ch04,qdac.ch05]

source = k2400 #swept source voltage




gate1 = qdac.ch01
for gate in gates:
    gate.dc_slew_rate_V_per_s(gate_ramp_slope)
    #ramp_QDAC_channel(gate, slew_rate = 1e-2,final_vg = start_vg, ramp_speed = gate_ramp_slope)
    gate.dc_constant_V(start_vg)
print('wait time')
#time.sleep(10)
print(f"going to sleep for the time it takes to ramp the gate({abs(start_vg-gate.dc_constant_V())/gate_ramp_slope + 30}) plus 30 seconds")
#time.sleep(20)
time.sleep(abs(start_vg-gate.dc_constant_V())/gate_ramp_slope + 30)
print("wake up, gates are")
for gate in gates:
    print(gate.dc_constant_V())




#freq = zurich.oscs.oscs1.freq
gate.label = '5gate voltage' # Change the label of the gate chanel
source.label = 'source_voltage' # Change the label of the source chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = start_vs)
exp_name = sample_name(prefix_name,exp_dict,postfix)
#----------- defined values------



# ------------------define sweep axes-------------------------

vgdc_sweep = gate1.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_vg_num)
source_sweep=source.volt.sweep(start=start_vs+offset, stop=stop_vs+offset, num = step_vs_num)
measured_parameter = k2400.curr   # keithley 2400 current

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT

# applied  voltages at the intrument level before atten
# uation
#vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #SF:this sets up the zurich LIA



#initialize swept contacts

#slow ramp and intial voltage
#gate.dc_slew_rate_V_per_s(gate_ramp_speed)
#gate.dc_constant_V(start_vg)

#source.dc_slew_rate_V_per_s(ramp_speed)
#source.dc_constant_V(start_vs)
k2.ramp_k2400(ramp_param=source,final_vg=start_vs+offset, step_size = step_vs, ramp_speed=source_ramp_speed)

#print("sleep while ramping gate")
#time.sleep(max([abs(start_vg/gate_ramp_speed),abs(start_vs/source_ramp_speed)])+1)  #wait for the time it takes to do both ramps plus one second

#set fast ramp speeds
#gate.dc_slew_rate_V_per_s(gate_ramp_speed)
#source.dc_slew_rate_V_per_s(step_ramp_speed)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(source_sweep.parameter)  # 
meas.register_parameter(vgdc_sweep.parameter)  # 
meas.register_custom_parameter('GIV', 'G_IV', unit='S', basis=[], setpoints=[source_sweep.parameter,vgdc_sweep.parameter])
meas.register_custom_parameter('GIVzero', 'G_IV_zero', unit='S', basis=[], setpoints=[source_sweep.parameter,vgdc_sweep.parameter])
meas.register_custom_parameter('I', 'current', unit='I', basis=[], setpoints=[source_sweep.parameter,vgdc_sweep.parameter])
meas.register_custom_parameter('R', 'R', unit='Ohm', basis=[], setpoints=[source_sweep.parameter,vgdc_sweep.parameter])
#meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[gate_sweep.parameter,source_sweep.parameter])
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[source_sweep.parameter,vgdc_sweep.parameter])
meas.register_custom_parameter('GIV5', 'G_IV5', unit='S', basis=[], setpoints=[source_sweep.parameter,vgdc_sweep.parameter])



# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:
    #save constant gates
    #datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V',qdac.ch01.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch02_dc_constant_V',qdac.ch02.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V',qdac.ch03.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch04_dc_constant_V',qdac.ch04.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch05_dc_constant_V',qdac.ch05.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V',qdac.ch06.dc_constant_V())

    fast_axis_unreversible_list = list(source_sweep) #(to deal with snake)
    reversed_sweep=False

    for gate_value in tqdm(vgdc_sweep, leave=False, desc='Gate Sweep', colour = 'green'): #slow axis loop (gate)
        #print('temperature')
        #Triton.MC()
        for gate in gates:
            gate.dc_constant_V(gate_value)

        time.sleep(tc+step_vg/gate_ramp_slope) # Wait 3 times the time contanst of the lock-in plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
        #init lists (to deal with snake)
        Ilist=[]
        Rlist=[]
        Glist=[]
       
        #VRlist=[]
        #PHASElist=[]
        # temp_fast_axis_list=[] #I think this line is unnecessary
        #following is necessary only if no snake #COMMENT OUT ONCE SNAKE IS WORKING
        #source.dc_slew_rate_V_per_s(ramp_speed)
        #source.dc_constant_V(start_vs)
        #time.sleep(abs((start_vs-stop_vs)/ramp_speed)) 
        #source.dc_slew_rate_V_per_s(step_ramp_speed)

        for source_value in tqdm(source_sweep, leave=False, desc='Source Sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
            source_sweep.set(source_value)
            #time.sleep(tc+step_vs/(source_ramp_speed/1e3)) # Wait keithley, plus the time it takes for the voltage to settle, remember source ramp is in V per ms - doesn't quite work! #SF FIX SLEEP TIMES!
            time.sleep(tc)
            measured_value = measured_parameter()-offset_i-1e-15 #1e-15 is to avoid division by zero
            #x = measured_value['x'][0] #SF: COMMENTED OUT 
            #y = measured_value['y'][0]#SF: COMMENTED OUT
            #xy_complex = measured_value #complex(x,y)
            #v_r_calc = np.absolute(xy_complex)
            #theta_calc = np.angle(xy_complex)
                    
            #G calculation
            R = source_value/measured_value
            G=1/R

            #add to lists (to deal with snake)
            Ilist=Ilist+[measured_value]
            Rlist=Rlist+[R]
            Glist=Glist+[G]
            #VRlist=VRlist+[v_r_calc]
            #PHASElist=PHASElist+[theta_calc]
            
        #Rlist.reverse
        
        #temp_fast_axis_list.reverse()
        if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Ilist.reverse()
            Rlist.reverse()
            Glist.reverse()
            #VRlist.reverse()
            #PHASElist.reverse()
        #InextList=[]+[0]
        GIVlist=[]
        GIVlist5=[]
        GIVlistzero=[]
        n=0
        for I in Ilist:
            if n>0:
                GIVlist=GIVlist+[(Ilist[n]-Ilist[n-1])/step_vs]
                GIVlistzero=GIVlistzero+[max(0,(Ilist[n]-Ilist[n-1]))/step_vs]
            n=n+1
        n=0
        for I in Ilist:
            if n>4:
                GIVlist5=GIVlist5+[(Ilist[n]-Ilist[n-5])/step_vs]
            n=n+1
        #del Ilist[-1:]
        #del Rlist[-1:]
        #del Glist[-1:]
        #GIVlist=GIVlist+[1e-8]
        GIVlist=GIVlist+[GIVlist[-1]]
        GIVlist5=GIVlist5+[GIVlist[-1]]+[GIVlist[-1]]+[GIVlist[-1]]+[GIVlist[-1]]+[GIVlist[-1]]
        GIVlistzero=GIVlistzero+[GIVlistzero[-1]]
        datasaver.add_result(('I',Ilist),
                            ('R',Rlist),
                            ('G',Glist),
                            ('GIV',GIVlist),
                            ('GIV5',GIVlist5),
                            ('GIVzero',GIVlistzero),
                            (vgdc_sweep.parameter,gate_value),
                            (source_sweep.parameter,fast_axis_unreversible_list))
        source_sweep.reverse() 
        reversed_sweep= not reversed_sweep

# Ramp down everything
#gate.dc_slew_rate_V_per_s(ramp_speed)
#source.dc_slew_rate_V_per_s(ramp_speed)

#gate(0)
#auxgate1(0)
#auxgate2(0)
k2.ramp_k2400(source,final_vg=0, step_size = step_source , ramp_speed=source_ramp_speed)

#source.dc_constant_V(0)

#qdac.ch03.dc_constant_V(0)
#qdac.ch04.dc_constant_V(0)
#qdac.ch05.dc_constant_V(0)
#qdac.ch06.dc_constant_V(0)
#qdac.ch07.dc_constant_V(0)


print("going to sleep for the time it takes to ramp the gate plus 10 seconds")
#time.sleep(10)

#time.sleep(abs(stop_vg)/ramp_speed/1000 + 10)
time.sleep(10)
print("wake up, gate and source are")
print(gate.dc_constant_V())
print(k2400.volt())