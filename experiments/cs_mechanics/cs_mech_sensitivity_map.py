# this code sweeps both gate voltages and for each point does a frequency sweep of the zurich, sweeping both gate and source together with a frequency difference given by the RLC
#the signal is then demodulated by the RLC frequency, the purpose is to measdure the sensitivity to mechanics as a fucntion of gate-voltages
#this code outputs two datasets-one 4d dataset showing Voltage and phase@zurich and Current@device as a function of both gate voltages and frequency
#the second dataset is 3d and shows maximum of voltage and current as well as frequency of this maximum as fkt of gate voltages
# Stefan Forstner


import numpy as np
import time
from tqdm import tqdm
import scipy as scp

from instruments import   station, bilt, Triton, zurich
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
#import drivers.k2400 as k2


def Average(lst): 
    return sum(lst) / len(lst) 

k2400 = keithley2400
#------User input----------------
ramp_mode = 'ramp'  #gate ramp mode
tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd_dB = 45 # attenuation at the source in dB
vsdac = 84e-6 # source AC voltage in volt
gate_ramp_speed = 1.2e-6 # V/ms
tg = 5e-3#settling time of gate

device_name = 'CD08_chipG7_devF6'

prefix_name = '_k2400_chargesensing_mechanics_setpoint_map'
#prefix_name = 'test'
postfix = '1K'



mix_down_f = 1.25e6 #

#Temp=Triton.MC()
#postfix = f"{Temp}K"

#frequency sweep
#####################
start_f = 84e6 #Hz unit
stop_f =  89e6 #Hz unit
step_num_f = 50+1 #100kHz

f_rlc = 1.25e6

#outer gate voltage range (slow axis, 5gate)
#####################
start_vgo = -0.876
stop_vgo = -0.872
step_vgo_num =4*5+1#0.2mV
step_vgo=np.absolute((start_vgo-stop_vgo)/step_vgo_num)
################
cutoff=0.2e-6#V, adjust so that no erroneous data is taken if signal is interrupted 
#inner gate voltage range (fast axis, CS)
#####################
start_vgi = -0.754
stop_vgi = -0.748
step_vgi_num = 6*2 #0.5mV
step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)

#--------Definitions-------------

#swept contacts def
outer_gate=bilt.ch01.v  # swept gate voltage
outer_auxgate1=bilt.ch04.v
outer_auxgate2=bilt.ch02.v
inner_gate=bilt.ch03.v


#zurich oscillators assignment
freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq

#construct name
#instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = vsdac)
exp_name = sample_name(prefix_name,exp_dict,postfix)


#gain and Z on drain line
#####################
gain_RT = 200       #
gain_HEMT = 5.64   #
Z_tot = 7521        #
###################

#initialize constant gates

bilt.channels.output_mode('ramp')
bilt.channels.ramp_slope(gate_ramp_speed)
outer_gate(start_vgo)
outer_auxgate1(start_vgo)
outer_auxgate2(start_vgo)
inner_gate(start_vgi)
print('wait time')
#time.sleep(10)
sleeptime=max(abs(start_vgo-outer_gate()),abs(start_vgo-outer_auxgate2()),abs(start_vgo-outer_auxgate2()),abs(start_vgi-inner_gate()))/gate_ramp_speed/1000+10
print(sleeptime)
time.sleep(sleeptime)
print("wake up, gates are")
print(outer_gate())
print(outer_auxgate1())
print(outer_auxgate2())
print(inner_gate())

#gate labels
outer_gate.label = '5g(outer)' # Change the label of the gate chanel
inner_gate.label = 'CS(inner)' # Change the label of the source chaneel


# ------------------define sweep axes and measuered parameter-------------------------
outer_gate_sweep=outer_gate.sweep(start=start_vgo, stop=stop_vgo, num = step_vgo_num)
inner_gate_sweep=inner_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)
freq_sweep = freq_rf.sweep(start=start_f, stop=stop_f, num = step_num_f)
fsweeplist=list(freq_sweep)
measured_parameter = zurich.demods.demods2.sample   # 


# ----------------Create a measurement-------------------------
#measurement1 - 4d data
experiment1 = new_experiment(name=exp_name+"3d", sample_name=device_name)
meas = Measurement(exp=experiment1)
meas.register_parameter(outer_gate_sweep.parameter)  # 
meas.register_parameter(inner_gate_sweep.parameter)
meas.register_custom_parameter('V_rfmax', 'AmplitudeMax', unit='V', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('maxf', 'maxf', unit='Hz', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('I_rfmax', 'currentMax', unit='I', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])

#measurement2 - 3d data
experiment2 = new_experiment(name=exp_name+"4d", sample_name=device_name)
meas2 = Measurement(exp=experiment2)
meas2.register_parameter(outer_gate_sweep.parameter)  # 
meas2.register_parameter(inner_gate_sweep.parameter)
meas2.register_parameter(freq_sweep.parameter)  # 
meas2.register_custom_parameter('V_rf', 'Amplitude', unit='V', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter,freq_sweep.parameter])
meas2.register_custom_parameter('I_rf', 'current', unit='I', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter,freq_sweep.parameter])
meas2.register_custom_parameter('Phase', 'phase', unit='rad', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter,freq_sweep.parameter])


# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function
with meas2.run() as datasaver2:
    with meas.run() as datasaver1:
    

        VRsweeplist=[2*cutoff]#initialize so the program doesnt get stuck in the beginning

        for outer_gate_value in tqdm(outer_gate_sweep, leave=False, desc='outer Gate Sweep', colour = 'green'): #slow axis loop (gate)
            #print('temperature')
            #Triton.MC()
            outer_gate_sweep.set(outer_gate_value)
            outer_auxgate1(outer_gate_value)
            outer_auxgate2(outer_gate_value)
            inner_gate(start_vgi)
            time.sleep(step_vgo/gate_ramp_speed/1000+abs((stop_vgi-start_vgi)/gate_ramp_speed/1000)) # Wait  the time it takes for the voltage to settle 
            #print("wait time")
            #print(step_vgo/gate_ramp_speed/1000+abs((stop_vgi-start_vgi)/gate_ramp_speed/1000))
            
            for inner_gate_value in tqdm(inner_gate_sweep, leave=False, desc='inner gate Sweep', colour = 'blue'): #fast axis loop 
                inner_gate_sweep.set(inner_gate_value)
                
                if Average(VRsweeplist) >= cutoff:#ie if measurement result makes sense
                    #init lists for 3d data
                    Irfsweeplist=[]
                    VRsweeplist=[]
                    PHASEsweeplist=[]
                    #frequency sweep
                    for f_value in (freq_sweep):
                        #adjust both source and gate frequencies
                        freq_rf(f_value-f_rlc)
                        freq_mech(f_value)
                        time.sleep(1.1*tc)#wait for lock-in integration time
                        measured_value = measured_parameter()
                        x = measured_value['x'][0] # 
                        y = measured_value['y'][0]#
                        xy_complex = complex(x,y)
                        V_rf = np.absolute(xy_complex)
                        phase = np.angle(xy_complex)
                        I_rf = V_rf/(gain_RT*gain_HEMT*Z_tot)
                        #write data in lists for 3d data
                        Irfsweeplist=Irfsweeplist+[I_rf]
                        VRsweeplist=VRsweeplist+[V_rf]
                        PHASEsweeplist=PHASEsweeplist+[phase]
                        #save 4d data
                        datasaver2.add_result(('V_rf',V_rf),
                                    ('I_rf',I_rf),
                                    ('Phase',phase),
                                    (outer_gate_sweep.parameter,outer_gate_value),
                                    (freq_sweep.parameter,f_value),
                                    (inner_gate_sweep.parameter,inner_gate_value))

                else:#ie if result of frequency sweep makes no sense, so signal has dropped
                    while Average(VRsweeplist) <= cutoff:#loop to wait until signal recovers
                        print(f"average={Average(VRsweeplist)}<=cutoff")#print waiting/error message
                        time.sleep(600)
                        Irfsweeplist=[]
                        VRsweeplist=[]
                        PHASEsweeplist=[]
                        for f_value in (freq_sweep):
                            freq_rf(f_value-f_rlc)
                            freq_mech(f_value)
                            time.sleep(1.1*tc)
                            measured_value = measured_parameter()
                            x = measured_value['x'][0]
                            y = measured_value['y'][0]
                            xy_complex = complex(x,y)
                            V_rf = np.absolute(xy_complex)
                            phase = np.angle(xy_complex)
                            I_rf = V_rf/(gain_RT*gain_HEMT*Z_tot)
                            Irfsweeplist=Irfsweeplist+[I_rf]
                            VRsweeplist=VRsweeplist+[V_rf]
                            PHASEsweeplist=PHASEsweeplist+[phase]
                
                max_index = VRsweeplist.index(max(VRsweeplist))
           
                datasaver1.add_result(('I_rfmax',max(Irfsweeplist)),
                                ('V_rfmax',max(VRsweeplist)),
                                ('maxf',fsweeplist[max_index]),
                                (outer_gate_sweep.parameter,outer_gate_value),
                                (inner_gate_sweep.parameter,inner_gate_value))  
    

print("wake up, gates are")
print(outer_gate())
print(outer_auxgate1())
print(outer_auxgate2())
print(inner_gate())
