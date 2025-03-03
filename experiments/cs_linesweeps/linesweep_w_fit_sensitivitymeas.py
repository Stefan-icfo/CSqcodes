# sweeps plunger gates diagonally while doing full GVg of charge sensor for each setpoint
# sweeps between two given points
# Stefan Forstner
#suggested name:linesweep_w_fit


import numpy as np


from instruments import   station, qdac,  Triton, zurich
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
#from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk

import time
from tqdm import tqdm
import scipy as scp
import matplotlib.pyplot as plt
from utils.CS_utils import breit_wigner_fkt, zurich_phase_voltage_current_conductance_compensate
import os




#------User input----------------
slew_rate=1e-2
x_avg=+4.6e-6
y_avg=-6.56e-6

tc = 0.1   # in seconds.
tg = 5e-3 
tc = 100e-3   # in seconds.
vsd_dB = 45 # attenuation at the source in dB
vsdac = 15e-6 # source AC voltage in volt
device_name = 'CD11_D7_C1'
#device_name =  'CD05_G6_E3_'# 
prefix_name = 'linescan_w_fit_g1g50.1g31.1ITC'#

#postfix = '20mK'

# exp_name = 'Test 50 K'

#mix_down_f = 1.25e6 #

Temp=Triton.MC()

#vsdkT=Temp/11604
#vsd=vsdkT

mix_down_f = 1.25e6 # RLC frequency
zurich.oscs.oscs0.freq(mix_down_f)
#outer gate voltage range (slow axis, 5gate)
#####################
start_vgo1 =  -1.9861#y#-1.96 gate2
stop_vgo1 =   -1.9873# #-1.94
start_vgo2 =  -1.9319#x#-1.9601 gate4
stop_vgo2 =   -1.93168# #-1.960

step_vgo_num =10*20 #50uV

step_vgo1=np.absolute((start_vgo1-stop_vgo1)/step_vgo_num)
step_vgo2=np.absolute((start_vgo2-stop_vgo2)/step_vgo_num)

#vars_to_save=[]

#inner gate voltage range (fast axis, CS)
#####################
start_vgi = -1.718#-0.788

stop_vgi = -1.716#-0.7763
step_vgi_num = 2*50#20uV
#step_vgi_num = round((stop_vgi-start_vgi)/vsd*upper_bound_lever_arm)
#print(f"step i num={step_vgi_num}")
step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)

initial_guess = [(start_vgi+stop_vgi)/2, 1e-4, 6e-6]#initial guess for peakV, Gamma,height for first GVg

postfix = f"mod_freq={zurich.demods.demods1.freq()},T={Temp}K,g1={round(qdac.ch01.dc_constant_V(),2)},g3={round(qdac.ch03.dc_constant_V(),2)},g5={round(qdac.ch05.dc_constant_V(),2)},gcsstart={round(start_vgi,2)}"

#--------Definitions-------------

#swept contacts
inner_gate=qdac.ch06.dc_constant_V  # swept gate voltage

outer_gate1=qdac.ch02.dc_constant_V
outer_gate2=qdac.ch04.dc_constant_V

#constant gate voltages, labelled by the channels they are connected to; 
#gate_V_ch3=+1
#gate_V_ch1=-3
#gate_V_ch5=-3

#initialize constant gates, comment out for single-gate device

#qdac.ch03.dc_slew_rate_V_per_s(slew_rate)
#qdac.ch03.dc_constant_V(gate_V_ch3)
#qdac.ch05.dc_slew_rate_V_per_s(slew_rate)
#qdac.ch05.dc_constant_V(gate_V_ch5)
#qdac.ch01.dc_slew_rate_V_per_s(slew_rate)
#qdac.ch01.dc_constant_V(gate_V_ch1)



outer_gate1(start_vgo1)
outer_gate2(start_vgo2)

inner_gate(start_vgi)
print('wait time')
#time.sleep(10)
sleeptime=max(abs(start_vgo1-outer_gate1()),abs(start_vgo2-outer_gate2()),abs(start_vgi-inner_gate()))/slew_rate+10
print(sleeptime)
time.sleep(sleeptime)
print("wake up, gates are")
print(outer_gate1())
#print(outer_auxgate1())
print(outer_gate2())
print(inner_gate())





#freq = zurich.oscs.oscs1.freq
outer_gate1.label = 'gate2' # Change the label of the gate chanel
inner_gate.label = 'CS(inner)' # Change the label of the source chaneel
instr_dict = dict(gate=[outer_gate1])
exp_dict = dict(mV = vsdac*1000)
exp_name = sample_name(prefix_name,exp_dict,postfix)
#----------- defined values------
#----------- defined values------
#####################
#gain_RT = 200       #
#gain_HEMT = 5.64   #
#Z_tot = 7521        #
###################

# ------------------define sweep axes-------------------------

outer_gate1_sweep=outer_gate1.sweep(start=start_vgo1, stop=stop_vgo1, num = step_vgo_num)
outer_gate2_sweep=outer_gate2.sweep(start=start_vgo2, stop=stop_vgo2, num = step_vgo_num)
outer_gate1_list=list(outer_gate1_sweep)
outer_gate2_list=list(outer_gate2_sweep)
inner_gate_sweep=inner_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)
measured_parameter0 = zurich.demods.demods0.sample
measured_parameter1 = zurich.demods.demods1.sample

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT
vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
# applied  voltages at the intrument level before attenuation
#vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #SF:this sets up the zurich LIA



# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(outer_gate1_sweep.parameter)  # 
meas.register_parameter(inner_gate_sweep.parameter)  # 
meas.register_custom_parameter('G0', 'G0', unit='S', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('V_r0', 'Amplitude0', unit='V', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('Phase0', 'Phase0', unit='rad', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('other_gate', 'other_gate', unit='V', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])

meas.register_custom_parameter('V_r1', 'Amplitude1', unit='V', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('Phase1', 'Phase1', unit='rad', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
#meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
#ManualParameter(name='fittedpeaklist', initial_value=[])
#meas.register_parameter(fittedpeaklist)

# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:
    #save barrier gates
    datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V',qdac.ch01.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch02_dc_constant_V',qdac.ch02.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V',qdac.ch03.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch04_dc_constant_V',qdac.ch04.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch05_dc_constant_V',qdac.ch05.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V',qdac.ch06.dc_constant_V())

    fast_axis_unreversible_list = list(inner_gate_sweep) #(to deal with snake)
    fast_axis_unreversible_array=np.array(fast_axis_unreversible_list)
    reversed_sweep=False
    i=0
    peakfitlist=[]
    First_run=True
    for outer_gate_value in tqdm(outer_gate1_sweep, leave=False, desc='outer Gate Sweep', colour = 'green'): #slow axis loop (gate)
        i=i+1#outergatesweepcounter
        #print('temperature')
        #Triton.MC()
        outer_gate1_sweep.set(outer_gate_value)
        outer_gate2_sweep.set(outer_gate2_sweep[i-1])
        time.sleep(max(abs(step_vgo1/slew_rate),abs(step_vgo2/slew_rate))) # Wait  the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
        Glist0=[]
        Vlist0=[]
        Phaselist0=[]
        Vlist1=[]
        Phaselist1=[]
        
        
        for inner_gate_value in tqdm(inner_gate_sweep, leave=False, desc='inner gate Sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
            inner_gate_sweep.set(inner_gate_value)
            time.sleep(1.1*tc+step_vgi/slew_rate) # Wait 3 times the time contanst of the lock-in plus gate ramp speed
            measured_value0 = measured_parameter0()
            measured_value1 = measured_parameter1()
            theta_calc0, v_r_calc0, I0,  G0 = zurich_phase_voltage_current_conductance_compensate(measured_value0, vsdac,x_avg, y_avg)
            theta_calc1, v_r_calc1, I1,  G1 = zurich_phase_voltage_current_conductance_compensate(measured_value1, vsdac,x_avg, y_avg)

            Glist0=Glist0+[G0]
            Vlist0=Vlist0+[v_r_calc0]
            Phaselist0=Phaselist0+[theta_calc0]
            Vlist1=Vlist1+[v_r_calc1]
            Phaselist1=Phaselist1+[theta_calc1]
        #temp_fast_axis_list.reverse()
        if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Glist0.reverse()
            Vlist0.reverse()
            Phaselist0.reverse()
            Vlist1.reverse()
            Phaselist1.reverse()
        
        if First_run==False: #ie if it's not the first run and popt has been measured, then redefine the initial guess by using the last fitted values
            initial_guess = [popt[0], popt[1], popt[2]]
        First_run=False
        popt, pcov = scp.optimize.curve_fit(breit_wigner_fkt, fast_axis_unreversible_list, Glist0, p0=initial_guess)
        peak_fit, hgamma_fit, B_fit=popt
        #peak_fit, hgamma_fit, B_fit, pcov = find_peak_fit(fast_axis_unreversible_array,np.array(Glist),initial_guess)
        B_0=initial_guess[2]
        hgamma_0=initial_guess[1]
        peak_0=initial_guess[0]
        peakfitlist.append(peak_fit)

        #plt.figure(1)
        #plt.plot(fast_axis_unreversible_list, Glist)
        #y_fit = breit_wigner_fkt(fast_axis_unreversible_array, peak_fit, hgamma_fit, B_fit)
        #y_fit_initial = breit_wigner_fkt(fast_axis_unreversible_array, peak_0, hgamma_0, B_0)
        #plt.plot(fast_axis_unreversible_list, y_fit)
        #plt.plot(fast_axis_unreversible_list, y_fit_initial)
        #plt.title('fit_initial')
        #plt.show()
        


        datasaver.add_result(('G0', Glist0),
                            ('V_r0', Vlist0),
                            ('Phase0', Phaselist0),
                            ('V_r1', Vlist1),
                            ('Phase1', Phaselist1),
                            ('other_gate', [outer_gate2_sweep[i-1]]*len(fast_axis_unreversible_list)),
                            (outer_gate1_sweep.parameter,outer_gate_value),
                            (inner_gate_sweep.parameter,fast_axis_unreversible_list))
        
      
        inner_gate_sweep.reverse() 
        reversed_sweep= not reversed_sweep

    #fittedpeaklist.set(peakfitlist)
    #datasaver.add_result(fittedpeaklist,fittedpeaklist.get())

# Ramp down everything
#gate.dc_slew_rate_V_per_s(ramp_speed)
#source.dc_slew_rate_V_per_s(ramp_speed)

#gate(0)
#auxgate1(0)


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
print(outer_gate1())
print(outer_gate2())
print(inner_gate())
#print("and source is")
#print(k2400.volt())
plt.plot(outer_gate2_list,peakfitlist)
plt.show()

foldername='C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_'
if not os.path.exists(foldername):
    os.makedirs(foldername) 
run_id = datasaver.run_id

filename=f'meas{run_id}_peakpos_V.npy'
path = os.path.join(foldername, filename)
np.save(path, np.array(peakfitlist))

filename=f'meas{run_id}_gate1_V.npy'
path = os.path.join(foldername, filename)
np.save(path, np.array(outer_gate1_list))

filename=f'meas{run_id}_gate2_V.npy'
path = os.path.join(foldername, filename)
np.save(path, np.array(outer_gate2_list))


