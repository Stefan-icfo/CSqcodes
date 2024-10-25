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
from utils.CS_utils import breit_wigner_fkt, zurich_phase_voltage_current_conductance_compensate,save_metadata_var, get_var_name
import os
from qcodes import Parameter



code_start_time=time.time_ns()/1e9
#------User input----------------
slew_rate=1e-2
x_avg=+4.38e-6
y_avg=-4.41e-6

tc = 0.1   # in seconds.
tg = 5e-3 
tc = 100e-3   # in seconds.
vsd_dB = 45 # attenuation at the source in dB
vsdac = 15e-6 # source AC voltage in volt
device_name = 'CD11_D7_C1'
#device_name =  'CD05_G6_E3_'# 
prefix_name = 'linescan_w_fit_g1at0.5_test_noslepingtime'#

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
vgo1_init_list=[-1.9389]
vgo2_init_list=[-1.9347]

vgo1 =  -1.939#y#-1.96 gate2
#stop_vgo1 =   -1.93901# #-1.94
vgo2 =  -1.9805#x#-1.9601 gate4
#stop_vgo2 =   -1.9805# #-1.960

stabilisation_time=10#s
cutoff_time=1800#s

vars_to_save=[vgo1,vgo2,cutoff_time,stabilisation_time]
#vars_to_save.append(vgo1_init_list)
#vars_to_save.append(vgo2_init_list)


#step_vgo_num =10*20 #50uV

#step_vgo1=np.absolute((start_vgo1-stop_vgo1)/step_vgo_num)
#step_vgo2=np.absolute((start_vgo2-stop_vgo2)/step_vgo_num)

#vars_to_save=[]

#inner gate voltage range (fast axis, CS)
#####################
start_vgi = -1.716#-0.788

stop_vgi = -1.715#-0.7763
step_vgi_num = 50#20uV
#step_vgi_num = round((stop_vgi-start_vgi)/vsd*upper_bound_lever_arm)
#print(f"step i num={step_vgi_num}")
step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)

initial_guess = [(start_vgi+stop_vgi)/2, 1e-4, 6e-6]#initial guess for peakV, Gamma,height for first GVg

vars_to_save.extend([start_vgi,stop_vgi,step_vgi_num,step_vgi])

#postfix = f"hysteresis_scan_{Temp}K,g1={round(qdac.ch01.dc_constant_V(),2)},g3={round(qdac.ch03.dc_constant_V(),2)},g5={round(qdac.ch05.dc_constant_V(),2)},gcsstart={round(start_vgi,2)}"
postfix = "hysteresis_scan"
#--------Definitions-------------

#swept contacts
inner_gate=qdac.ch06.dc_constant_V  # swept gate voltage

outer_gate1=qdac.ch02.dc_constant_V
outer_gate2=qdac.ch04.dc_constant_V

#initialisating

inner_gate(start_vgi)

sleeptime=abs(start_vgi-inner_gate())/slew_rate
print(f"ramping inner gate from {round(inner_gate(),4)}V to {start_vgi}V, sleeptime={round(sleeptime,3)}s")
time.sleep(sleeptime)
print("starting outer gate initialization ")
print(f"starting_values: x={round(outer_gate2(),4)}, y={round(outer_gate1(),4)}")

for x,y in zip(vgo2_init_list,vgo1_init_list):
    outer_gate2(x)
    outer_gate1(y)
    sleeptime=max(abs(x-outer_gate2()),abs(y-outer_gate1()))/slew_rate
    print(f"setting to x={x}V, y={y}V, waiting {round(sleeptime,3)}s")
    time.sleep(sleeptime)

print(f"waiting for additional {stabilisation_time}s for configuration to stabilize")
time.sleep(stabilisation_time)
outer_gate1(vgo1)
outer_gate2(vgo2)
finalramptime=max(abs(vgo2-outer_gate2()),abs(vgo1-outer_gate1()))/slew_rate
print(f"ramping outer gates to final values x={vgo2}V, y={vgo1}V")#, sleeptime = {round(sleeptime,3)}s")
#time.sleep(sleeptime)
meas_start_time=time.time_ns()/1e9-finalramptime
meas_time_passed=-finalramptime
#inner_gate(start_vgi)
#print('wait time')
#time.sleep(10)
#sleeptime=max(abs(vgo1-outer_gate1()),abs(vgo2-outer_gate2()),abs(start_vgi-inner_gate()))/slew_rate
#print(sleeptime)
#time.sleep(sleeptime)
print(f"initialisation finished, gates are, time since start={time.time_ns()-code_start_time}s")
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

#outer_gate1_sweep=outer_gate1.sweep(start=start_vgo1, stop=stop_vgo1, num = step_vgo_num)
#outer_gate2_sweep=outer_gate2.sweep(start=start_vgo2, stop=stop_vgo2, num = step_vgo_num)
#outer_gate1_list=list(outer_gate1_sweep)
#outer_gate2_list=list(outer_gate2_sweep)
inner_gate_sweep=inner_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)
measured_parameter = zurich.demods.demods0.sample

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT
vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
# applied  voltages at the intrument level before attenuation
#vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #SF:this sets up the zurich LIA

time_param = Parameter('time_param',
                            label='time',
                            unit='s',  # Or the appropriate unit
                            set_cmd=None,  # If there is no instrument to communicate with directly
                            get_cmd=None)  # Define get_cmd if you need to read a value


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(time_param)  # 
meas.register_parameter(inner_gate_sweep.parameter)  # 
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[time_param,inner_gate_sweep.parameter])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[time_param,inner_gate_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[time_param,inner_gate_sweep.parameter])
meas.register_custom_parameter('other_gate', 'other_gate', unit='V', basis=[], setpoints=[time_param,inner_gate_sweep.parameter])
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
    varnames=[]
    #for i in range(len(vars_to_save)):
    #    varnames.append(get_var_name(vars_to_save[i]))
    #save_metadata_var(datasaver.dataset,varnames,vars_to_save)

    y_names=[]
    x_names=[]
    for i in range(len(vgo2_init_list)):
        x_names.append(f"x{i}")
        y_names.append(f"y{i}")
        


    #save_metadata_var(datasaver.dataset,x_names,vgo2_init_list[i])
    #save_metadata_var(datasaver.dataset,y_names,vgo1_init_list[i])



    fast_axis_unreversible_list = list(inner_gate_sweep) #(to deal with snake)
    fast_axis_unreversible_array=np.array(fast_axis_unreversible_list)
    reversed_sweep=False
    i=0
    peakfitlist=[]
    timeslist=[]
    First_run=True
   
    while meas_time_passed<cutoff_time:
        i=i+1#outergatesweepcounter
        #print('temperature')
        #Triton.MC()
        current_time=time.time_ns()/1e9
        meas_time_passed=current_time-meas_start_time
        Glist=[]
        Vlist=[]
        Phaselist=[]
        
        
        
        for inner_gate_value in tqdm(inner_gate_sweep, leave=False, desc='inner gate Sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
            inner_gate_sweep.set(inner_gate_value)
            time.sleep(1.1*tc+step_vgi/slew_rate) # Wait 3 times the time contanst of the lock-in plus gate ramp speed
            measured_value = measured_parameter()
            theta_calc, v_r_calc, I,  G = zurich_phase_voltage_current_conductance_compensate(measured_value, vsdac,x_avg, y_avg)

            Glist=Glist+[G]
            Vlist=Vlist+[v_r_calc]
            Phaselist=Phaselist+[theta_calc]
        #temp_fast_axis_list.reverse()
        if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Glist.reverse()
            Vlist.reverse()
            Phaselist.reverse()
        
        if First_run==False: #ie if it's not the first run and popt has been measured, then redefine the initial guess by using the last fitted values
            initial_guess = [popt[0], popt[1], popt[2]]
        First_run=False
        popt, pcov = scp.optimize.curve_fit(breit_wigner_fkt, fast_axis_unreversible_list, Glist, p0=initial_guess)
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
        

        print(f"meastimepassed={meas_time_passed}s")
        timeslist.append(meas_time_passed)
        datasaver.add_result(('G', Glist),
                            ('V_r', Vlist),
                            ('Phase', Phaselist),
                            (time_param,meas_time_passed),
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
plt.plot(timeslist,peakfitlist)
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


