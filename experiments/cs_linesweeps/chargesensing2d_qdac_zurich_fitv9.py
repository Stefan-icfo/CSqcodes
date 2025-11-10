# sweeps plunger gates diagonally while doing full GVg of charge sensor for each setpoint
# calculates sweep axis based on triple points and delta
# Stefan Forstner



import numpy as np


from instruments import   station, qdac,  Triton, zurich,exp
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
from utils.CS_utils import breit_wigner_fkt, idt_perpendicular_angle, make_detuning_axis_noncenter, make_detuning_axis_noncenterM
import os
from qcodes import Parameter


from instruments import exp

#------User input----------------
slew_rate=1e-2
x_avg=+4.38e-6
y_avg=-4.41e-6

tc = 0.1   # in seconds.
tg = 5e-3 
tc = 30e-3   # in seconds.
vsd_dB = 42.3 # attenuation at the source in dB
vsdac = 10.9e-6# source AC voltage in volt
device_name = exp.device_name#'CD12'
#device_name =  'CD05_G6_E3_'# 
prefix_name ='test'#






# exp_name = 'Test 50 K'

#mix_down_f = 1.25e6 #

#Temp=Triton.MC()
#postfix = f"{Temp}K"
#vsdkT=Temp/11604
#vsd=vsdkT

#mix_down_f = 1.25e6 # RLC frequency
#zurich.oscs.oscs0.freq(mix_down_f)
#outer gate voltage range (slow axis, 5gate)
#####################


#idt_point1_x=-1.5234
#idt_point1_y=-2.3726
#idt_point2_x=-1.5161
#idt_point2_y=-2.3674
#delta=2e-3

idt_point2_x=0.94113
idt_point2_y=1.0352
idt_point1_x=0.93758
idt_point1_y=1.0496
delta=1.5e-3
step_vgo_num = 30 +1

#step_vgo_num =90+1 #
xi=0#move along ict (take traces not through centerbut closer to  triple pt)
epsilon_0 =1e-3#-900e-6#move prependicular to ict (compensate for drift)
start_vgo2,start_vgo1,stop_vgo2,stop_vgo1=make_detuning_axis_noncenterM(idt_point1_x,idt_point1_y,idt_point2_x,idt_point2_y,delta,xi,epsilon_0) 

print(f"start_vgo2={start_vgo2},start_vgo2={start_vgo2},stop_vgo2={stop_vgo2},stop_vgo1={stop_vgo1}")

step_vgo1=np.absolute((start_vgo1-stop_vgo1)/step_vgo_num)
step_vgo2=np.absolute((start_vgo2-stop_vgo2)/step_vgo_num)

vars_to_save=[slew_rate,tc,x_avg,y_avg,idt_point1_x,idt_point1_y,idt_point2_x,idt_point2_y,delta,step_vgo_num]



postfix = f"re={zurich.output1_amp1():4g},epsilon_0={epsilon_0},g1={round(qdac.ch01.dc_constant_V(),2)},g3={round(qdac.ch03.dc_constant_V(),2)},g5={round(qdac.ch05.dc_constant_V(),2)}"
#step_vgo_num=30+1 #sqrt(100^2+200^2)uV



step_vgo1=np.absolute((start_vgo1-stop_vgo1)/step_vgo_num)
step_vgo2=np.absolute((start_vgo2-stop_vgo2)/step_vgo_num)


#inner gate voltage range (fast axis, CS)
#####################
#start_vgi = -1.224#-0.788
#stop_vgi = -1.222#-0.776
#step_vgi_num = 30*2

start_vgi = 0.845
stop_vgi = 0.8465
step_vgi_num = 15*5#40uV
#step_vgi_num = round((stop_vgi-start_vgi)/vsd*upper_bound_lever_arm)
#print(f"step i num={step_vgi_num}")
step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)

#qdac.ramp_multi_ch_slowly([1,3,5],[0,-0.7,-0.2])
qdac.ramp_multi_ch_slowly([2,4,6],[idt_point1_y,idt_point1_x,start_vgi])

initial_guess = [(start_vgi+stop_vgi)/2, 1e-4, 5e-6]#initial guess for peakV, Gamma,height for first GVg
if start_vgi>initial_guess[0] or stop_vgi<initial_guess[0]:
    print("WARNIG:INITIAL GUESS OUT OF RANGE")

#--------Definitions-------------

#swept contacts
inner_gate=qdac.ch06.dc_constant_V  # swept gate voltage

outer_gate1=qdac.ch02.dc_constant_V
outer_gate2=qdac.ch04.dc_constant_V

#constant gate voltages, labelled by the channels they are connected to; 


#initialize constant gates, comment out for single-gate device

#qdac.ch03.dc_slew_rate_V_per_s(slew_rate)
#qdac.ch03.dc_constant_V(gate_V_ch3)
#qdac.ch05.dc_slew_rate_V_per_s(slew_rate)




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
outer_gate1.label = '5g(outer)' # Change the label of the gate chanel
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

g1_array=np.array(outer_gate1_list)
g2_array=np.array(outer_gate2_list)
delta_array=np.sqrt(abs((g1_array-start_vgo1)**2+(g2_array-start_vgo2)**2))
delta_array-=delta

inner_gate_sweep=inner_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)
measured_parameter = zurich.demods.demods0.sample

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT
vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
# applied  voltages at the intrument level before attenuation
#vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #SF:this sets up the zurich LIA
delta_current=0
delta_param = Parameter('delta', label='delta', unit='V',
                       get_cmd=lambda: delta_now)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(delta_param)
meas.register_parameter(inner_gate_sweep.parameter)   # 
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[delta_param,inner_gate_sweep.parameter])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[delta_param,inner_gate_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[delta_param,inner_gate_sweep.parameter])
meas.register_custom_parameter('outer_gate1', 'outer_gate1', unit='V', basis=[], setpoints=[delta_param,inner_gate_sweep.parameter])
meas.register_custom_parameter('outer_gate2', 'outer_gate2', unit='V', basis=[], setpoints=[delta_param,inner_gate_sweep.parameter])



# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:

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
        Glist=[]
        Vlist=[]
        Phaselist=[]
        
        
        for inner_gate_value in tqdm(inner_gate_sweep, leave=False, desc='inner gate Sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
            inner_gate_sweep.set(inner_gate_value)
            time.sleep(1.1*tc+step_vgi/slew_rate) # Wait 3 times the time contanst of the lock-in plus gate ramp speed
            measured_value = measured_parameter()
            theta_calc, v_r_calc, I,  G = zurich.phase_voltage_current_conductance_compensate(vsdac)

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
        


        datasaver.add_result(('G', Glist),
                            ('V_r', Vlist),
                            ('Phase', Phaselist),
                            (delta_param,delta_array[i-1]),
                            ('outer_gate1', [outer_gate1_sweep[i-1]]*len(fast_axis_unreversible_list)),
                            ('outer_gate2', [outer_gate2_sweep[i-1]]*len(fast_axis_unreversible_list)),
                            (inner_gate_sweep.parameter,fast_axis_unreversible_list))
        
      
        inner_gate_sweep.reverse() 
        reversed_sweep= not reversed_sweep



#print("going to sleep for the time it takes to ramp the gate plus 10 seconds")
#time.sleep(10)

#time.sleep(abs(stop_vg)/ramp_speed/1000 + 10)
print("wake up, gates are")
print(outer_gate1())
print(outer_gate2())
print(inner_gate())
#print("and source is")
#print(k2400.volt())
#plt.plot(outer_gate1_list,peakfitlist)
#plt.show()

foldername='C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part3'
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

filename=f'meas{run_id}_delta_array.npy'
path = os.path.join(foldername, filename)
np.save(path, delta_array)
