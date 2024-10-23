# Sweeps gate voltages while reading out with Zurich LIA 
# Sweeps two gates while keeping others constant
# Stefan Forstner using template by Parameshwar Prasad

# snake is working now. watch out - there might still be issues when setting step_ramp_speed to slow and the step_size too large!

import numpy as np
import os

from instruments import  manual, station, qdac,zurich
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name


from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
import time
from tqdm import tqdm


#------User input----------------
ramp_source = 1e-6 # V/ms
step_v = 1e-6 # source steps
ramp_speed = 0.01 # V/s for large ramps
step_ramp_speed=0.1 # between steps, V/s
tc = 10e-3   # in seconds. 
vsdac = 100e-6 # source DC voltage in volt
vsd_dB = 45 
device_name = 'CD11_D7_C1_all5g'
prefix_name = 'Charge_stability_k2400_QDev'
postfix = '20mK_constantgates135at-3+1-3_zoom'
#offset = -10e-6 #voltage offset of k2400
#offset_i=-44e-12



#outer voltage range (slow axis)
#####################
start_vg1 = -1.48    #
stop_vg1 = -1.465      #
step_vg1_num = 30+1     #
step_vg1=np.absolute((start_vg1-stop_vg1)/step_vg1_num)


#inner voltage range (fast axis)
#####################
start_vg2 = -1.50    #
stop_vg2 = -1.40     #
step_vg2_num = 50*10+1    #
step_vg2=np.absolute((start_vg2-stop_vg2)/step_vg2_num)

start_vgcs=-0.1295 #
step_cs_num=50*10+1#
delta=5e-3#10mV
#constant gate voltages, labelled by the channels they are connected to; 
gate_V_ch3=+1
gate_V_ch1=-3
gate_V_ch5=-3

#initialize constant gates, comment out for single-gate device

qdac.ch03.dc_slew_rate_V_per_s(ramp_speed)
qdac.ch03.dc_constant_V(gate_V_ch3)
qdac.ch05.dc_slew_rate_V_per_s(ramp_speed)
qdac.ch05.dc_constant_V(gate_V_ch5)
qdac.ch01.dc_slew_rate_V_per_s(ramp_speed)
qdac.ch01.dc_constant_V(gate_V_ch1)
# qdac.ch06.dc_slew_rate_V_per_s(ramp_speed)
# qdac.ch06.dc_constant_V(gate_V_ch6)
# qdac.ch07.dc_slew_rate_V_per_s(ramp_speed)
# qdac.ch07.dc_constant_V(gate_V_ch7)

#--------Definitions-------------

#swept contacts
gate1=qdac.ch02
  # swept outer gate voltage
gate2=qdac.ch04 #swept inner gate voltage
#source = k2400 # source 
csgate=qdac.ch06

gate1.label = 'gate2' # Change the label of the gate1 chanel
gate2.label = 'gate4' # Change the label of the gate2 chaneel
instr_dict = dict(gate1=[gate1])
exp_dict = dict(vsdac = vsdac)
exp_name = sample_name(prefix_name,exp_dict,postfix)

#----------- defined values------#----------- defined values------
#####################
gain_RT = 200       #
gain_HEMT = 5.64   #
Z_tot = 7521        #
###################


# ------------------define sweep axes-------------------------

gate1_sweep=gate1.dc_constant_V.sweep(start=start_vg1, stop=stop_vg1, num = step_vg1_num)
gate2_sweep=gate2.dc_constant_V.sweep(start=start_vg2, stop=stop_vg2, num = step_vg2_num)
g1sweeplist=list(gate1_sweep)
g2sweeplist=list(gate2_sweep)
  # lock-in amplitude measured # SF:CHANGED FROM 'zurich.source.demod_complex'

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT

# applied  voltages at the intrument level before attenuation




#initialize swept contacts

#slow ramp and intial voltage
gate1.dc_slew_rate_V_per_s(ramp_speed)
gate1.dc_constant_V(start_vg1)

gate2.dc_slew_rate_V_per_s(ramp_speed)
gate2.dc_constant_V(start_vg2)

csgate.dc_slew_rate_V_per_s(ramp_speed)
csgate.dc_constant_V(start_vgcs-delta)


#time.sleep(max([abs(start_vg1/ramp_speed),abs(start_vg2/ramp_speed)])+1)  #wait for the time it takes to do both ramps plus one second
measured_parameter = zurich.demods.demods0.sample
vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))  
time.sleep(30)
print("gate channels")
print(qdac.ch01.dc_constant_V())
print(qdac.ch02.dc_constant_V())
print(qdac.ch03.dc_constant_V())
print(qdac.ch04.dc_constant_V())
print(qdac.ch05.dc_constant_V())
print(qdac.ch06.dc_constant_V())
print(qdac.ch07.dc_constant_V())
#set fast ramp speeds
gate1.dc_slew_rate_V_per_s(step_ramp_speed)
gate2.dc_slew_rate_V_per_s(step_ramp_speed)
csgate.dc_slew_rate_V_per_s(step_ramp_speed)


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(gate1_sweep.parameter)  # 
meas.register_parameter(gate2_sweep.parameter)  # 
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('R', 'R', unit='Ohm', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])


meas.register_custom_parameter('peak_position', 'V_peak', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('peak_Value', 'G_peak', unit='S', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])

meas.register_custom_parameter('sit_position', 'V', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('top_bound', 'G_top', unit='S', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('bottom_bound', 'G_bottom', unit='S', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])

AllData3d_conductance=np.zeros((step_vg1_num,step_vg2_num,step_cs_num))
AllData3d_pp=np.zeros((step_vg1_num,step_vg2_num,step_cs_num))



# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:

    fast_axis_unreversible_list = list(gate2_sweep) #(to deal with snake)
    reversed_sweep=False
    peakpos=start_vgcs
    
    for gate1_value in tqdm(gate1_sweep, leave=False, desc='outer gate sweep', colour = 'green'): #slow axis loop (gate)
        gate1_sweep.set(gate1_value)
        time.sleep(1*tc+step_vg1/step_ramp_speed) # Wait 3 times the time contanst of the lock-in plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
        #init lists (to deal with snake)
        Glist=[]
        Vlist=[]
        Rlist=[]
        Phaselist=[]
        peakpos_list=[]
        peakG_list=[]
        
 
        # temp_fast_axis_list=[] #I think this line is unnecessary
        #following is necessary only if no snake #COMMENT OUT ONCE SNAKE IS WORKING
        #gate2.dc_slew_rate_V_per_s(ramp_speed)
        #gate2.dc_constant_V(start_vg2)
        #time.sleep(abs((start_vg2-stop_vg2)/ramp_speed)) 
        #gate2.dc_slew_rate_V_per_s(step_ramp_speed)

        for gate2_value in tqdm(gate2_sweep, leave=False, desc='inner gate sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
            gate2_sweep.set(gate2_value)
            time.sleep(1.1*tc+step_vg2/step_ramp_speed) # Wait 3 times the time constant of the lock-in, plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
            cs_sweep=csgate.dc_constant_V.sweep(start=peakpos-delta, stop=peakpos+delta, num = step_cs_num)
            cssweepcondlist=[]
            csgate.dc_constant_V(start_vgcs-delta)
            time.sleep(2*delta/step_ramp_speed+1)
            
            for gatecs_value in (cs_sweep):
                csweeplist=list(cs_sweep)
                cs_sweep.set(gatecs_value)
                time.sleep(1.1*tc+2*delta/step_cs_num/step_ramp_speed)
                measured_value = measured_parameter()
                x = measured_value['x'][0] #SF: COMMENTED OUT 
                y = measured_value['y'][0]#SF: COMMENTED OUT
                #xy_complex = measured_value
                xy_complex = complex(x,y)
                v_r_calc = np.absolute(xy_complex)
                theta_calc = np.angle(xy_complex)
                    
            #G calculation
                I = v_r_calc/(gain_RT*gain_HEMT*Z_tot)
                G = 1/((vsdac/I)-Z_tot)
                R = 1/G

                AllData3d_conductance[g1sweeplist.index(gate1_value),g2sweeplist.index(gate2_value),csweeplist.index(gatecs_value)]=G
                AllData3d_conductance[g1sweeplist.index(gate1_value),g2sweeplist.index(gate2_value),csweeplist.index(gatecs_value)]=gatecs_value
                cssweepcondlist=cssweepcondlist+[G]

            peak_G=max(cssweepcondlist)
            peakG_list=peakG_list+[peak_G]
            
            peakpos=csweeplist[cssweepcondlist.index(max(cssweepcondlist))]

            Glist=Glist+[G]
            Vlist=Vlist+[v_r_calc]
            Rlist=Rlist+[R]
            Phaselist=Phaselist+[theta_calc]
            peakpos_list=peakpos_list+[peakpos]
            
            
        #Rlist.reverse
        if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Glist.reverse()
            Vlist.reverse()
            Rlist.reverse()
            Phaselist.reverse()
            peakG_list.reverse()
            peakpos_list.reverse()
            #GIVlist.reverse()
            #VRlist.reverse()
            #PHASElist.reverse()
        datasaver.add_result(('R', Rlist),
                            ('G', Glist),
                            ('V_r', Vlist),
                            ('Phase', Phaselist),
                            ('peak_position',peakpos_list),
                            ('peak_Value',peakG_list),
                            (gate1_sweep.parameter,gate1_value),
                            (gate2_sweep.parameter,fast_axis_unreversible_list))
          
        #temp_fast_axis_list = list(gate2_sweep) #(to deal with snake)
        #temp_fast_axis_list.reverse
        #datasaver.add_result(('Conductance',Glist),('Resistance', Rlist),('Current',Ilist),
         #                   (gate1_sweep.parameter,gate1_value),
         #                   (gate2_sweep.parameter,fast_axis_unreversible_list))
        gate2_sweep.reverse()
        reversed_sweep= not reversed_sweep 

# Ramp down everything
gate1.dc_slew_rate_V_per_s(ramp_speed)
gate2.dc_slew_rate_V_per_s(ramp_speed)



