# sweeps plunger gates diagonally while doing full GVg of charge sensor for each setpoint
# calculates sweep axis based on triple points and delta
# Stefan Forstner

#IN CONSTRUCTION



import numpy as np


#from instruments import   station, qdac,  Triton, zurich, exp
import instruments
qdac=instruments.qdac
zurich=instruments.zurich
Triton=instruments.Triton

exp=instruments.exp

from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk

import time
from tqdm import tqdm
from utils.CS_utils import *

from qcodes import Parameter

#------User input----------------

#naming
device_name = exp.device_name
prefix_name ='pt_linesweep'#

#sensitivity measurement parameters 
mod_amplitude=75e-3
mod_frequency=100e6

#s-d rf voltage at cnt
vsdac = exp.source_amplitude_CNT # source AC voltage in volt

#outer gate voltage range (slow axis, 5gate)
#####################

#idt_point2_y=-2.26365





start_vgo1=-1
stop_vgo1=1
start_vgo2=0
stop_vgo2=0.001

step_vgo_num =200+1
#or:just manually input start and stop
#start_vgo1=
#start_vgo2=
#stop_vgo1=
#stop_vgo2=
print(f"start_vgo1={start_vgo1},start_vgo2={start_vgo2},stop_vgo1={stop_vgo1},stop_vgo2={stop_vgo2}")
#inner gate voltages, CS
stop_vgi=-1.5
start_vgi=-1
step_vgi_num=500*2#exp.step_num_cs

#start_vgi = -1.645#-0.788
#stop_vgi = -1.650#41-0.776
#step_vgi_num = 500

# standard params / user input
Temp=Triton.MC()
tc = 1.1*exp.tc  # in seconds.
tg = exp.tg 
zurich.freq0(exp.freq_RLC)

step_vgo1=np.absolute((start_vgo1-stop_vgo1)/step_vgo_num)
step_vgo2=np.absolute((start_vgo2-stop_vgo2)/step_vgo_num)

vars_to_save=[Temp,tc,start_vgi,stop_vgi,step_vgo_num]

postfix = f"g1={round(qdac.ch01.dc_constant_V(),2)},g3={round(qdac.ch03.dc_constant_V(),2)},g5={round(qdac.ch05.dc_constant_V(),2)},temp={Temp:3g}"

step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)


#--------Definitions-------------

#swept contacts
inner_gate=qdac.ch06.dc_constant_V  # swept gate voltage
outer_gate1=qdac.ch01.dc_constant_V
outer_gate2=qdac.ch05.dc_constant_V

print("preramping")
qdac.ramp_multi_ch_slowly([qdac.ch01,qdac.ch05,qdac.ch06],[start_vgo1,start_vgo2,start_vgi])

print("gates are")
print(outer_gate1())
print(outer_gate2())
print(inner_gate())

#freq = zurich.oscs.oscs1.freq
outer_gate1.label = '5g(outer)' # Change the label of the gate chanel
inner_gate.label = 'CS(inner)' # Change the label of the source chaneel
#instr_dict = dict(gate=[outer_gate1])
#exp_dict = dict(mV = vsdac*1000)
exp_name = "sensitivity_RF_sweep"#sample_name(prefix_name,exp_dict,postfix)



# ------------------define sweep axes-------------------------

outer_gate1_sweep=outer_gate1.sweep(start=start_vgo1, stop=stop_vgo1, num = step_vgo_num)
outer_gate2_sweep=outer_gate2.sweep(start=start_vgo2, stop=stop_vgo2, num = step_vgo_num)
inner_gate_sweep=inner_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)
#as lists
outer_gate1_list=list(outer_gate1_sweep)
outer_gate2_list=list(outer_gate2_sweep)
#as np arrays
g1_array=np.array(outer_gate1_list)
g2_array=np.array(outer_gate2_list)

#delta array
delta_array=np.sqrt(abs((g1_array-start_vgo1)**2+(g2_array-start_vgo2)**2))
#shift delta origin to center



#------------define delta param--------------------

delta_param = Parameter('delta', label='delta', unit='V')


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(delta_param)
meas.register_parameter(inner_gate_sweep.parameter)   # 
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('I_sens', 'I_sens', unit='A', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
#meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[delta_param,inner_gate_sweep.parameter])
meas.register_custom_parameter('outer_gate1', 'outer_gate1', unit='V', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('outer_gate2', 'outer_gate2', unit='V', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])


# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:
    qdac.add_dc_voltages_to_metadata(datasaver)
    zurich.save_config_to_metadata(datasaver)
    #saving metadata variables
    varnames=[]
    for i in range(len(vars_to_save)):
        varnames.append(get_var_name(vars_to_save[i]))
    save_metadata_var(datasaver.dataset,varnames,vars_to_save)

   # fast_axis_unreversible_list = list(inner_gate_sweep) #(to deal with snake)
   # fast_axis_unreversible_array=np.array(fast_axis_unreversible_list)
    reversed_sweep=False
    i=0
    peakfitlist=[]
    First_run=True
    for outer_gate_value in tqdm(outer_gate1_sweep, leave=False, desc='outer Gate Sweep', colour = 'green'): #slow axis loop (gate)
        i=i+1#outergatesweepcounter
        qdac.ramp_multi_ch_fast([qdac.ch02,qdac.ch04], [outer_gate_value,outer_gate2_sweep[i-1]])
        outer_gate1_sweep.set(outer_gate_value)
        outer_gate2_sweep.set(outer_gate2_sweep[i-1])
        qdac.ch06.ramp_ch(start_vgi)
        time.sleep(max(abs(step_vgo1/exp.slew_rate),abs(step_vgo2/exp.slew_rate))) # Wait  the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
        
        Vglist, Glist, I_sens_list = exp.GVG_fun_sensitivity(
        start_vg=start_vgi,
        stop_vg=stop_vgi,
        step_num=step_vgi_num,
        save_in_database=False,
        return_data=True,
        return_only_Vg_G_and_Isens=True,
        reverse=False,
        sens_demod=zurich.demod2,
        RF_sens_osc=zurich.freq2,
        mod_gate=qdac.ch06,
        mod_amplitude=mod_amplitude,
        mod_frequency=mod_frequency,
        RF_meas_osc=zurich.freq0,
        RF_drive_osc=zurich.freq1,
        drive_type="RF"#LF for qdac sine wave, RF for zurich
        )
        #print(f"(delta_array[i-1])={delta_array[i-1]}")
        datasaver.add_result(('G', Glist),
                            ('I_sens', I_sens_list),
                            ('outer_gate1', [outer_gate1_sweep[i-1]]*len(Vglist)),
                            ('outer_gate2', [outer_gate2_sweep[i-1]]*len(Vglist)),
                            (delta_param,delta_array[i-1]),
                            (inner_gate_sweep.parameter,Vglist))
        
      
       # inner_gate_sweep.reverse() 
       # reversed_sweep= not reversed_sweep


#time.sleep(abs(stop_vg)/ramp_speed/1000 + 10)
print("wake up, gates are")
print(outer_gate1())
print(outer_gate2())
print(inner_gate())

