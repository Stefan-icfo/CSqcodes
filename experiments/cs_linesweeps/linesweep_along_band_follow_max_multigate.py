# weeps gates of both sections for charge sensing
# Stefan Forstner



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


#k2450 = Keithley2450
#------User input----------------
slew_rate=1e-2


   # in seconds.
tg = 5e-3 
tc = 30e-3   # in seconds.
vsd_dB = 45 # attenuation at the source in dB
vsdac = 10.9e-6 # source AC voltage in volt
device_name = 'CD11_D7_C1'
#device_name =  'CD05_G6_E3_'# 
prefix_name = "_linesweep_followmax_nobias_findgap_forall_upfrom0_withguessed_ccap0.15and0.1"

postfix = ''

# exp_name = 'Test 50 K'


mix_down_f = 1.25e6 # RLC frequency
#outer gate voltage range (slow axis, 5gate)
#####################
start_vgo1 =  0#y
stop_vgo1 =   0.5#
step_vgo_num = 200#20mV

step_vgo1=np.absolute((start_vgo1-stop_vgo1)/step_vgo_num)



#inner gate voltage range (fast axis, CS)
#####################
start_vgi = -1.76#-0.788
stop_vgi = -1.46#-0.776
step_vgi_num = 300#40uV
#step_vgi_num = round((stop_vgi-start_vgi)/vsd*upper_bound_lever_arm)
#print(f"step i num={step_vgi_num}")
step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)

start_vgi_scan=-1.52#first guess for peak
scan_range=30e-3
lower_boundary=start_vgi_scan-scan_range/2
upper_boundary=start_vgi_scan+scan_range/2
#scan_slope=-4.5e-2 #approx crosscapacitance

print(f'Scanning over {step_vgi_num*scan_range/(stop_vgi-start_vgi)} points in vgi')

#--------Definitions-------------

#swept contacts
inner_gate=qdac.ch06.dc_constant_V  # swept gate voltage

main_gate=qdac.ch01.dc_constant_V
aux_gates=[qdac.ch02.dc_constant_V,qdac.ch03.dc_constant_V,qdac.ch04.dc_constant_V,qdac.ch05.dc_constant_V]
crosscap_g1_guess=-0.15
crosscap_g2_guess=-0.1
compensation_maingate=crosscap_g1_guess*inner_gate()
aux_compensations=[crosscap_g2_guess*inner_gate(),0,0,0]
#outer_gate2=qdac.ch05.dc_constant_V
postfix = f'_crosscap_g1_guess={crosscap_g1_guess}_and_crosscap_g2_guess={crosscap_g2_guess}'


main_gate(start_vgo1)
#outer_gate2(start_vgo2)

inner_gate(start_vgi_scan-scan_range/2)
print('wait time')
#time.sleep(10)
sleeptime=max(abs(start_vgo1-main_gate()),abs(start_vgi_scan-scan_range/2-inner_gate()))/slew_rate+2
print(sleeptime)
time.sleep(sleeptime)
print("wake up, gates are")
print(main_gate())
#print(outer_auxgate1())
#print(outer_gate2())
print(inner_gate())


#freq = zurich.oscs.oscs1.freq
main_gate.label = 'main_gate_wo_cc' # Change the label of the gate chanel
inner_gate.label = 'CS(inner)' # Change the label of the source chaneel
instr_dict = dict(gate=[main_gate])
exp_dict = dict(mV = vsdac*1000)
exp_name = prefix_name+postfix
#----------- defined values------
#----------- defined values------
#####################
gain_RT = 200       #
gain_HEMT = 5.64   #
Z_tot = 7521        #
###################

# ------------------define sweep axes-------------------------

outer_gate1_sweep=main_gate.sweep(start=start_vgo1, stop=stop_vgo1, num = step_vgo_num)
#outer_gate2_sweep=outer_gate2.sweep(start=start_vgo2, stop=stop_vgo2, num = step_vgo_num)
inner_gate_sweep=inner_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)
measured_parameter = zurich.demods.demods0.sample


#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT
vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?




# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(outer_gate1_sweep.parameter)  # 
meas.register_parameter(inner_gate_sweep.parameter)  # 
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
meas.register_custom_parameter('R', 'R', unit='Ohm', basis=[], setpoints=[outer_gate1_sweep.parameter,inner_gate_sweep.parameter])
#meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])


# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:
    qdac.add_dc_voltages_to_metadata(datasaver=datasaver)

    fast_axis_unreversible_list = list(inner_gate_sweep) #(to deal with snake)
    reversed_sweep=False
    i=0
    #n=0#outer sweep count
    for outer_gate_value in tqdm(outer_gate1_sweep, leave=False, desc='outer Gate Sweep', colour = 'green'): #slow axis loop (gate)
        i=i+1#outergatesweepcounter
        #print('temperature')
        #Triton.MC()
        outer_gate1_sweep.set(outer_gate_value+compensation_maingate)
        for aux_gate,comp in zip(aux_gates,aux_compensations):
            aux_gate(outer_gate_value+comp)
        time.sleep(abs(step_vgo1/slew_rate)) # Wait  the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
        Glist=[]
        Vlist=[]
        Rlist=[]
        Phaselist=[]
        
        #print(f"lb={lower_boundary},ub={upper_boundary}")
        for inner_gate_value in tqdm(inner_gate_sweep, leave=False, desc='inner gate Sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
            if (inner_gate_value >= lower_boundary and inner_gate_value <= upper_boundary):
                inner_gate_sweep.set(inner_gate_value)
                time.sleep(1.1*tc+step_vgi/slew_rate) # Wait 3 times the time contanst of the lock-in plus gate ramp speed
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

                Glist=Glist+[G]
                Vlist=Vlist+[v_r_calc]
                Rlist=Rlist+[R]
                Phaselist=Phaselist+[theta_calc]
            else:
                Glist=Glist+[-1e-15]
                Vlist=Vlist+[-1e-15]
                Rlist=Rlist+[-1e-15]
                Phaselist=Phaselist+[-1e-15]
        #temp_fast_axis_list.reverse()
        Glist_np=np.array(Glist)
        maxid=np.argmax(Glist_np)
        V_of_max=list(inner_gate_sweep)[maxid]
        #print(f"maxid={maxid}")
        #print(f"V_of_max{V_of_max}")
        lower_boundary=V_of_max-scan_range/2
        upper_boundary=V_of_max+scan_range/2
        if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Glist.reverse()
            Vlist.reverse()
            Rlist.reverse()
            Phaselist.reverse()
            #GIVlist.reverse()
            #VRlist.reverse()
            #PHASElist.reverse()
        datasaver.add_result(('R', Rlist),
                            ('G', Glist),
                            ('V_r', Vlist),
                            ('Phase', Phaselist),
                            (outer_gate1_sweep.parameter,outer_gate_value),
                            (inner_gate_sweep.parameter,fast_axis_unreversible_list))
        
        
        inner_gate_sweep.reverse() 
        reversed_sweep= not reversed_sweep
  


#time.sleep(abs(stop_vg)/ramp_speed/1000 + 10)
print("wake up, main gate is")
print(main_gate())

print(inner_gate())
#print("and source is")
#print(k2400.volt())