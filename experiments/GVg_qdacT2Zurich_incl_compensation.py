# The program does gate sweep with Zurich LIA at RLC frequency; can use multiple gates
# Stefan Forstner 
# conducance vs gate voltage using zurich
#IN CONSTRUCTION

import numpy as np


from instruments import station, zurich, qdac,Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
from utils.CS_utils import zurich_phase_voltage_current_conductance_compensate, save_metadata_var, get_var_name

import time
from tqdm import tqdm




#------User input----------------
gate_ramp_slope = 1e-2# V/s
tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd_dB = 39 # attenuation at the source in dB
source_amplitude_instrumentlevel_GVg = 20e-3
vsdac=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_GVg)-vsd_dB)/10

print(f"source amp at CNT for GVg:{vsdac*1e6} uV")
device_name = 'CD11_D7_C1'
prefix_name = 'Conductance_rf_'
#postfix = '25mK_'
# exp_name = 'Test 50 K'
#Temp=Triton.MC()
postfix = f"_g1={round(qdac.ch01.dc_constant_V(),2)},g2={round(qdac.ch02.dc_constant_V(),2)},g3={round(qdac.ch03.dc_constant_V(),2)},g4={round(qdac.ch04.dc_constant_V(),2)},g5={round(qdac.ch05.dc_constant_V(),2)}"
#postfix = f"_5g={round(qdac.ch05.dc_constant_V(),7)}"
#postfix = f"_{zurich.sigouts.sigouts0.amplitudes.amplitudes0.value()}mVrf_pk"
#postfix = f"_100mVrf_pk"
x_avg=+3.4e-6  #+1.51e-5@75#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-5.4e-6  #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100

vars_to_save=[gate_ramp_slope,tc,vsd_dB,source_amplitude_instrumentlevel_GVg,vsdac,x_avg,y_avg]

mix_down_f = 1.25e6 # RLC frequency
#####################
start_vg = -2.232

stop_vg = -2.229
step_num= 30*20
#

step_vg=np.absolute((start_vg-stop_vg)/step_num)
#####################



gate=qdac.ch06




freq = zurich.oscs.oscs0.freq
freq(mix_down_f)
source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
source_amplitude_param(source_amplitude_instrumentlevel_GVg)#set source aplitude
gate.label = 'cs_gate' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(mV = vsdac*1000)
exp_name = sample_name(prefix_name,exp_dict,postfix)

#----------- defined values------
#####################
gain_RT = 200       #
gain_HEMT = 5.64   #%
Z_tot = 7521        #
###################


# ------------------Create a new Experiment-------------------------

#vgdc_sweep = gate.sweep(start=start_vg, stop=stop_vg, num = step_num)  # gate parameter and the instrument will go here
vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_num)
measured_parameter = zurich.demods.demods0.sample   # lock-in amplitude measured # SF:CHANGED FROM 'VOLTAGE'

#------------init--------------------
#manual.vsd_attn(vsd_dB) # SF: COMMENTED OUT
# applied  voltages at the intrument level before attenuation
vsdac0 = rms2pk(d2v(v2d(vsdac)+vsd_dB))   #what is this?
#zi_uhfli_GVg_setup(vsdac0,mix_down_f,tc)  #what is this?
#bilt.channels.output_mode('ramp')
#bilt.channels.ramp_slope(ramp_speed)


gate.dc_slew_rate_V_per_s(gate_ramp_slope)
gate.dc_constant_V(start_vg)
print(f"going to sleep for the time it takes to ramp the gate({abs(start_vg-gate.dc_constant_V())/gate_ramp_slope}) plus 30 seconds")
#time.sleep(20)
time.sleep(abs(start_vg-gate.dc_constant_V())/gate_ramp_slope + 30)
print("wake up, gate is")
print(gate.dc_constant_V())



# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(vgdc_sweep.parameter)  # register1the 1st1indepe n 10e-3ent parameter
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('R', 'R', unit='Ohm', basis=[], setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('I', 'Current', unit='A', basis=[], setpoints=[vgdc_sweep.parameter])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V_start',qdac.ch01.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch02_dc_constant_V_start',qdac.ch02.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V_start',qdac.ch03.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch04_dc_constant_V_start',qdac.ch04.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch05_dc_constant_V_start',qdac.ch05.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V_start',qdac.ch06.dc_constant_V())
    varnames=[]
    for i in range(len(vars_to_save)):
        varnames.append(get_var_name(vars_to_save[i]))
    save_metadata_var(datasaver.dataset,varnames,vars_to_save)
    # for i in range(2):
    for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='gate voltage Sweep', colour = 'green'):
        
        gate.dc_constant_V(vgdc_value)
        time.sleep(1.1*tc+step_vg/gate_ramp_slope) # Wait 3 times the time contanst of the lock-in plus gate ramp speed
        measured_value = measured_parameter()
        theta_calc, v_r_calc, I,  G = zurich_phase_voltage_current_conductance_compensate(measured_value, vsdac,x_avg, y_avg)
        R = 1/G
        datasaver.add_result(('R', R),
                            ('G', G),
                            ('V_r', v_r_calc),
                            ('Phase', theta_calc),
                            ('I', I),
                            (vgdc_sweep.parameter,vgdc_value))
        # vgdc_sweep.reverse()

# Ramp down everything
print(gate.dc_constant_V())
#gate.dc_constant_V(0)

#for i in range(8):
#        param1 = getattr(zurich.sigouts.sigouts0.enables, f'enables{i}')
#        param2 = getattr(zurich.sigouts.sigouts1.enables, f'enables{i}')
#        param1.value(0)
#        param2.value(0)
