# The program does gate sweep with Zurich LIA at RLC frequency; can use multiple gates
# Stefan Forstner 
# conducance vs gate voltage using zurich

import numpy as np


from instruments import station, zurich, qdac
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk

import time
from tqdm import tqdm


#------User input----------------
gate_ramp_slope = 1e-2# V/s
tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd_dB = 45 # attenuation at the source in dB
vsdac = 200e-6 # source AC voltage in volt, for now set manually
device_name = 'CD11_D7_C1_20mK'
prefix_name = 'Conductance_qdac_zurich_testpostpcchange_wg12345at0.6,-2,1.1,-2,0.6test_3mvACongate2'
postfix = '20mK'
# exp_name = 'Test 50 K'

mix_down_f = 1.25e6 # RLC frequency
#####################
start_vg = -0.5#
stop_vg = 0#
step_num =5001 #50uV

step_vg=np.absolute((start_vg-stop_vg)/step_num)
#####################



gate=qdac.ch06



freq = zurich.oscs.oscs1.freq
gate.label = 'cs_gate' # Change the label of the gate chaneel
instr_dict = dict(gate=[gate])
exp_dict = dict(vsdac = vsdac)
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

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    # for i in range(2):
    for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='gate voltage Sweep', colour = 'green'):
        
        gate.dc_constant_V(vgdc_value)
        time.sleep(1.1*tc+step_vg/gate_ramp_slope) # Wait 3 times the time contanst of the lock-in plus gate ramp speed
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
        datasaver.add_result(('R', R),
                            ('G', G),
                            ('V_r', v_r_calc),
                            ('Phase', theta_calc),
                            (vgdc_sweep.parameter,vgdc_value))
        # vgdc_sweep.reverse()

# Ramp down everything
print(gate.dc_constant_V())
#gate(0)

#for i in range(8):
#        param1 = getattr(zurich.sigouts.sigouts0.enables, f'enables{i}')
#        param2 = getattr(zurich.sigouts.sigouts1.enables, f'enables{i}')
#        param1.value(0)
#        param2.value(0)