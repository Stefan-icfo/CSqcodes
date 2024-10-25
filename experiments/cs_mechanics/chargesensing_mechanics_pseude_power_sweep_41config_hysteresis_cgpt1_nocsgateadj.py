# weeps gates of both sections for charge sensing
# Stefan Forstner

import numpy as np
from instruments import station,  Triton, zurich
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
import time
from tqdm import tqdm
import scipy as scp
import matplotlib.pyplot as plt
from utils.CS_utils import (
    breit_wigner_fkt,
    breit_wigner_detuning,
    zurich_phase_voltage_current_conductance,
    zurich_phase_voltage_current_conductance_compensate,
    idt_perpendicular_angle,
    make_detuning_axis,
    save_metadata_var,
    get_var_name,
    zurich_working,
)
import os
from qcodes import Parameter
import copy
import math

script_path = __file__
print("Full path of the script:", script_path)
start_time = time.time()

debug = False
# ------User input----------------
slew_rate = 1e-2

tc = 100e-3  # in seconds.
att_source_dB = 39  # attenuation at the source in dB
att_gate_dB = 46 
device_name = "CD11_D7_C1"
prefix_name = "_cs_mechanics_detune_hysteresis__now_no_attenuator!10uVdrive"

postfix = "22mK"

# Compensation values
x_avg = +3.4e-6  # +1.51e-5@75#+4.38e-6#@20mVpk -2.41e-5@100
y_avg = -5.4e-6

mix_down_f = 1.25e6  # RLC frequency
# Outer gate voltage range (slow axis, 5gate)
#####################

sit_point_g2 = -2.5126  # -1.9204
sit_point_g4 = -1.5475  # -1.8785

print(sit_point_g2, sit_point_g4)

vars_to_save = [
    slew_rate,
    tc,
    att_source_dB,
    att_gate_dB,
    x_avg,
    y_avg,
    mix_down_f,
    sit_point_g2,
    sit_point_g4,
]

# Inner gate voltage range (fast axis, CS)
#####################
start_vgi = -1.232  # -0.788
stop_vgi = -1.229  # -0.776
step_vgi_num = 3 * 100  # 20uV

step_vgi = np.absolute((start_vgi - stop_vgi) / step_vgi_num)

initial_guess = [-1.231, 1e-4, 3e-6]  # initial guess for peakV, Gamma,height for first GVg
sitfraction = 0.6  # where to sit on Coulomb peak. For now on left side

vars_to_save.extend([start_vgi, stop_vgi, step_vgi_num])
#####################
stop_f = 121.904e6  # 122e6 #Hz unit
start_f = 121.909e6  # 121.94e6 #Hz unit
step_num_f = 5 * 50  # 20Hz
num_reps = 16  # Number of repetitions
instr_magVrms = 1e-6  # Constant value for instr_magVrms


vars_to_save.extend([start_f, stop_f, step_num_f])

source_amplitude_instrumentlevel_GVg = 20e-3
source_amplitude_CNT_GVg = d2v(v2d(np.sqrt(1 / 2) * source_amplitude_instrumentlevel_GVg) - att_source_dB)
print(f"source amp at CNT for GVg:{source_amplitude_CNT_GVg*1e6} uV")
source_amplitude_instrumentlevel_mech = 20e-3
source_amplitude_CNT_mech = d2v(v2d(np.sqrt(1 / 2) * source_amplitude_instrumentlevel_mech) - att_source_dB)
print(f"source amp at CNT for mech:{source_amplitude_CNT_mech*1e6} uV")

vars_to_save.extend(
    [
        source_amplitude_instrumentlevel_GVg,
        source_amplitude_CNT_GVg,
        source_amplitude_instrumentlevel_mech,
        source_amplitude_CNT_mech,
    ]
)

# --------Definitions-------------

# Swept contacts
#cs_gate = qdac.ch06.dc_constant_V  # swept gate voltage

#outer_gate1 = qdac.ch02.dc_constant_V

#outer_gate2 = qdac.ch04.dc_constant_V

#outer_gate1(sit_point_g2)
#outer_gate2(sit_point_g4)
#cs_gate(start_vgi)

#sleeptime = (
#    max(
#        abs(sit_point_g2 - outer_gate1()),
#        abs(sit_point_g4 - outer_gate2()),
#        abs(start_vgi - cs_gate()),
#    )
#    / slew_rate
#    + 10
#)
#print(sleeptime)
#time.sleep(sleeptime)
#print("wake up, gates are")
#print(outer_gate1())
#print(outer_gate2())
#print(cs_gate())

#cs_gate.label = "CS_gate(inner)"  # Change the label of the source channel
exp_dict = dict(uVrfsource=source_amplitude_CNT_mech * 1e6)
exp_name = sample_name(prefix_name, exp_dict, postfix)

#cs_sweep = cs_gate.sweep(start=start_vgi, stop=stop_vgi, num=step_vgi_num)

freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq
freq_rlc(mix_down_f)
freq_mech(start_f)
freq_rf(start_f - mix_down_f)
time.sleep(tc)

freq_sweep = freq_rf.sweep(start=start_f, stop=stop_f, num=step_num_f)
freq_sweep_list = list(freq_sweep)
measured_parameter = zurich.demods.demods2.sample
measured_aux_parameter = zurich.demods.demods0.sample

source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value
gate_rf_enabled_param = getattr(zurich.sigouts.sigouts1.enables, f"enables{1}")

postfix = f"gate_amp_at_instr:{gate_amplitude_param()*1e6} uV"

# Introduce repetition parameter
repetition_param = Parameter(
    "repetition", label="Repetition", unit="", get_cmd=lambda: rep_num
)

# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(repetition_param)
meas.register_parameter(freq_sweep.parameter)
meas.register_custom_parameter(
    "V_rf",
    "Amplitude",
    unit="V",
    basis=[],
    setpoints=[repetition_param, freq_sweep.parameter],
)
meas.register_custom_parameter(
    "Phase",
    "Phase",
    unit="rad",
    basis=[],
    setpoints=[repetition_param, freq_sweep.parameter],
)
meas.register_custom_parameter(
    "I_rf",
    "current",
    unit="I",
    basis=[],
    setpoints=[repetition_param, freq_sweep.parameter],
)
meas.register_custom_parameter(
    "time_",
    "time_",
    unit="t",
    basis=[],
    setpoints=[repetition_param, freq_sweep.parameter],
)


with meas.run() as datasaver:
    # Saving metadata parameters

    datasaver.dataset.add_metadata("script_file", script_path)
    # Add the constant drive amplitude to metadata
    datasaver.dataset.add_metadata("start_drive", gate_amplitude_param())
    #datasaver.dataset.add_metadata("final_drive", instr_magVrms)

    # Saving metadata variables
    varnames = []
    for idx in range(len(vars_to_save)):
        varnames.append(get_var_name(vars_to_save[idx]))
    save_metadata_var(datasaver.dataset, varnames, vars_to_save)

    
    i = 0  # Outer sweep counter
    First_run = True
        
        
    last_gate_amplitude_CNT = 1  # Random value
    fast_axis_unreversible_list = list(freq_sweep)
    reversed_sweep = False
       
    for rep_num in tqdm(
            range(num_reps), leave=False, desc="Repetitions", colour="green"
        ):
        i = rep_num + 1

        time.sleep(
                tc
            )  # Wait the time it takes for the voltage to settle - doesn't quite work!

            # Lists for GVgs
        G_list = []
       
            # Lists for frequency sweeps
        I_list = []
        V_list = []
        Phase_list = []
        timelist = []

            # Run GVg to find sitpos
        if reversed_sweep == False:
                freq_rf(mix_down_f)  # Set to RLC frequency
                gate_rf_enabled_param.value(0)  # Switch off gate output
                source_amplitude_param(
                    source_amplitude_instrumentlevel_GVg
                )  # Set source amplitude
                
                

                

                #gate_amplitude_param(instr_magVrms)
                print(f"__gate at instr is {round(gate_amplitude_param()*1e6,4)} uV")
                

            #last_gate_amplitude_CNT = copy.copy(gate_amplitude_CNT)
        gate_rf_enabled_param.value(1)  # Switch on gate rf

        source_amplitude_param(source_amplitude_instrumentlevel_mech)
        freq_sweep_start_time = time.time()
        for f_value in tqdm(
                freq_sweep, leave=False, desc="Frequency Sweep", colour="green"
            ):
                freq_rf(f_value - freq_rlc())
                freq_mech(f_value)
                time.sleep(
                    1.1 * tc
                )  # Wait 1.1 times the time constant of the lock-in
                measured_value = measured_parameter()
                (
                    theta,
                    v_r,
                    I,
                    G,
                ) = zurich_phase_voltage_current_conductance(
                    measured_value, source_amplitude_CNT_mech
                )
                I_list.append(I)
                V_list.append(v_r)
                Phase_list.append(theta)
                current_time = time.time() - freq_sweep_start_time
                timelist.append(current_time)

        if reversed_sweep:
                I_list.reverse()
                V_list.reverse()
                Phase_list.reverse()
                timelist.reverse()

        datasaver.add_result(
                ("I_rf", I_list),
                ("V_rf", V_list),
                ("Phase", Phase_list),
                ("time_", timelist),
                (repetition_param, rep_num),
                (freq_sweep.parameter, fast_axis_unreversible_list),
            )

        freq_sweep.reverse()
        print(f"__reversed? {reversed_sweep}")
        reversed_sweep = not reversed_sweep
       
        First_run = False

    total_measurement_time = time.time() - start_time
    print(f"total_measurement_time {total_measurement_time}s")
    datasaver.dataset.add_metadata("total_measurement_time", total_measurement_time)
# End of measurement

gate_rf_enabled_param.value(0)


run_id = datasaver.run_id
foldername = (
    f"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_\\meas{run_id}"
)
if not os.path.exists(foldername):
    os.makedirs(foldername)
