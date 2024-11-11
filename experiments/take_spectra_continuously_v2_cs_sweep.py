#takes zurich demod fft in bursts
#plots and saves psd


import numpy as np


from instruments import station, zurich, Triton, qdac
from qcodes.dataset import Measurement, new_experiment

from utils.CS_utils import zurich_phase_voltage_current_conductance_compensate, save_metadata_var, get_var_name

import time
from tqdm import tqdm
from qcodes import Parameter

from utils.zurich_data_fkt import *


def voltage_to_psd(v_rms, rbw, impedance=50):
  
    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd


exp_name="spectrum_vs_time_50avg_10Kfilter_208mHzBW_thermal700mK"
exp_name="spectrum_vs_time_50avg_10Kfilter_208mHzBW_drive_11uV20db_att_700mK"
#exp_name="spectrum_30mK_crosscap_g2_for_last_thermomech_at120MHz_1mVpk@instr"
device_name = 'CD11_D7_C1'

filter_bw=10e3
rbw=209.584e-3
BURST_DURATION = 4.772
#SAMPLING_RATE=13730
nr_bursts=8
reps=20
demod_ch=3

#mix_down_f = 1.25e6 # RLC frequency
#####################
gate=qdac.ch06.dc_constant_V()
start_vg = -1.231

stop_vg = -1.227
step_num= 4*100
vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_num)

freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq

device_name = f'CD11_D7_C1{freq_mech},{freq_rf}'


#vars_to_save=[gate_ramp_slope,tc,vsd_dB,source_amplitude_instrumentlevel_GVg,vsdac,x_avg,y_avg]

time_param = Parameter('time_param',
                            label='time',
                            unit='s',  # Or the appropriate unit
                            set_cmd=None,  # If there is no instrument to communicate with directly
                            get_cmd=None)  # Define get_cmd if you need to read a value
freq_param = Parameter('freq_param',
                            label='freq',
                            unit='Hz',  # Or the appropriate unit
                            set_cmd=None,  # If there is no instrument to communicate with directly
                            get_cmd=None)  # Define get_cmd if you need to read a value


gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value

# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(vgdc_sweep.parameter)  
meas.register_parameter(freq_param) 
meas.register_custom_parameter('Voltage_fft_avg', 'V_fft_avg', unit='V', basis=[], setpoints=[vgdc_sweep.parameter,freq_param])
# meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
meas.register_custom_parameter('psd', 'psd', unit='W/Hz', basis=[], setpoints=[vgdc_sweep.parameter,freq_param])



experiment_aux = new_experiment(name=exp_name+'_aux', sample_name=device_name)
meas_aux = Measurement(exp=experiment_aux)
meas_aux.register_parameter(vgdc_sweep.parameter)  
meas_aux.register_parameter(freq_param)
meas_aux.register_custom_parameter('Voltage_fft', 'V_fft', unit='V', basis=[], setpoints=[vgdc_sweep.parameter,freq_param])
meas_aux.register_custom_parameter('Voltage_fft_log', 'V_fft_log', unit='logV', basis=[], setpoints=[vgdc_sweep.parameter,freq_param])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)


# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:
    datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V_start',qdac.ch01.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch02_dc_constant_V_start',qdac.ch02.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V_start',qdac.ch03.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch04_dc_constant_V_start',qdac.ch04.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch05_dc_constant_V_start',qdac.ch05.dc_constant_V())
    datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V_start',qdac.ch06.dc_constant_V())
    datasaver.dataset.add_metadata('probe_freq',freq_rf())
    datasaver.dataset.add_metadata('rlc_freq',freq_rlc())
    datasaver.dataset.add_metadata('center_freq',freq_mech())
    datasaver.dataset.add_metadata('filter_bw',filter_bw)
    datasaver.dataset.add_metadata('rbw',rbw)
    datasaver.dataset.add_metadata('nr_bursts',nr_bursts)
    datasaver.dataset.add_metadata('nr_reps',reps)
    gate_amplitude_value=gate_amplitude_param()
    datasaver.dataset.add_metadata('gateampatinstr',gate_amplitude_value)

    with meas_aux.run() as datasaver_aux:
        varnames=[]
    #or i in range(len(vars_to_save)):
    #    varnames.append(get_var_name(vars_to_save[i]))
    #save_metadata_var(datasaver.dataset,varnames,vars_to_save)
    # for i in range(2):
        for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='gate voltage Sweep', colour = 'green'):
            full_data, averaged_data_per_burst, averaged_data, freq,filter_data  = take_spectrum(demod_ch)    
            meas_time=0
            for data,avg_data,filter in zip(full_data,averaged_data_per_burst,filter_data):
                logdata=np.log(data)
                #compensated_data=data/filter
                #compensated_avg_data=data/filter
                datasaver_aux.add_result(('Voltage_fft', data),
                                    ('Voltage_fft_log', logdata),
                                    (vgdc_sweep.parameter,vgdc_value),
                                    (freq_param,freq))
                
                target_size = np.shape(avg_data)[0]
                factor = len(freq) // target_size  # Factor by which to compress

                # Reshape the array and compute the mean along the compressed axis
                compressed_freq = np.mean(freq[:target_size*factor].reshape(-1, factor), axis=1)
                avg_data_psd=voltage_to_psd(avg_data, rbw)
                datasaver.add_result(('Voltage_fft_avg', avg_data),
                                     ('psd', avg_data_psd),
                                    (time_param,meas_time),
                                    (freq_param,compressed_freq))

                meas_time+=BURST_DURATION*nr_bursts



