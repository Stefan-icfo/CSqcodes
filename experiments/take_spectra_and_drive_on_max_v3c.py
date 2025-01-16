#takes zurich demod fft in bursts
#plots and saves psd


import numpy as np
import scipy as scp


from instruments import station, zurich, Triton, qdac
from qcodes.dataset import Measurement, new_experiment

from utils.CS_utils import *

import time
from tqdm import tqdm
from qcodes import Parameter
import copy

from utils.zurich_data_fkt import *

from experiments.Do_GVg_and_adjust_sitpos import do_GVg_and_adjust_sitpos

#exp_name="spectrum_vs_time_50avg_10Kfilter_208mHzBW_thermal30mK"
#exp_name="spectrum_vs_time_50avg_10Kfilter_208mHzBW_drive_1.1uV20db_att_30mK"
#exp_name="spectrum_30mK_crosscap_g2_for_last_thermomech_at120MHz_1mVpk@instr"

device_name = 'CD11_D7_C1_100mK_'

filter_bw=10e3
rbw=209.584e-3
BURST_DURATION = 4.772
#SAMPLING_RATE=13730
nr_bursts=8
reps=4
reps_nodrive=5
reps_drive=5
demod_ch=3
drive_offset=0



freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq


def voltage_to_psd(v_rms, rbw, impedance=50):
  
    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd

def take_long_spectra(reps=reps,demod_ch=demod_ch):
   # zurich.set_frequencies_to_json_config("160MHz_squeezed_singledot2")
   # print("JUST SET BACK FREQUENCIES")
    meas_time=0
    datas,avg_datas,avg_datas_psd,meas_times=[],[],[],[]
    for n in tqdm(range(reps)):
            full_data, averaged_data_per_burst, averaged_data, freq,compressed_freq,filter_data  = take_spectrum(demod_ch)  


            for data,avg_data in zip(full_data,averaged_data_per_burst):
                datas.append(data)
                avg_datas.append(avg_data)
                avg_data_psd=voltage_to_psd(avg_data, rbw)#calculate psd
                avg_datas_psd.append(avg_data_psd)
                meas_times.append(meas_time)

                
                meas_time+=BURST_DURATION*nr_bursts

            #calculate compressed frequency axis, assuming freq is always the same  
            #target_size = np.shape(avg_data)[0]
            #factor = len(freq) // target_size  # Factor by which to compress
            #compressed_freq = np.mean(freq[:target_size*factor].reshape(-1, factor), axis=1)  # Reshape the array and compute the mean along the compressed axis

    values_to_return={'Voltage_fft': np.array(datas),'Voltage_fft_avg' : np.array(avg_datas), 'avg_psd' : np.array(avg_datas_psd), "freq": np.array(freq), "compressed_freq" : np.array(compressed_freq), "meas_times" : np.array(meas_times)}
    
    return values_to_return

#vars_to_save=[gate_ramp_slope,tc,vsd_dB,source_amplitude_instrumentlevel_GVg,vsdac,x_avg,y_avg]

from experiments.cs_experiment import CSExperiment
temp_meas_100mK=CSExperiment()

def run_thermomech_temp_meas():
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
    exp_name=f"1dot_drive_nodrive_spectrum_{gate_amplitude_param()*1e6}uVdriveatinstr"
    # ----------------Create a measurement-------------------------
    experiment = new_experiment(name=exp_name, sample_name=device_name)
    meas = Measurement(exp=experiment)
    meas.register_parameter(time_param)  
    meas.register_parameter(freq_param) 
    #meas.register_custom_parameter('Voltage_fft_avg', 'V_fft_avg', unit='V', basis=[], setpoints=[time_param,freq_param])
    # meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
    meas.register_custom_parameter('avg_psd', 'avg_psd', unit='W/Hz', basis=[], setpoints=[time_param,freq_param])



    #experiment_fulldata = new_experiment(name=exp_name+'_full_data', sample_name=device_name)
    #meas_aux = Measurement(exp=experiment_fulldata)
    #meas_aux.register_parameter(time_param)  
    #meas_aux.register_parameter(freq_param)
    #meas_aux.register_custom_parameter('Voltage_fft', 'V_fft', unit='V', basis=[], setpoints=[time_param,freq_param])
    #meas_aux.register_custom_parameter('Voltage_fft_log', 'V_fft_log', unit='logV', basis=[], setpoints=[time_param,freq_param])

    # meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)
    experiment_1D = new_experiment(name=exp_name+'_1D', sample_name=device_name)
    meas_aux_aux = Measurement(exp=experiment_1D)
    #meas_aux.register_parameter(time_param)  
    meas_aux_aux.register_parameter(freq_param)
    meas_aux_aux.register_custom_parameter('avg_avg_psd_nodrive', 'avg_avg_psd_nodrive', unit='W/Hz', basis=[], setpoints=[freq_param])
    meas_aux_aux.register_custom_parameter('avg_avg_psd_drive', 'avg_avg_psd_drive', unit='W/Hz', basis=[], setpoints=[freq_param])
    meas_aux_aux.register_custom_parameter('avg_avg_psd_nodrive_scaled', 'avg_avg_psd_nodrive_scaled', unit='a.u.', basis=[], setpoints=[freq_param])
    meas_aux_aux.register_custom_parameter('avg_avg_psd_nodrive_w_driven_value', 'avg_avg_psd_nodrive_w_driven_value', unit='W/Hz', basis=[], setpoints=[freq_param])

    # # -----------------Start the Measurement-----------------------

    with meas.run() as datasaver:
        #saving metadata parameters
        qdac.add_dc_voltages_to_metadata(datasaver)
        zurich.save_config_to_metadata(datasaver)
        datasaver.dataset.add_metadata('probe_freq',freq_rf())
        datasaver.dataset.add_metadata('rlc_freq',freq_rlc())
        datasaver.dataset.add_metadata('center_freq',freq_mech())
        datasaver.dataset.add_metadata('filter_bw',filter_bw)
        datasaver.dataset.add_metadata('rbw',rbw)
        datasaver.dataset.add_metadata('nr_bursts',nr_bursts)
        datasaver.dataset.add_metadata('nr_reps',reps)
        datasaver.dataset.add_metadata('nr_reps_nodrive',reps_nodrive)
        datasaver.dataset.add_metadata('nr_reps_drive',reps_drive)
        gate_amplitude_value=gate_amplitude_param()
        datasaver.dataset.add_metadata('gateampatinstr',gate_amplitude_value)

        with meas_aux_aux.run() as datasaver_aux_aux:
            varnames=[]
            zurich.set_frequencies_to_json_config("160MHz_squeezed_singledot2")
            print("JUST SET BACK FREQUENCIES")

            zurich.sigout1_amp1_enabled_param.value(0)
            returned_values_nodrive=take_long_spectra(reps=reps_nodrive,demod_ch=demod_ch)
            
            

            #read vslues needed for saving anc calculation
            avg_psd_array_nodrive=returned_values_nodrive['avg_psd']
            compressed_freq_array=returned_values_nodrive["compressed_freq"]
            meas_times_nodrive=returned_values_nodrive['meas_times']


            #now save these values
            for m_time,avg_psd in zip(meas_times_nodrive,avg_psd_array_nodrive):
                datasaver.add_result(#('Voltage_fft_avg', returned_values_nodrive['Voltage_fft_avg']),
                                        ('avg_psd', avg_psd),
                                        (freq_param,compressed_freq_array),
                                        (time_param,m_time))
            #time the nodrive measurement ended, in burst time -needed for next datasave
            switch_time=meas_times_nodrive[-1]

            
            avg_avg_psd_nodrive=np.mean(avg_psd_array_nodrive,axis=0)
            
            #now calculate peak frequency
            #max_relative_freq = freq[np.argmax(averaged_data)]
            max_relative_freq=compressed_freq_array[np.argmax(avg_avg_psd_nodrive)]#offset from zero frequency of demodulator
            #empiric:
            #max_relative_freq=max_relative_freq/2

            freq_rlc_value=freq_rlc()#mixdown frequency
            freq_rf_value=freq_rf()#source drive frequency, mech frequency minus (convention) mixdown frequency
            print(f"frequency shift for drive:{max_relative_freq}")
            freq_mech(freq_rf_value+freq_rlc_value+max_relative_freq)#drive at maximum
            freq_rf(freq_rf_value+max_relative_freq)#drive at maximum
            zurich.sigout1_amp1_enabled_param.value(1)
            print("drive on")

            returned_values_drive=take_long_spectra(reps=reps_drive,demod_ch=demod_ch)
            #now save average values
            meas_times_drive=returned_values_drive['meas_times']+meas_times_nodrive+BURST_DURATION
            avg_psd_array_drive=returned_values_drive['avg_psd']

            for m_time,avg_psd in zip(meas_times_nodrive,avg_psd_array_drive):
                datasaver.add_result(#('Voltage_fft_avg', returned_values_drive['Voltage_fft_avg']),
                                        ('avg_psd', avg_psd),
                                        (time_param,m_time),
                                        (freq_param,compressed_freq_array))#if the center frequency of demod 3 was changed, this would have to be adapted
            
            #now focus on non-averaged data for sharp drive peak
            avg_v_array_driven=returned_values_drive['Voltage_fft']
            #average along time axis, make sure axis is choden right
            avg_avg_v_driven=np.mean(avg_v_array_driven,axis=0)

            driven_value_narrowband=voltage_to_psd(max(avg_avg_v_driven), rbw)
            drive_difference_narrowband=driven_value_narrowband-max(avg_avg_psd_nodrive)

            #for testing,etc
            avg_driven_psd_array=np.array(returned_values_drive['avg_psd'])
            avg_avg_driven_psd=np.mean(avg_driven_psd_array,axis=0)
            #now plot for testing purposes
            """
            X, Y = np.meshgrid(compressed_freq_array, meas_times_nodrive, indexing='ij')
            plt.pcolor(X,Y,avg_psd_array_nodrive.T)
            plt.title("driven psd vs time")
            plt.show()#now plot for testing purposes

            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.plot(max_relative_freq,1.1*max(avg_avg_psd_nodrive),'g*')
            plt.title("nondriven psd avg and positon of maximum")
            plt.show()




            X, Y = np.meshgrid(compressed_freq_array, meas_times_drive, indexing='ij')
            plt.pcolor(X,Y,avg_driven_psd_array.T)
            plt.title("driven psd vs time")
            plt.show()
            
            
            plt.plot(compressed_freq_array,avg_avg_driven_psd)
            plt.title("driven psd and max pos/driving pos")
            plt.plot(max_relative_freq,1.1*max(avg_avg_psd_nodrive),'g*')
            plt.show()
            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.title("nondriven avg psd")
            plt.show()
            
            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.title("narrowband driven avg psd together with non-driven psd")
            plt.plot(max_relative_freq,driven_value_narrowband,'g*')
            plt.show()
            """
            #now fit lorentzian to scaled value
            Gamma_guess=0.5e3
            offset_approx=1.5e-15
            initial_guess=[max_relative_freq,Gamma_guess,max(avg_avg_psd_nodrive),min(avg_avg_psd_nodrive)]
            popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, compressed_freq_array, avg_avg_psd_nodrive, p0=initial_guess)
            """
            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.title("Lorentzian fit initial guess")
            plt.plot(compressed_freq_array,lorentzian_fkt(compressed_freq_array,initial_guess[0],initial_guess[1],initial_guess[2],initial_guess[3]))
            plt.show()
            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.title("Lorentzian fit")
            plt.plot(compressed_freq_array,lorentzian_fkt(compressed_freq_array,popt[0],popt[1],popt[2],popt[3]))
            plt.show()
            """
            lorentzian, area_under_lorentzian=lorentzian_fkt_w_area(compressed_freq_array,popt[0],popt[1],popt[2],popt[3])

            avg_avg_psd_nodrive_with_driven_value=copy.copy(avg_avg_psd_nodrive)
            closest_index = np.argmin(np.abs(compressed_freq_array - max_relative_freq))
            avg_avg_psd_nodrive_with_driven_value[closest_index]=driven_value_narrowband

            datasaver_aux_aux.add_result(('avg_avg_psd_nodrive',avg_avg_psd_nodrive),
                                        ('avg_avg_psd_nodrive_w_driven_value',avg_avg_psd_nodrive_with_driven_value),
                                        ('avg_avg_psd_drive',avg_avg_driven_psd),
                                        ('avg_avg_psd_nodrive_scaled',avg_avg_driven_psd/drive_difference_narrowband),
                                        (freq_param,compressed_freq_array))
            
            


            zurich.sigout1_amp1_enabled_param.value(0)
            print("drive off")

            print(f"driven_value_narrowband {driven_value_narrowband}")
            print(f"drive_difference_narrowband {drive_difference_narrowband}")
            print(f"max(avg_avg_psd) {max(avg_avg_psd_nodrive)}")
            print(f"area_under_lorentzian {area_under_lorentzian}")
            print(f"area_under_lorentzian scaled {area_under_lorentzian/drive_difference_narrowband}")

            datasaver.dataset.add_metadata('driven_value_narrowband',driven_value_narrowband)
            datasaver.dataset.add_metadata('drive_difference_narrowband',drive_difference_narrowband)
            datasaver.dataset.add_metadata('max_avg_avg_psd_',max(avg_avg_psd_nodrive))
            datasaver.dataset.add_metadata('area_under_lorentzian',area_under_lorentzian)
            datasaver.dataset.add_metadata('area_under_lorentzian_scaled',area_under_lorentzian/drive_difference_narrowband)
            datasaver.dataset.add_metadata('freq_mech_corrected',freq_mech())
            datasaver.dataset.add_metadata('freq_rf_',freq_rf_value)

#now import instance of exp_class

#from experiments.cs_experiment import CSExperiment

#temp_meas_180mK=CSExperiment()\
        

        temp_meas_100mK.area_values_scaled.append(area_under_lorentzian/drive_difference_narrowband)
        temp_meas_100mK.area_values_unscaled.append(area_under_lorentzian)
do_GVg_and_adjust_sitpos()
for n in range(10):
    #do_GVg_and_adjust_sitpos()
    run_thermomech_temp_meas()
    print(f"overall rep nr {n}")