#takes zurich demod fft in bursts
#plots and saves psd


import numpy as np
import scipy as scp
import os


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
local_time=time.localtime()#test
formatted_time = f"{local_time.tm_mday:02}:{local_time.tm_hour:02}:{local_time.tm_min:02}"
print("Current time (dd:hh:mm):", formatted_time)


time.sleep(10) 
device_name = 'CD11_D7_C1_mapi_sensitivity_meas'
from experiments.cs_experiment import *
temp_meas_fluctuating_base_mK2=thermomech_measurement()

filter_bw=10e3
rbw=209.584e-3#0.808190#209.584e-3
BURST_DURATION =4.772#1.193#0.569523# 4.772
SAMPLING_RATE = 13730#54.93e3
#SAMPLING_RATE=13730
nr_bursts=7
#reps=4
reps_nodrive=4
#reps_drive=20
demod_ch=3
drive_offset=0

slope=1#if not measuring it

freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq


def voltage_to_psd(v_rms, rbw, impedance=50):
  
    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd

def take_long_spectra(reps,demod_ch=demod_ch):
   # zurich.set_frequencies_to_json_config("160MHz_squeezed_singledot2")
   # print("JUST SET BACK FREQUENCIES")
    meas_time=0
    datas,avg_datas,avg_datas_psd,meas_times=[],[],[],[]
    for n in tqdm(range(reps)):
            full_data, averaged_data_per_burst, averaged_data, freq,compressed_freq,filter_data  = take_spectrum(demod_ch,SAMPLING_RATE = SAMPLING_RATE)  
            #freq_real=freq+freq_mech()
            compressed_freq_real=compressed_freq+freq_mech()

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

    values_to_return={'Voltage_fft': np.array(datas),'Voltage_fft_avg' : np.array(avg_datas), 'avg_psd' : np.array(avg_datas_psd), "freq": np.array(freq), "compressed_freq" : np.array(compressed_freq), "meas_times" : np.array(meas_times),"compressed_freq_real" : np.array(compressed_freq_real),'filter':np.array(filter_data)}
    
    return values_to_return

#vars_to_save=[gate_ramp_slope,tc,vsd_dB,source_amplitude_instrumentlevel_GVg,vsdac,x_avg,y_avg]


gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value
gate_amplitude_value = gate_amplitude_param()
def run_thermomech_temp_meas(reps_nodrive=reps_nodrive):

    ###########################################3
    local_time=time.localtime()
    formatted_time = f"{local_time.tm_mday:02}:{local_time.tm_hour:02}:{local_time.tm_min:02}"
    print("day and time before gvg (dd:hh:mm):", formatted_time)
    

    #slope,sitpos=do_GVg_and_adjust_sitpos(testplot=True)
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


    
    gate_amp_uV=gate_amplitude_param()*1e6
    exp_name="1dot_nodrive_spectrum"
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
    #meas_aux_aux.register_custom_parameter('avg_avg_psd_drive', 'avg_avg_psd_drive', unit='W/Hz', basis=[], setpoints=[freq_param])
    #meas_aux_aux.register_custom_parameter('avg_avg_psd_nodrive_scaled', 'avg_avg_psd_nodrive_scaled', unit='a.u.', basis=[], setpoints=[freq_param])
    #meas_aux_aux.register_custom_parameter('avg_avg_psd_nodrive_w_driven_value', 'avg_avg_psd_nodrive_w_driven_value', unit='W/Hz', basis=[], setpoints=[freq_param])

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
        #datasaver.dataset.add_metadata('nr_reps',reps)
        datasaver.dataset.add_metadata('nr_reps_nodrive',reps_nodrive)
        #datasaver.dataset.add_metadata('nr_reps_drive',reps_drive)
        gate_amplitude_value=gate_amplitude_param()
        datasaver.dataset.add_metadata('gateampatinstr',gate_amplitude_value)

        with meas_aux_aux.run() as datasaver_aux_aux:
            varnames=[]
            #zurich.set_frequencies_to_json_config("160MHz_squeezed_singledot2")
            #print("JUST SET BACK FREQUENCIES")
            ############################################################################
            #datasaver.dataset.add_metadata('time_before_gvg',time_before_gvg)
            local_time=time.localtime()
            formatted_time = f"{local_time.tm_mday:02}:{local_time.tm_hour:02}:{local_time.tm_min:02}"
            print("day and time before spectrum (dd:hh:mm):", formatted_time)
            #datasaver.dataset.add_metadata('time_before_spectrum',time_before_spectrum)
            ############################################################################

            zurich.sigout1_amp1_enabled_param.value(0)
            returned_values_nodrive=take_long_spectra(reps=reps_nodrive,demod_ch=demod_ch)

            ###############################################################################
            local_time=time.localtime()
            formatted_time = f"{local_time.tm_mday:02}:{local_time.tm_hour:02}:{local_time.tm_min:02}"
            print("day and time after spectrum (dd:hh:mm):", formatted_time)
           # datasaver.dataset.add_metadata('time_after_spectrum',time_after_spectrum)
            ##############################################################################
            
            

            #read vslues needed for saving anc calculation
            avg_psd_array_nodrive=returned_values_nodrive['avg_psd']
            compressed_freq_array=returned_values_nodrive["compressed_freq"]
            meas_times_nodrive=returned_values_nodrive['meas_times']
            compressed_freq_array_real= returned_values_nodrive["compressed_freq_real"]
            filter=returned_values_nodrive["filter"]


            #now save these values
            for m_time,avg_psd in zip(meas_times_nodrive,avg_psd_array_nodrive):
                datasaver.add_result(#('Voltage_fft_avg', returned_values_nodrive['Voltage_fft_avg']),
                                        ('avg_psd', avg_psd),
                                        (freq_param,compressed_freq_array_real),
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
            print(f"max pos :{max_relative_freq}")
           
            """
            X, Y = np.meshgrid(compressed_freq_array, meas_times_nodrive, indexing='ij')
           # plt.ion()
            plt.pcolor(X,Y,avg_psd_array_nodrive.T)
            plt.title("driven psd vs time")
            plt.show()#now plot for testing purposes
            plt.pause(0.001)

            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.plot(max_relative_freq,1.1*max(avg_avg_psd_nodrive),'g*')
            plt.title("nondriven psd avg and positon of maximum")
            plt.show()
            plt.pause(0.001)

            

            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.title("nondriven avg psd")
            plt.show()
            plt.pause(0.001)

            
            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.title("narrowband driven avg psd together with non-driven psd")
            plt.plot(max_relative_freq,driven_value_narrowband,'g*')
            plt.show()
            plt.pause(0.001)
            """


            
            mask = (compressed_freq_array >= -30e3) & (compressed_freq_array <= 30e3)
            compressed_freq_array=compressed_freq_array[mask]
            avg_avg_psd_nodrive=avg_avg_psd_nodrive[mask]
            #avg_avg_psd_nodrive_filtercomp=avg_avg_psd_nodrive/filter
            #now fit lorentzian to scaled value
            Gamma_guess=5e3
            offset_approx=1.5e-15
            initial_guess=[max_relative_freq,Gamma_guess,max(avg_avg_psd_nodrive),min(avg_avg_psd_nodrive)]
            freq_span=max(compressed_freq_array)-min(compressed_freq_array)
            # Define bounds for the parameters
            #lower_bounds = [min(compressed_freq_array)+0.25*freq_span, 0, max(avg_avg_psd_nodrive/1.5, 0]  # Replace with appropriate lower bounds
            #upper_bounds = [max(compressed_freq_array)-0.25*freq_span, 1e3, np.inf, np.inf]  # Gamma is bounded to <= 1e3

            popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, compressed_freq_array, avg_avg_psd_nodrive, p0=initial_guess)
            
            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.title("Lorentzian fit initial guess")
            plt.plot(compressed_freq_array,lorentzian_fkt(compressed_freq_array,initial_guess[0],initial_guess[1],initial_guess[2],initial_guess[3]))
            plt.show()
            
            
            
            
            lorentzian, area_under_lorentzian=lorentzian_fkt_w_area(compressed_freq_array,popt[0],popt[1],popt[2],popt[3])

            #avg_avg_psd_nodrive_with_driven_value=copy.copy(avg_avg_psd_nodrive)
            #closest_index = np.argmin(np.abs(compressed_freq_array - max_relative_freq))
            #avg_avg_psd_nodrive_with_driven_value[closest_index]=driven_value_narrowband

            datasaver_aux_aux.add_result(('avg_avg_psd_nodrive',avg_avg_psd_nodrive),
             #                           ('avg_avg_psd_nodrive_w_driven_value',avg_avg_psd_nodrive_with_driven_value),
              #                          ('avg_avg_psd_drive',avg_avg_driven_psd[mask]),
               #                         ('avg_avg_psd_nodrive_scaled',avg_avg_driven_psd[mask]/drive_difference_narrowband),
                                        (freq_param,compressed_freq_array_real[mask]))
            

            #zurich.sigout1_amp1_enabled_param.value(0)
            #print("drive off")

            #print(f"driven_value_narrowband {driven_value_narrowband}")
            #print(f"drive_difference_narrowband {drive_difference_narrowband}")
            print(f"max(avg_avg_psd) {max(avg_avg_psd_nodrive)}")
            print(f"area_under_lorentzian {area_under_lorentzian}")
            #print(f"area_under_lorentzian scaled by drive: {area_under_lorentzian/drive_difference_narrowband}")
            #print(f"area_under_lorentzian scaled by slope: {area_under_lorentzian/slope^2}")
            print(f"width of lorentzian {popt[1]}")
            print(f"slope:{slope}")

            #datasaver.dataset.add_metadata('driven_value_narrowband',driven_value_narrowband)
            #datasaver.dataset.add_metadata('drive_difference_narrowband',drive_difference_narrowband)
            datasaver.dataset.add_metadata('max_avg_avg_psd_',max(avg_avg_psd_nodrive))
            datasaver.dataset.add_metadata('area_under_lorentzian',area_under_lorentzian)
            #datasaver.dataset.add_metadata('area_under_lorentzian_scaled_by_drive',area_under_lorentzian/drive_difference_narrowband)
            #datasaver.dataset.add_metadata('area_under_lorentzian_scaled_by_slope',area_under_lorentzian/slope^2)
            datasaver.dataset.add_metadata('freq_mech_corrected',freq_mech())
            datasaver.dataset.add_metadata('freq_rf_',freq_rf_value)
            datasaver.dataset.add_metadata('width_of_lorentzian',popt[1])
            datasaver.dataset.add_metadata('slope',slope)

            foldername='C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_'
            if not os.path.exists(foldername):
                os.makedirs(foldername) 
            run_id = datasaver.run_id

            filename=f'meas{run_id}_thermal_lo_fit.png'
            path = os.path.join(foldername, filename)

            plt.plot(compressed_freq_array,avg_avg_psd_nodrive)
            plt.title("Lorentzian fit")
            plt.plot(compressed_freq_array,lorentzian_fkt(compressed_freq_array,popt[0],popt[1],popt[2],popt[3]))##change back if not working!
            plt.savefig(path)
            plt.close()
            


#now import instance of exp_class

#from experiments.cs_experiment import CSExperiment

#temp_meas_180mK=CSExperiment()\
        

        temp_meas_fluctuating_base_mK2.area_values_scaled_by_slope.append(area_under_lorentzian/slope)
        temp_meas_fluctuating_base_mK2.area_values_unscaled.append(area_under_lorentzian)
        temp_meas_fluctuating_base_mK2.slopes.append(slope)


#for n in range(5):
   # do_GVg_and_adjust_sitpos()
run_thermomech_temp_meas()
#print(f"done round {n}")
#time.sleep(20)
    #print(n)

#do_GVg_and_adjust_sitpos()
#run_thermomech_temp_meas(reps_nodrive=100)
    #if n==10:
        #gate_amplitude_param(gate_amplitude_value/2)
     #   print("resetting gate amp")
    #print(f"overall rep nr {n}")