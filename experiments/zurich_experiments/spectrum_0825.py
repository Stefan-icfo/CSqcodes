#takes zurich demod fft in bursts
#plots and saves psd


import numpy as np
import scipy as scp
import os


from instruments import station, zurich, Triton, qdac, exp
from qcodes.dataset import Measurement, new_experiment

from utils.CS_utils import *

import time
from tqdm import tqdm
from qcodes import Parameter
import copy

from utils.zurich_data_fkt import *



Temp=0.035
time.sleep(10) 
device_name = 'CD12_B5_F4'
#exp_name=f"1dot_nodrive_spectrum_temp={Temp:4g}_zurichrange_divide_freq_by_half_nomask"#_cs_at_{sweet_CS_spot}
exp_name=f"Spectrum_{Temp:4g}_140mKg4dot"
from experiments.cs_experiment import *


filter_bw=200e3

###########################values for 219k data transfer######################

rbw=3.353#3.353#1.676#0.83819#3.353#2756972058123#0.808190#209.584e-3
BURST_DURATION = 298.262e-3##298.262e-3#298.262e-3#596.523e-3#1.193#1.1943 0.569523# 4.772 2.386#
SAMPLING_RATE = 219.72656250e3#13.73e3#219.72656250e3#219.72656250e3##109.86328125e3#54.93e3#109.86328125e3
#SAMPLING_RATE=13730

###########################values for 109k data transfer######################


nr_bursts=7
#reps=4
reps_nodrive=5
#reps_drive=20
demod_ch=3
drive_offset=0
#mode_freq=552.03e6
mask_boundary=100e3
avg_num=21



freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq

freq_rlc(1.25e6)

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
            full_data, averaged_data_per_burst, averaged_data, freq,compressed_freq,filter_data  = take_spectrum(demod_ch=demod_ch,SAMPLING_RATE = SAMPLING_RATE,BURST_DURATION=BURST_DURATION,nr_burst=nr_bursts, avg_num=avg_num)  
            #freq_real=freq+freq_mech()
            compressed_freq_real=compressed_freq+freq_mech()

            for data,avg_data in zip(full_data,averaged_data_per_burst):
                datas.append(data)
                avg_datas.append(avg_data)
                avg_data_psd=voltage_to_psd(avg_data, rbw)#calculate psd
                avg_datas_psd.append(avg_data_psd)
                meas_times.append(meas_time)

                
                meas_time+=BURST_DURATION*nr_bursts

            
    values_to_return={'Voltage_fft': np.array(datas),'Voltage_fft_avg' : np.array(avg_datas), 'avg_psd' : np.array(avg_datas_psd), "freq": np.array(freq), "compressed_freq" : np.array(compressed_freq), "meas_times" : np.array(meas_times),"compressed_freq_real" : np.array(compressed_freq_real),'filter':np.array(filter_data)}
    
    return values_to_return

#vars_to_save=[gate_ramp_slope,tc,vsd_dB,source_amplitude_instrumentlevel_GVg,vsdac,x_avg,y_avg]


gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes1.value
gate_amplitude_value = gate_amplitude_param()

def run_thermomech_temp_meas(reps_nodrive=reps_nodrive,exp_name=exp_name,fit_lorentzian=False,take_time_resolved_spectrum=False):
#    zurich.set_mixdown(mode_freq)

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


    
    #gate_amp_uV=gate_amplitude_param()*1e6
    
    # ----------------Create a measurement-------------------------
    experiment = new_experiment(name=exp_name+f"cs_at_{round(qdac.ch06.dc_constant_V(),4)}_outputsource={round(zurich.output0_amp0(),4)}", sample_name=device_name)
    meas = Measurement(exp=experiment)
    meas.register_parameter(time_param)  
    meas.register_parameter(freq_param) 
   # meas.register_custom_parameter('V_fft_avg', 'V_fft_avg', unit='V', basis=[], setpoints=[time_param,freq_param])
    # meas.register_parameter(measured_parameter, setpoints=[vgdc_sweep.parameter])  # register the 1st dependent parameter
    meas.register_custom_parameter('avg_psd', 'avg_psd', unit='W/Hz', basis=[], setpoints=[time_param,freq_param])



    experiment_1D = new_experiment(name=exp_name+'_1D'+f"cs_at_{round(qdac.ch06.dc_constant_V(),4)}_outputsource={round(zurich.output0_amp0(),4)}", sample_name=device_name)
    meas_aux_aux = Measurement(exp=experiment_1D)
    #meas_aux.register_parameter(time_param)  
    meas_aux_aux.register_parameter(freq_param)
    meas_aux_aux.register_custom_parameter('avg_avg_psd_nodrive', 'avg_avg_psd_nodrive', unit='W/Hz', basis=[], setpoints=[freq_param])
    meas_aux_aux.register_custom_parameter('V_fft_avg_avg', 'V_fft_avg_avg', unit='V', basis=[], setpoints=[freq_param])
    
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
            
            returned_values_nodrive=take_long_spectra(reps=reps_nodrive,demod_ch=demod_ch)

            ###############################################################################
            #read vslues needed for saving anc calculation
            avg_psd_array_nodrive=returned_values_nodrive['avg_psd']
            avg_V_array_nodrive=returned_values_nodrive['Voltage_fft_avg']
            compressed_freq_array=returned_values_nodrive["compressed_freq"]
            meas_times_nodrive=returned_values_nodrive['meas_times']
            compressed_freq_array_real= returned_values_nodrive["compressed_freq_real"]
            filter=returned_values_nodrive["filter"]


            #now save these values
            if take_time_resolved_spectrum:
                for m_time,avg_psd in zip(meas_times_nodrive,avg_psd_array_nodrive):
                    datasaver.add_result(#('V_fft_avg', avg_V_array_nodrive),
                                        ('avg_psd', avg_psd),
                                        (freq_param,compressed_freq_array_real),
                                        (time_param,m_time))
            #time the nodrive measurement ended, in burst time -needed for next datasave
            #switch_time=meas_times_nodrive[-1]

            
            avg_avg_psd_nodrive=np.mean(avg_psd_array_nodrive,axis=0)
            avg_avg_V_nodrive=np.mean(returned_values_nodrive['Voltage_fft_avg'],axis=0)
            
            #now calculate peak frequency
            #max_relative_freq = freq[np.argmax(averaged_data)]
            max_relative_freq=compressed_freq_array[np.argmax(avg_avg_psd_nodrive)]#offset from zero frequency of demodulator
            exp.max_thermomech_freq=compressed_freq_array_real[np.argmax(avg_avg_psd_nodrive)]
            #empiric:
            #max_relative_freq=max_relative_freq/2

            freq_rlc_value=freq_rlc()#mixdown frequency
            freq_rf_value=freq_rf()#source drive frequency, mech frequency minus (convention) mixdown frequency
            print(f"max pos :{max_relative_freq}")


            
            #mask = (compressed_freq_array >= -mask_boundary) & (compressed_freq_array <= mask_boundary)
            compressed_freq_array=compressed_freq_array#[mask]
            avg_avg_psd_nodrive=avg_avg_psd_nodrive#[mask]

            
            datasaver_aux_aux.add_result(('avg_avg_psd_nodrive',avg_avg_psd_nodrive),
                                        ('V_fft_avg_avg',avg_avg_V_nodrive),
              #                          ('avg_avg_psd_drive',avg_avg_driven_psd[mask]),
               #                         ('avg_avg_psd_nodrive_scaled',avg_avg_driven_psd[mask]/drive_difference_narrowband),
                                        #(freq_param,compressed_freq_array_real[mask]))
                                        (freq_param,compressed_freq_array_real))
            
            datasaver.dataset.add_metadata('max_avg_avg_psd_',max(avg_avg_psd_nodrive))
            datasaver.dataset.add_metadata('freq_mech_corrected',freq_mech())
            datasaver.dataset.add_metadata('freq_rf_',freq_rf_value)
            print(f"max(avg_avg_psd) {max(avg_avg_psd_nodrive)}")

            if fit_lorentzian:
                Gamma_guess=2e3
                offset_approx=1.5e-15
                initial_guess=[max_relative_freq,Gamma_guess,max(avg_avg_psd_nodrive),min(avg_avg_psd_nodrive)]
                freq_span=max(compressed_freq_array)-min(compressed_freq_array)
                
                try:
                    popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, compressed_freq_array, avg_avg_psd_nodrive, p0=initial_guess)
                except:
                    popt=initial_guess
                    print("fitting error")
          
                lorentzian, area_under_lorentzian=lorentzian_fkt_w_area(compressed_freq_array,popt[0],popt[1],popt[2],popt[3])


                #zurich.sigout1_amp1_enabled_param.value(0)
                #print("drive off")
                
                #print(f"driven_value_narrowband {driven_value_narrowband}")
                #print(f"drive_difference_narrowband {drive_difference_narrowband}")
                
                print(f"area_under_lorentzian {area_under_lorentzian}")
                #print(f"area_under_lorentzian scaled by drive: {area_under_lorentzian/drive_difference_narrowband}")
                #print(f"area_under_lorentzian scaled by slope: {area_under_lorentzian/slope^2}")
                print(f"width of lorentzian {popt[1]}")
                #print(f"slope:{slope}")

                #datasaver.dataset.add_metadata('driven_value_narrowband',driven_value_narrowband)
                #datasaver.dataset.add_metadata('drive_difference_narrowband',drive_difference_narrowband)
                
                datasaver.dataset.add_metadata('area_under_lorentzian',area_under_lorentzian)
                #datasaver.dataset.add_metadata('area_under_lorentzian_scaled_by_drive',area_under_lorentzian/drive_difference_narrowband)
                #datasaver.dataset.add_metadata('area_under_lorentzian_scaled_by_slope',area_under_lorentzian/slope^2)
                
                datasaver.dataset.add_metadata('width_of_lorentzian',popt[1])
                #sdatasaver.dataset.add_metadata('slope',slope)
                #sdatasaver.dataset.add_metadata('Temperature',Triton.MC())

                foldername='C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4'
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
            

run_thermomech_temp_meas(take_time_resolved_spectrum=False)
   





