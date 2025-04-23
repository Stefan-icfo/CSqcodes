# The program does gate sweep with Zurich LIA at RLC frequency; can use multiple gates
# Stefan Forstner 
# conducance vs gate voltage using zurich


import numpy as np


from instruments import station, zurich, Triton
from qcodes.dataset import Measurement, new_experiment

from utils.CS_utils import centered_moving_average, save_metadata_var, get_var_name

import time
from tqdm import tqdm
from qcodes import Parameter

from utils.zurich_data_fkt import *
#from utils.CS_utils import centered_moving_average




exp_name=f"ringupringdown_20k_75mV_singledot_temp={Triton.MC()}"
#exp_name="crosscap120MHz_g2_13Hz_1mV@instr50mK"
#device_name = 'CD11_D7_C1'
#exp_name="test_"
device_name = 'CD11_D7_C1'
demod_ch=4
proxy_ch=5


filter_bw=20e3
#rbw=13
#rbw=200e-3


#BURST_DURATION = (on_time+off_time)/bursts_per_cycle
BURST_DURATION = 1
SAMPLING_RATE = 27470#54.93e3#27470#13730#27470#
nr_burst=5

#on_times=[4,8,12,16]
#off_times=[6,10,14,18]
nr_avg=21


freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq


#vars_to_save=[gate_ramp_slope,tc,vsd_dB,source_amplitude_instrumentlevel_GVg,vsdac,x_avg,y_avg]

time_param = Parameter('time_param',
                            label='time',
                            unit='s',  # Or the appropriate unit
                            set_cmd=None,  # If there is no instrument to communicate with directly
                            get_cmd=None)  # Define get_cmd if you need to read a value





# ----------------Create a measurement-------------------------





experiment= new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(time_param)  
#meas_aux.register_parameter(freq_param)
meas.register_custom_parameter('x', 'x', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('y', 'y', unit='V', basis=[], setpoints=[time_param])

meas.register_custom_parameter('v_r', 'v_r', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('Phase', 'Phase', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('v_r_avg', 'v_r_avg', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('code_time', 'code_time', unit='s', basis=[], setpoints=[time_param])
meas.register_custom_parameter('time_since_burst_start','time_since_burst_start', unit='s', basis=[], setpoints=[time_param])

meas.register_custom_parameter('x_proxy', 'x_proxy', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('y_proxy', 'y_proxy', unit='V', basis=[], setpoints=[time_param])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)

#connecting to device
session = Session("localhost")
device = session.connect_device("DEV20039")  # Replace with your actual device ID
device.demods[demod_ch].enable(True)
print("Connected to the device successfully.")

daq_module = session.modules.daq
daq_module.device(device)
daq_module.type(1)  # triggered acquisition
daq_module.grid.mode(2)
time.sleep(2)

sample_nodes = [
    device.demods[demod_ch].sample.x,
    device.demods[demod_ch].sample.y,
    device.demods[proxy_ch].sample.x,
    device.demods[proxy_ch].sample.y
    ]

TOTAL_DURATION = BURST_DURATION * nr_burst
num_cols = int(np.ceil(SAMPLING_RATE * BURST_DURATION))
#num_bursts = int(np.ceil(TOTAL_DURATION / BURST_DURATION))
daq_module.count(nr_burst)  # Set number of bursts to collect
daq_module.duration(BURST_DURATION)  # Set duration of each burst
daq_module.grid.cols(num_cols)  # Set number of columns

    # Subscribe to sample node
for node in sample_nodes:
    daq_module.subscribe(node)
    #daq_module.subscribe(sample_nodes[1])
    # Retrieve clock base for timing

clockbase = device.clockbase()

    #   Start data acquisition
daq_module.execute()
time.sleep(2)  # Allow some time for the system to warm up

# # -----------------Start the Measurement-----------------------
gate_rf_enabled_param = getattr(zurich.sigouts.sigouts1.enables, f'enables{1}')
gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes0.value
start_time=time.time()
with meas.run() as datasaver:
    #datasaver.dataset.add_metadata('qdac_ch01_dc_constant_V_start',qdac.ch01.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch02_dc_constant_V_start',qdac.ch02.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch03_dc_constant_V_start',qdac.ch03.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch04_dc_constant_V_start',qdac.ch04.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch05_dc_constant_V_start',qdac.ch05.dc_constant_V())
    #datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V_start',qdac.ch06.dc_constant_V())
        datasaver.dataset.add_metadata('probe_freq',freq_rf())
        datasaver.dataset.add_metadata('rlc_freq',freq_rlc())
        datasaver.dataset.add_metadata('center_freq',freq_mech())
        datasaver.dataset.add_metadata('filter_bw',filter_bw)
        datasaver.dataset.add_metadata('gate_amp_at_instr',gate_amplitude_param())

    #datasaver.dataset.add_metadata('rbw',rbw)
    #with meas_aux.run() as datasaver_aux:
        #varnames=[]
    #or i in range(len(vars_to_save)):
    #    varnames.append(get_var_name(vars_to_save[i]))
    #save_metadata_var(datasaver.dataset,varnames,vars_to_save)
    # for i in range(2):
        time_offset=0
        timestamp0=np.nan
        pre_trig_times=[]
        post_trig_times=[]
        gate_on_times=[]
        gate_off_times=[]
        for burst_idx in range(nr_burst):
            try:
            # Read data from the DAQ
                daq_data = daq_module.read(raw=False, clk_rate=clockbase)
                pre_trig_times.append(time.time()-start_time)
                daq_module.forcetrigger(1)
                trig_time=time.time()-start_time
                post_trig_times.append(trig_time)
                #print(f" switched trigger in burst {burst_idx + 1},current time: {trig_time}")
                
                if burst_idx % 2==0: 
                    
                    gate_rf_enabled_param.value(1)
                    current_time=time.time()-start_time
                    gate_on_times.append(current_time)
                    #triggerdelay=trig_time-current_time
                    #print(f" switched on gate  {burst_idx + 1},current time: {current_time}")
                else:
                    gate_rf_enabled_param.value(0)
                    current_time=time.time()-start_time
                    gate_off_times.append(current_time)
                    #triggerdelay=trig_time-current_time
                    #print(f" switched off gate in {burst_idx + 1},current time: {current_time}")

                
                
                saveandaddtime=True
                for node in sample_nodes:
                    if node in daq_data.keys():
                        print("node in keys")
                        for sig_burst in daq_data[node]:
                            value = sig_burst.value[0, :]  # Get the value of the signal
                            #t0_burst = sig_burst.header['createdtimestamp'][0] / clockbase
                            #if np.any(np.isnan(timestamp0)):
                                # Set our first timestamp to the first timestamp we obtain.
                                #timestamp0 = sig_burst["timestamp"][0, 0]
                            #t = (sig_burst["timestamp"][0, :] - timestamp0) / clockbase
                            
                            t = (sig_burst.time)+time_offset  
                            #  print(f"t0={t0_burst}")
                            if value.size > 0:  # Check if value is not empty
                            # Apply averaging every 100 points
                                if node==device.demods[demod_ch].sample.x:
                                    x_data=value
                                    print("reading x")
                                if node==device.demods[demod_ch].sample.y:
                                    y_data=value
                                    print("reading y")
                                if node==device.demods[proxy_ch].sample.x:
                                    x_proxy=value
                                    print("reading proxy x")
                                if node==device.demods[proxy_ch].sample.y:
                                    y_proxy=value
                                    print("reading proxy y")
                           
                    else:
                        print(f"Burst {burst_idx + 1}: No data available for node {node}")
                        saveandaddtime=False
                xy_complex = x_data + 1j * y_data
                v_r = np.absolute(xy_complex)
                v_r_avg=centered_moving_average(v_r,n=nr_avg)
                theta = np.angle(xy_complex)
          
                
                if saveandaddtime:
                    current_time=time.time()-start_time
                    print(f" saving burst {burst_idx + 1},current time: {current_time}")
                    datasaver.add_result(('x', x_data),
                                ('y', y_data),
                                ('x_proxy', x_proxy),
                                ('y_proxy', y_proxy),
                                ('v_r', v_r),
                                ('v_r_avg', v_r_avg),
                                ('Phase', theta),
                                ('time_since_burst_start',t-time_offset),
                               (time_param,t)        
                               )
                
            except Exception as e:
                print(f"Error reading data for burst {burst_idx + 1}: {e}")
            
            
            
            time.sleep(BURST_DURATION)
            current_time=time.time()-start_time
            
            #time.sleep(BURST_DURATION)
            time_offset+=BURST_DURATION
            print(time_offset)
        i=0
        for pre_trig_time in pre_trig_times:
            i+=1
            datasaver.dataset.add_metadata(f'pre_trig_time_{i}',pre_trig_time)

        i=0
        for post_trig_time in post_trig_times:
            i+=1
            datasaver.dataset.add_metadata(f'post_trig_time_{i}',post_trig_time)   
        
        i=0
        for gate_on_time in gate_on_times:
            i+=1
            datasaver.dataset.add_metadata(f'gate_on_time_{i}',gate_on_time) 

        i=0
        for gate_off_time in gate_off_times:
            i+=1
            datasaver.dataset.add_metadata(f'gate_off_time_{i}',gate_off_time)   

gate_rf_enabled_param.value(0)

            #full_time_data,full_x_data,full_y_data = demod_xy_timetrace(sample_nodes=sample_nodes, daq_module=daq_module, device=device, demod_ch=0)    
        #meas_time=0
        #for data_x,data_y,t in zip(full_x_data,full_y_data,full_time_data):
        #    datasaver.add_result(('x', data_x),
        #                        ('y', data_y),
        #                        (time_param,t))
                
                #target_size = np.shape(avg_data)[0]
               # factor = len(freq) // target_size  # Factor by which to compress

                # Reshape the array and compute the mean along the compressed axis
                

                #meas_time+=BURST_DURATION*nr_bursts



