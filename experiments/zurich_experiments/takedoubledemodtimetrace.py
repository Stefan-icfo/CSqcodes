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


exp_name="autocorrelation_short_noiseforRuggi"

#exp_name="autocorrelation_20s_150mK_onICT"
#exp_name="crosscap120MHz_g2_13Hz_1mV@instr50mK"
device_name = 'CD11_D7_C1'

demod_ch1=3
demod_ch2=4


filter_bw=100
#rbw=13
#rbw=200e-3


#BURST_DURATION = (on_time+off_time)/bursts_per_cycle
BURST_DURATION = 1
SAMPLING_RATE = 13730#13730#27470#109900
nr_burst=6

#on_times=[4,8,12,16]
#off_times=[6,10,14,18]
nr_avg=41


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
meas.register_custom_parameter('x1', 'x1', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('y1', 'y1', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('v_r1', 'v_r1', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('Phase1', 'Phase1', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('v_r1_avg', 'v_r1_avg', unit='V', basis=[], setpoints=[time_param])

meas.register_custom_parameter('x2', 'x2', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('y2', 'y2', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('v_r2', 'v_r2', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('Phase2', 'Phase2', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('v_r2_avg', 'v_r2_avg', unit='V', basis=[], setpoints=[time_param])

meas.register_custom_parameter('code_time', 'code_time', unit='s', basis=[], setpoints=[time_param])
meas.register_custom_parameter('time_since_burst_start','time_since_burst_start', unit='s', basis=[], setpoints=[time_param])

# meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)

#connecting to device
session = Session("localhost")
device = session.connect_device("DEV20039")  # Replace with your actual device ID
device.demods[demod_ch1].enable(True)
device.demods[demod_ch2].enable(True)
print("Connected to the device successfully.")

daq_module = session.modules.daq
daq_module.device(device)
daq_module.type(1)  # triggered acquisition
daq_module.grid.mode(2)
time.sleep(2)

sample_nodes = [
    device.demods[demod_ch1].sample.x,
    device.demods[demod_ch1].sample.y,
    device.demods[demod_ch2].sample.x,
    device.demods[demod_ch2].sample.y
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
#gate_rf_enabled_param = getattr(zurich.sigouts.sigouts1.enables, f'enables{1}')
#gate_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes0.value
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
        #datasaver.dataset.add_metadata('gate_amp_at_instr',gate_amplitude_param())

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
              
                current_time=time.time()-start_time
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
                                if node==device.demods[demod_ch1].sample.x:
                                    x1_data=value
                                    print("reading x1")
                                if node==device.demods[demod_ch1].sample.y:
                                    y1_data=value
                                    print("reading y1")     
                                if node==device.demods[demod_ch2].sample.x:
                                    x2_data=value
                                    print("reading x2")
                                if node==device.demods[demod_ch2].sample.y:
                                    y2_data=value
                                    print("reading y2")        
                    else:
                        print(f"Burst {burst_idx + 1}: No data available for node {node}")
                        saveandaddtime=False
                xy1_complex = x1_data + 1j * y1_data
                xy2_complex = x2_data + 1j * y2_data
                v_r1 = np.absolute(xy1_complex)
                v_r2 = np.absolute(xy2_complex)
                v_r1_avg=centered_moving_average(v_r1,n=nr_avg)
                v_r2_avg=centered_moving_average(v_r2,n=nr_avg)
                theta1 = np.angle(xy1_complex)
                theta2 = np.angle(xy2_complex)
          
                
                if saveandaddtime:
                    current_time=time.time()-start_time
                    print(f" saving burst {burst_idx + 1},current time: {current_time}")
                    datasaver.add_result(('x1', x1_data),
                                         ('x2', x2_data),
                                ('y1', y1_data),
                                ('y2', y2_data),
                                ('v_r1', v_r1),
                                ('v_r2', v_r2),
                                ('v_r1_avg', v_r1_avg),
                                ('v_r2_avg', v_r2_avg),
                                ('Phase1', theta1),
                                ('Phase2', theta2),
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
     

 



