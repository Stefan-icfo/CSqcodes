#syncing qdac and zurich via triggered aquisition

import numpy as np
from instruments import station, zurich
from qcodes.dataset import Measurement, new_experiment
from utils.CS_utils import *
import time
from tqdm import tqdm
from qcodes import Parameter
from utils.zurich_data_fkt import *
#from utils.CS_utils import centered_moving_average




exp_name="sync_qdac_zurich"
#exp_name="crosscap120MHz_g2_13Hz_1mV@instr50mK"
#device_name = 'CD11_D7_C1'
device_name = 'test'
demod_ch=0

BURST_DURATION = 1
SAMPLING_RATE = 13730#27470#
nr_burst=5

from instruments import qdac
channel_6 = qdac.channel(6)

# Explicitly set to sweep mode
channel_6.dc_mode("sweep") 

sweep_time_s = (1 - 0) / 0.01  # Total time in seconds for 10 mV/s from 0V to 1V
sweep = channel_6.dc_sweep(start_V=1, stop_V=0, points=100, repetitions=1, dwell_s=sweep_time_s / 100, stepped=False)

# Allocate internal trigger
internal_trigger = qdac.allocate_trigger()

# Start the sweep on trigger
sweep.start_on(internal_trigger)

# Configure external trigger output
external_trigger = qdac.external_triggers[3]
external_trigger.source_from_trigger(internal_trigger)
external_trigger.width_s(0.1)  # Set pulse width to 100µs (adjust as needed)
#external_trigger.amplitude(1)     # Set amplitude of external trigger pulse to 1V
qdac.trigger(internal_trigger)
freq_rlc = zurich.oscs.oscs0.freq


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
    device.demods[demod_ch].sample.y
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
# Set DAQ module for triggered acquisition
daq_module.type(1)  # Triggered acquisition mode
# Specify trigger source
#daq_module.set('trigger.source', 0)  # Set to 0 for auxiliary input trigger (external)
# Set trigger level (e.g., threshold in volts)
#daq_module.set('trigger.level', 0.5)  # Trigger level at 0.5 V
# Set trigger slope or edge type (0 for rising edge, 1 for falling edge, 2 for both)
#daq_module.set('trigger.edge', 0)  # Rising edge trigger
# Set trigger hysteresis to avoid noise-triggered events (e.g., in volts)
#daq_module.set('trigger.holdoff.time', 0.001)  # 1 ms hold-off time

    #   Start data acquisition
daq_module.execute()
time.sleep(2)  # Allow some time for the system to warm up

# # -----------------Start the Measurement-----------------------
time_offset=0
start_time=time.time()
with meas.run() as datasaver:

        timestamp0=np.nan
        for burst_idx in range(nr_burst):
            try:
            # Read data from the DAQ
                daq_data = daq_module.read(raw=False, clk_rate=clockbase)
                qdac.trigger(internal_trigger)
                #daq_module.forcetrigger(1)
                #print(f" switched trigger in burst {burst_idx + 1},current time: {trig_time}")
                

                saveandaddtime=True
                for node in sample_nodes:
                    if node in daq_data.keys():
                        print("node in keys")
                        for sig_burst in daq_data[node]:
                            value = sig_burst.value[0, :]  # Get the value of the signal
                            t = (sig_burst.time) 
                            #  print(f"t0={t0_burst}")
                            if value.size > 0:  # Check if value is not empty
                            # Apply averaging every 100 points
                                if node==device.demods[demod_ch].sample.x:
                                    x_data=value
                                    print("reading x")
                                if node==device.demods[demod_ch].sample.y:
                                    y_data=value
                                    print("reading y")
                           
                    else:
                        print(f"Burst {burst_idx + 1}: No data available for node {node}")
                        saveandaddtime=False
                xy_complex = x_data + 1j * y_data
                v_r = np.absolute(xy_complex)
                v_r_avg=centered_moving_average(v_r,n=10)
                theta = np.angle(xy_complex)
          
                
                if saveandaddtime:
                    current_time=time.time()-start_time
                    print(f" saving burst {burst_idx + 1},current time: {current_time}")
                    datasaver.add_result(('x', x_data),
                                ('y', y_data),
                                ('v_r', v_r),  
                               (time_param,t)        
                               )
                
            except Exception as e:
                print(f"Error reading data for burst {burst_idx + 1}: {e}")
            
            
            
            time.sleep(BURST_DURATION)
            current_time=time.time()-start_time
            
            #time.sleep(BURST_DURATION)
            time_offset+=BURST_DURATION
            print(time_offset)
     

 

