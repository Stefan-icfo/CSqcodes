#syncing qdac and zurich via triggered aquisition

import numpy as np
from instruments import station, zurich
from qcodes.dataset import Measurement, new_experiment
from utils.CS_utils import *
import time
from tqdm import tqdm
from qcodes import Parameter
from utils.zurich_data_fkt import *





exp_name="sync_qdac_zurich_execute_v2"
#exp_name="crosscap120MHz_g2_13Hz_1mV@instr50mK"
#device_name = 'CD11_D7_C1'
device_name = 'test1'
demod_ch=0

BURST_DURATION = 1
SAMPLING_RATE = 13730#27470#
nr_burst=2

from instruments import qdac

internal_trigger = qdac.allocate_trigger()
INT1 = internal_trigger

external_trigger = qdac.external_triggers[3]
external_trigger.source_from_trigger(internal_trigger)
external_trigger.width_s(0.1)  

qdac.trigger(internal_trigger)
freq_rlc = zurich.oscs.oscs0.freq




time_param = Parameter('time_param',
                            label='time',
                            unit='s',  
                            set_cmd=None,  
                            get_cmd=None)  



experiment= new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(time_param)  
#meas_aux.register_parameter(freq_param)
meas.register_custom_parameter('x', 'x', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('y', 'y', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('v_r', 'v_r', unit='V', basis=[], setpoints=[time_param])
meas.register_custom_parameter('Phase', 'Phase', unit='V', basis=[], setpoints=[time_param])


session = Session("localhost")
device = session.connect_device("DEV20039")  
device.demods[demod_ch].enable(True)
print("Connected to the device successfully.")

daq_module = session.modules.daq
daq_module.device(device)

daq_module.grid.mode(2)
time.sleep(2)

sample_nodes = [
    device.demods[demod_ch].sample.x,
    device.demods[demod_ch].sample.y
    ]
# '/dev20039/scopes/0/wave'

TOTAL_DURATION = BURST_DURATION * nr_burst
num_cols = int(np.ceil(SAMPLING_RATE * BURST_DURATION))

daq_module.count(nr_burst)  
daq_module.duration(BURST_DURATION)  
daq_module.grid.cols(num_cols)  
daq_module.type(0)
    
for node in sample_nodes:
    daq_module.subscribe(node)
   
    

clockbase = device.clockbase()


    #   Start data acquisition
daq_module.execute()
time.sleep(2)  # Allow some time for the system to warm up
qdac.write('sour8:sine:freq 1') 
#set voltage to 250 mv  
qdac.write('sour8:sine:span 0.5') 
qdac.write('sour8:sine:count inf')
#qdac.write('sour8:sine:trig:sour imm')
qdac.write('sour8:sine:trig:sour INT1')
qdac.write('sour8:sine:init:imm')

# trigger 
#qdac.write('sour8:sine:trig:sour INT1')
#qdac.write('sour8:sine:init') 
#qdac.write('sour8:sine:abor')
#
time_offset=0
start_time=time.time()
with meas.run() as datasaver:

        timestamp0=np.nan
        for burst_idx in range(nr_burst):
            try:
            # Read data from the DAQ
                
                time.sleep(0.01)
                qdac.write('tint 1')
                #qqdac.trigger(internal_trigger)
                daq_data = daq_module.read(raw=False, clk_rate=clockbase)
                #qdac.trigger(internal_trigger)
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
                    local_t = t
                    for ti, xi, yi, vri, phii in zip(local_t, x_data, y_data, v_r, theta):
                        datasaver.add_result(
                            (time_param, ti),
                            ('x',         xi),
                            ('y',         yi),
                            ('v_r',       vri),
                            ('Phase',     phii),
                        )
                
            except Exception as e:
                print(f"Error reading data for burst {burst_idx + 1}: {e}")
            
            
            
            time.sleep(BURST_DURATION)
            current_time=time.time()-start_time
            
            #time.sleep(BURST_DURATION)
            time_offset+=BURST_DURATION
            print(time_offset)
     


