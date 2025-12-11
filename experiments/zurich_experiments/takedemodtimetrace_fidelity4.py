# The program does gate sweep with Zurich LIA at RLC frequency; can use multiple gates
# Stefan Forstner 
# conducance vs gate voltage using zurich


import numpy as np


from instruments import zurich
from qcodes.dataset import Measurement, new_experiment

from utils.CS_utils import centered_moving_average, save_metadata_var, get_var_name

import time
from tqdm import tqdm
from qcodes import Parameter

from utils.zurich_data_fkt import *
#from utils.CS_utils import centered_moving_average
import os

# Read TC from environment (set by sweep notebook). None if not provided.
_tc_env = os.environ.get("TC_NOW_S", None)
TC_NOW = float(_tc_env) if _tc_env is not None else None


exp_name="demodtimetrace"
if TC_NOW is None:
    costum_prefix = "fidelitysignaltoSL20mv"
else:
    costum_prefix = f"LOOPSTART1.25MHZfidelitysignaltoSL20mvC{TC_NOW*1e6:.3f}us_"

#exp_name="autocorrelation_20s_150mK_onICT"
#exp_name="crosscap120MHz_g2_13Hz_1mV@instr50mK"
device_name = 'CD12_B5_F4'

demod_ch=3

filter_bw=1e3




#BURST_DURATION = (on_time+off_time)/bursts_per_cycle

SAMPLING_RATE = 1.75781250e6#13.73291016e3#109.86328125e3#219.72656250e3#27470#109900


def takedemodtimetrace(SAMPLING_RATE=SAMPLING_RATE,
                       filter_bw=filter_bw,
                       costum_prefix=costum_prefix):


    freq_rf = zurich.oscs.oscs0.freq


#vars_to_save=[gate_ramp_slope,tc,vsd_dB,source_amplitude_instrumentlevel_GVg,vsdac,x_avg,y_avg]

    time_param = Parameter('time_param',
                            label='time',
                            unit='s',  # Or the appropriate unit
                            set_cmd=None,  # If there is no instrument to communicate with directly
                            get_cmd=None)  # Define get_cmd if you need to read a value





# ----------------Create a measurement-------------------------



    experiment= new_experiment(name=costum_prefix+exp_name, sample_name=device_name)
    meas = Measurement(exp=experiment)
    meas.register_parameter(time_param)  
    meas.register_custom_parameter('x', 'x', unit='V', basis=[], setpoints=[time_param])
    meas.register_custom_parameter('y', 'y', unit='V', basis=[], setpoints=[time_param])
    meas.register_custom_parameter('v_r', 'v_r', unit='V', basis=[], setpoints=[time_param])
    meas.register_custom_parameter('Phase', 'Phase', unit='V', basis=[], setpoints=[time_param])
    meas.register_custom_parameter('code_time', 'code_time', unit='s', basis=[], setpoints=[time_param])

    # meas.add_after_run(end_game, args = [instr_dict]) # Runs the line after the run is finished, even if the code stops abruptly :)

    #connecting to device
    session = Session("localhost")
    device = session.connect_device("DEV20039")  # Replace with your actual device ID
    device.demods[demod_ch].enable(True)
    print("Connected to the device successfully.")

    daq_module = session.modules.daq
    daq_module.device(device)
        
    daq_module.clearhistory() # reset history
    daq_module.type(0) # 0: continuous acquisition mode, 6: trigger on line source
    daq_module.holdoff.time(10e-3) #time before next acquisition
    daq_module.delay(0) # offset the signal recorded with respect to the trigger
    #daq_module.triggernode('/dev20039/demods/0/sample.TrigIn1') # trigger source, comment if no trigger

    daq_module.grid.mode(4)
    time.sleep(2)

    sample_nodes = [
    device.demods[demod_ch].sample.x,
    device.demods[demod_ch].sample.y
    ]

    num_cols = 100000#int(np.ceil(SAMPLING_RATE * BURST_DURATION))
    num_rep = 1

  #  daq_module.count(nr_burst)  # Set number of bursts to collect
 #   daq_module.duration(BURST_DURATION)  # Set duration of each burst
    daq_module.grid.cols(num_cols)  # Set number of columns
    daq_module.grid.rows(num_rep)
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

    start_time=time.time()
    with meas.run() as datasaver:
        
        datasaver.dataset.add_metadata('probe_freq',freq_rf())
        datasaver.dataset.add_metadata('filter_bw',filter_bw)
        time_offset=0
        timestamp0=np.nan
        pre_trig_times=[]
        post_trig_times=[]
    
        # try:
        # Read data from the DAQ
        daq_data = daq_module.read(raw=False, clk_rate=clockbase)
        x_data,y_data= [],[]
        for k in range(num_rep):
            for node in sample_nodes:
                for sig_burst in daq_data[node]:
                    
                    value = sig_burst.value[k, :]  # Get the value of the signal
                    if value.size > 0:  # Check if value is not empty
                        if node==device.demods[demod_ch].sample.x:
                            x_data=value
                            print("reading x")
                        if node==device.demods[demod_ch].sample.y:
                            y_data=value
                            print("reading y")           
        x_data = np.array(x_data)  
        y_data =np.array(y_data)
        xy_complex = x_data + 1j * y_data
        v_r = np.absolute(xy_complex)
        #  v_r_avg=centered_moving_average(v_r,n=nr_avg)
        #theta = np.angle(xy_complex)
    
        t = sig_burst.time
        
        

        datasaver.add_result(('x', x_data),
                        ('y', y_data),
                        ('v_r', v_r),
                        #('v_r_avg', v_r_avg),
                        #('Phase', theta),
                        #('time_since_burst_start',t-time_offset), 
                        (time_param,t)        
                        )
        
        #except Exception as e:
           # print(f"Error reading data ")
        
        
    daq_module.finish()      # Stop execution
    daq_module.clearhistory()  # optional
    daq_module.unsubscribe("*")   # Remove all subscriptions   
    

 

takedemodtimetrace()

