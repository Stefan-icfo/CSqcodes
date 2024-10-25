from qcodes import load_or_create_experiment, Measurement
from instruments import station, bilt
#from database import backupDatabase
#from utils.read_temperature import read_T
from utils.bot import adaptive_card
import os
import time
from tqdm import tqdm
import numpy as np

def create_loops(define_loops, define_channels):
    fast_param_name = define_loops['fast_loop']['channel']
    fast_param = define_channels[fast_param_name]['channel']
    fast_axis = fast_param.sweep(**define_loops['fast_loop']['params'])

    slow_param_name = define_loops['slow_loop']['channel']
    slow_param = define_channels[slow_param_name]['channel']
    slow_axis = slow_param.sweep(**define_loops['slow_loop']['params'])

    return fast_axis, slow_axis


def init_experiment(define_loops, define_channels):
    
    for channel_key in define_channels:
        # set all labels
        define_channels[channel_key]['channel'].label = define_channels[channel_key]['label']
        # set all ramp_rates
        define_channels[channel_key]['channel'].instrument.ramp_slope(define_channels[channel_key]['safe_ramp'])
        # put everythong into ramp mode
        define_channels[channel_key]['channel'].instrument.output_mode('ramp')

        # set all software threshold parameters
        define_channels[channel_key]['channel'].instrument.software_theshold_upper_current(.015)
        define_channels[channel_key]['channel'].instrument.software_theshold_lower_current(-.015)
        define_channels[channel_key]['channel'].instrument.software_theshold_delay(1000)
        define_channels[channel_key]['channel'].instrument.software_threshold_state('on')

    initialize_voltages_and_block(define_loops, define_channels)

def initialize_voltages_and_block(define_loops, define_channels):

    loop_param_list = [define_loops[key]['channel'] for key in define_loops]
    
    #set constant voltage
    for channel_key in define_channels:
        # constant voltage
        if channel_key not in loop_param_list:
            define_channels[channel_key]['channel'].set(define_channels[channel_key]['constant_voltage'])

    # set init sweep voltages
    for loop_key in define_loops:
        chan_name = define_loops[loop_key]['channel']
        define_channels[chan_name]['channel'].set(define_loops[loop_key]['params']['start'])

    #wait to set
    try:
        for channel_key in define_channels:
            define_channels[channel_key]['channel'].instrument.block_until_set()
    except KeyboardInterrupt:
        print('channels stopped')
    finally:
        [define_channels[key]['channel'](define_channels[key]['channel']()) for key in define_channels]


def bild_title(define_loops, bild_title_dict):
    out_string = ""
    out_string += f"({define_loops['slow_loop']['channel']},{define_loops['fast_loop']['channel']})"
    for key in bild_title_dict:
        out_string += key + ":" + f"{bild_title_dict[key]():.2E}" + f"{bild_title_dict[key].unit},"

    return out_string[:-1]


def general_2D_map(measured_parameter, define_loops, define_channels, TC_ZI, bilt_settling_time, bild_title_dict, snake=True, SAMPLE_NAME=""):

    init_experiment(define_loops, define_channels)
    fast_axis, slow_axis = create_loops(define_loops, define_channels)

    #defining the exp and measure params
    exp = load_or_create_experiment(experiment_name=os.path.basename(__file__)[:-3], sample_name= SAMPLE_NAME + bild_title(define_loops, bild_title_dict))
    exp2 = load_or_create_experiment(experiment_name=os.path.basename(__file__)[:-3], sample_name= "current measurments")
    meas = Measurement(exp=exp, station=station)
    meas2 = Measurement(exp=exp2, station=station)


    meas.register_parameter(fast_axis.parameter,paramtype="array")  # register the first independent parameter
    meas.register_parameter(slow_axis.parameter,paramtype="array")  # register the first independent parameter


    meas.register_parameter(measured_parameter, setpoints=(fast_axis.parameter, slow_axis.parameter), paramtype="array")  # now register the dependent one

    meas2.register_parameter(slow_axis.parameter)
    for chan in bilt.channels:
        meas2.register_parameter(chan.i, setpoints=(slow_axis.parameter,))  # now register the dependent one

    meas.add_after_run(lambda: (time.sleep(0.1), [chan.v(chan.v()) for chan in bilt.channels]), args=[])
    #meas.add_after_run(backupDatabase, args=[])


    with meas.run() as datasaver,  meas2.run() as  datasaver2:

        meas.add_after_run(adaptive_card, args=[datasaver, [read_T()]])

        temp_fast_vec = [None] * len(fast_axis)  #np.empty(len(fast_axis),dtype=np.complex_)
        
        helper_string_slow = f"slow -> {define_loops['slow_loop']['channel']}({define_loops['slow_loop']['params']['start']}->{define_loops['slow_loop']['params']['stop']})"
        helper_string_fast = f"fast -> {define_loops['fast_loop']['channel']}({define_loops['fast_loop']['params']['start']}->{define_loops['fast_loop']['params']['stop']})"
        
        for slow_value in tqdm(slow_axis, desc=helper_string_slow, dynamic_ncols=True):
            with slow_axis.parameter.instrument.output_mode.set_to('exp'):
                currents_measurements = [chan.i() for chan in bilt.channels] # Measures the current drawn for all open channels)
                slow_axis.set(slow_value)
                slow_axis.parameter.instrument.block_until_set() #to make sure the slow axis is really set before beginning the fast axis.

            with fast_axis.parameter.instrument.output_mode.set_to('exp'):
                for i, fast_value in enumerate(tqdm(fast_axis, desc=helper_string_fast, leave=False, dynamic_ncols=True)):

                    #set
                    fast_axis.set(fast_value)
                    # fast_axis.parameter.instrument.block_until_set()

                    #get
                    time.sleep(bilt_settling_time)
                    time.sleep(3*TC_ZI)

                    get_v = measured_parameter()

                    temp_fast_vec[i] = get_v


            #to deal with snake
            temp_fast_vec = temp_fast_vec if fast_axis[0]<fast_axis[-1] else np.flip(temp_fast_vec)
            temp_fast_axis_list = list(fast_axis) if fast_axis[0]<fast_axis[-1] else list(fast_axis)[::-1]


            #save
            datasaver.add_result((fast_axis.parameter, temp_fast_axis_list),
                                (slow_axis.parameter, slow_value),
                                (measured_parameter, temp_fast_vec))
            
            curr_list = [(chan.i, currents_measurements[p]) for p, chan in enumerate(bilt.channels)]
            datasaver2.add_result((slow_axis.parameter, slow_value), *curr_list)
            
            if snake:
                fast_axis.reverse()
            else:
                fast_axis.parameter.instrument.set_voltage_and_block(fast_axis[0])

        dataset = datasaver.dataset
    
    return dataset