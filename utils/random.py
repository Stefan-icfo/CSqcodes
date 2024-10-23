import numpy as np

def calculate_measurement_time_2D_map(rangeslow, rangefast, number_of_points_fast, number_of_points_slow, ramp_rate, integration_time, snake, waiting_time):
    sum=0
    
    #moving part
    snake_coef = 1 if snake else 2
    slow_time = rangeslow/ramp_rate/1000
    fast_time = rangefast/ramp_rate/1000


    sum += slow_time + number_of_points_slow*snake_coef*fast_time

    #integration part
    sum += number_of_points_fast*number_of_points_slow*(integration_time+waiting_time)
    

    return sum

def V_array(start, end, v_step):
    if start > end:
        step = int((start-end)/v_step)
    else:
        step = int((end-start)/v_step)

    V_array = np.linspace(start, end, step)

    return V_array

def check_mode(mode, chans_fixed, chans_sweep):
    
    if len(chans_fixed)+len(chans_sweep) != 5:
        raise ValueError('You have chosen an incorrect number of channnels. Total number channels must be 5.')

    if mode == 'single_gate':
        if len(chans_sweep) != 1:
            raise ValueError('For single gate sweep mode the total number of channels_sweep must be 1.')

    if mode == 'all_gates':
        if len(chans_sweep) != 5:
            raise ValueError('For all gates sweep mode the total number of channels_sweep must be 5.')

    if mode == 'custom_gates':
        pass

    else:
        raise ValueError('The mode you chose does not exist. The mode argument can only take: "single_gate", "custom_gate" or "all_gates"')

def barrier_gen(slow_value,fast_value):
    barrier_value = slow_value/2 + fast_value/2
    
    return (barrier_value,)