import os
import time
from instruments import station
from qcodes import Measurement, load_or_create_experiment

# -----------------------------------------------------------------------------------------
        
def general_time_1D(time_axis, measurement_parameters, sample_name, sleep_time=0, function=None):
    exp = load_or_create_experiment(experiment_name=os.path.basename(__file__)[
        :-3], sample_name=sample_name)

    time_axis.reset_clock() 

    meas = Measurement(exp=exp, station=station)

    meas.register_parameter(time_axis, paramtype="array")

    for par in measurement_parameters:
        meas.register_parameter(par, setpoints=(
            time_axis,), paramtype="array")


    with meas.run() as datasaver:
        print('Time trace running...')
        while True:
            time.sleep(sleep_time)
            
            # this will construct list of tupples made for add_results
            measured_values = [(par, par()) for par in measurement_parameters]

            datasaver.add_result(
                (time_axis,time_axis()),
                *measured_values
            )
