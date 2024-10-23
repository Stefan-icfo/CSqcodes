from instruments import station
import os
#from utils.read_temperature import read_T
from qcodes import Measurement, load_or_create_experiment
import time
from drivers.bilt import SweepMultiParam
from utils.tqdm_local import tqdm_local
import numpy as np
from utils.random import barrier_gen
from database import backupDatabase
from utils.bot import adaptive_card


def general1D(sweep_object, measurement_parameters, sample_name, sleep_time=0, function=None):
    exp = load_or_create_experiment(experiment_name=os.path.basename(__file__)[
        :-3], sample_name=sample_name)
    meas = Measurement(exp=exp, station=station)

    meas.register_parameter(sweep_object.parameter, paramtype="array")

    for par in measurement_parameters:
        meas.register_parameter(par, setpoints=(
            sweep_object.parameter,), paramtype="array")

    meas.add_after_run(backupDatabase, args=[])

    with meas.run() as datasaver:

        meas.add_after_run(adaptive_card, args=[
                           datasaver, [read_T()]])

        for value in tqdm_local(sweep_object):

            sweep_object.set(value)

            if function:
                function()

            time.sleep(sleep_time)

            # this will construct list of tupples made for add_results
            measured_values = [(par, par()) for par in measurement_parameters]

            if isinstance(sweep_object, SweepMultiParam):
                value = value[0]

            datasaver.add_result(
                (sweep_object.parameter, value),
                *measured_values
            )

#------------------------------------------------------------------------------------------------------------------------------------------

def general_2D(fast_axis, slow_axis, measured_parameters, sleep, sample_name, slow_function=None, fast_function=None, snake=True):

    exp = load_or_create_experiment(experiment_name=os.path.basename(__file__)[
                                    :-3], sample_name=sample_name)
    meas = Measurement(exp=exp, station=station)

    # register the first independent parameter
    meas.register_parameter(fast_axis.parameter, paramtype="array")
    # register the first independent parameter
    meas.register_parameter(slow_axis.parameter, paramtype="array")

    for par in measured_parameters:
        meas.register_parameter(par, setpoints=(
            fast_axis.parameter, slow_axis.parameter), paramtype="array")  # now register the dependent one

    meas.add_after_run(backupDatabase, args=[])

    with meas.run() as datasaver:

        meas.add_after_run(adaptive_card, args=[
                           datasaver, [read_T()]])

        temp_fast_vecs = [[None] * len(fast_axis) for par in measured_parameters]

        for slow_value in tqdm_local(slow_axis):
            slow_axis.set(slow_value)

            if slow_function:
                slow_function()

            for i, fast_value in enumerate(tqdm_local(fast_axis, is_fast_axis=True)):
                fast_axis.set(fast_value)

                if fast_function:
                    fast_function()

                time.sleep(sleep)

                for par_num, par in enumerate(measured_parameters):
                    temp_fast_vecs[par_num][i] = par()

            # to deal with snake
            temp_fast_axis_list = list(fast_axis)
            if fast_axis[0] > fast_axis[-1]:  # this is the reversed case
                temp_fast_vecs = [np.flip(meas) for meas in temp_fast_vecs]
                temp_fast_axis_list = temp_fast_axis_list[::-1]

            # to deal with the multi sweep object
            if isinstance(slow_axis, SweepMultiParam):
                slow_value = slow_value[0]
            if isinstance(fast_axis, SweepMultiParam):
                temp_fast_axis_list = [elem[0] for elem in temp_fast_axis_list]

            #this will construct list of tupples for add_result
            for_add_results = [(par, meas) for par, meas in zip(measured_parameters, temp_fast_vecs)]

            datasaver.add_result((fast_axis.parameter, temp_fast_axis_list),
                                 (slow_axis.parameter, slow_value),
                                 *for_add_results)

            if snake:
                fast_axis.reverse()


if __name__ == "__main__":
    from instruments import bilt, zurich
    from utils.titles import build_simple_title

    sweep_slow = bilt.gate3.v.sweep(0, 1, 0.1)
    sweep_fast = bilt.sweep_multi_channel([
        {
            'channel': bilt.gate1.v,
            'sweep_values': {
                'start': 0,
                'stop': 1,
                'step': 0.1,
            },
        },
        {
            'channel': bilt.gate2.v,
            'sweep_values': {
                'start': 1,
                'stop': 2,
                'step': 0.1,
            }
        }
    ])

    # sweep_slow = zurich.sweep_multi_output([
    #     {
    #         'channel': zurich.source.freq,
    #         'sweep_values': {
    #             'start': 1.5e6,
    #             'stop': 1.6e6,
    #             'step': 10e3,
    #         },
    #     },
    #     {
    #         'channel': zurich.drain.freq,
    #         'sweep_values': {
    #             'start': 1.6e6,
    #             'stop': 1.7e6,
    #             'step': 10e3,
    #         }
    #     }
    # ])

    bild_title_dict = {
        'temp': read_T(),
        'voltage': zurich.sigouts.sigouts0.amplitudes3,
        'gates': sweep_slow,
        'gate2': sweep_fast,
    }

    a = build_simple_title('slaven_', bild_title_dict)
    print(a)

    bilt.set_rates([sweep_fast], sweep_rate=1.2e-4)

    # general_2D(sweep_fast, sweep_slow, [zurich.demods.demods0.sample], 0.1, "whatever2d",
    #            slow_function=bilt.gates.block_until_set, fast_function=bilt.gates.block_until_set, snake=True)
    # general1D(sweep_slow, [zurich.demods.demods0.sample], 'whatever', 1)

#------------------------------------------------------------------------------------------------------------------------------------------


def general_time_2D(fast_axis, time_axis, measured_parameters, sleep, sample_name, fast_function=None, snake=True):

    exp = load_or_create_experiment(experiment_name=os.path.basename(__file__)[
                                    :-3], sample_name=sample_name)
    meas = Measurement(exp=exp, station=station)

    # register the first independent parameter
    meas.register_parameter(fast_axis.parameter, paramtype="array")
    # register the first independent parameter
    meas.register_parameter(time_axis, paramtype="array")

    for par in measured_parameters:
        meas.register_parameter(par, setpoints=(
            fast_axis.parameter, time_axis), paramtype="array")  # now register the dependent one

    meas.add_after_run(backupDatabase, args=[])

    with meas.run() as datasaver:

        time_axis.reset_clock()

        while True:

            meas.add_after_run(adaptive_card, args=[
                            datasaver, [read_T()]])

            temp_fast_vecs = [[None] * len(fast_axis) for par in measured_parameters]

            for i, fast_value in enumerate(tqdm_local(fast_axis, is_fast_axis=True)):
                fast_axis.set(fast_value)

                if fast_function:
                    fast_function()

                time.sleep(sleep)

                for par_num, par in enumerate(measured_parameters):
                    temp_fast_vecs[par_num][i] = par()

            # to deal with snake
            temp_fast_axis_list = list(fast_axis)
            if fast_axis[0] > fast_axis[-1]:  # this is the reversed case
                temp_fast_vecs = [np.flip(meas) for meas in temp_fast_vecs]
                temp_fast_axis_list = temp_fast_axis_list[::-1]

            # to deal with the multi sweep object
            if isinstance(time_axis, SweepMultiParam):
                slow_value = slow_value[0]
            if isinstance(fast_axis, SweepMultiParam):
                temp_fast_axis_list = [elem[0] for elem in temp_fast_axis_list]

            #this will construct list of tupples for add_result
            for_add_results = [(par, meas) for par, meas in zip(measured_parameters, temp_fast_vecs)]

            datasaver.add_result((fast_axis.parameter, temp_fast_axis_list),
                                (time_axis, time_axis()),
                                *for_add_results)

            if snake:
                fast_axis.reverse()

# ------------------------------------------------------------------------------------------------------------------------------------------

def general_symmetric_2D(fast_axis, slow_axis, barrier_axis, measured_parameters, sleep, sample_name, slow_function=None, fast_function=None, snake=True):

    exp = load_or_create_experiment(experiment_name=os.path.basename(__file__)[
                                    :-3], sample_name=sample_name)
    meas = Measurement(exp=exp, station=station)

    # register the first independent parameters
    meas.register_parameter(fast_axis.parameter, paramtype="array")
    meas.register_parameter(slow_axis.parameter, paramtype="array")
    meas.register_parameter(barrier_axis.parameter, paramtype="array") # <-------------------

    for par in measured_parameters:
        meas.register_parameter(par, setpoints=(
            fast_axis.parameter, slow_axis.parameter), paramtype="array")  # now register the dependent one

    meas.add_after_run(backupDatabase, args=[])

    with meas.run() as datasaver:

        meas.add_after_run(adaptive_card, args=[
                           datasaver, [read_T()]])

        temp_fast_vecs = [[None] * len(fast_axis) for par in measured_parameters]

        for slow_value in tqdm_local(slow_axis):
            slow_axis.set(slow_value)

            if slow_function:
                slow_function()

            for i, fast_value in enumerate(tqdm_local(fast_axis, is_fast_axis=True)):
                fast_axis.set(fast_value)
                barrier_value = barrier_gen(slow_value[0],fast_value[0])
                barrier_axis.set(barrier_value)

                if fast_function:
                    fast_function()

                time.sleep(sleep)

                for par_num, par in enumerate(measured_parameters):
                    temp_fast_vecs[par_num][i] = par()

            # to deal with snake
            temp_fast_axis_list = list(fast_axis)
            if fast_axis[0] > fast_axis[-1]:  # this is the reversed case
                temp_fast_vecs = [np.flip(meas) for meas in temp_fast_vecs]
                temp_fast_axis_list = temp_fast_axis_list[::-1]

            # to deal with the multi sweep object
            if isinstance(slow_axis, SweepMultiParam):
                slow_value = slow_value[0]
            if isinstance(fast_axis, SweepMultiParam):
                temp_fast_axis_list = [elem[0] for elem in temp_fast_axis_list]

            #this will construct list of tupples for add_result
            for_add_results = [(par, meas) for par, meas in zip(measured_parameters, temp_fast_vecs)]

            datasaver.add_result((fast_axis.parameter, temp_fast_axis_list),
                                 (slow_axis.parameter, slow_value),
                                 (barrier_axis.parameter, barrier_value),
                                 *for_add_results)

            if snake:
                fast_axis.reverse()