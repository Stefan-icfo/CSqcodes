from instruments import bilt, station, keithly
from experiment_functions.general_functions import general1D


def check_leakage(gate_params, other_params):
    for gate in bilt.gates:
        
        #sweep object
        gate_sweep_object = gate.v.sweep(**gate_params)
        
        #set sweep rate
        bilt.set_rates([gate_sweep_object], sweep_rate=other_params['sweep_rate'])
        
        #init constant gates and source
        bilt.set_constant_gates([gate_sweep_object], 0)

        #init sweep gate
        gate_sweep_object.set(gate_sweep_object[0])


        #wait until everything is there
        bilt.block_until_set()


        #make the sweep and measure current
        general1D(gate_sweep_object,
                  keithly.curr,
                  f'leakage-{gate.name}',
                  sleep_time=bilt.SETTLING_TIME,
                  function=gate.block_until_set)


if __name__ == '__main__':

    gate_params = {
        'start': 0,
        'stop': 1,
        'step': 0.1
    }
    other_params={
        'sweep_rate': 10e-2
    }

    check_leakage(gate_params)
