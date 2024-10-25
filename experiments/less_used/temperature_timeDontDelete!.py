import time
from instruments import k2400, Triton
from experiment_functions.roger1D_time import general_time_1D
from qcodes.instrument.specialized_parameters import ElapsedTimeParameter

#------user input----------------
SAMPLE_NAME = 'Actual cool down data 2'
sleep = 150 # time it takes to loop over all the thermometers + 5 seconds
volt_range = 200e-3
current_range = 10e-6
current_out = 10e-6
compliance_voltage = 1

# Fast axis
time_axis = ElapsedTimeParameter('time_axis')

# Parameters definition
measured_parameter1 = k2400.volt
measured_parameter2 = k2400.T_cernox
measured_parameter3 = Triton.T5

# ================ Initialise Keithley2400 ========================

k2400.mode('CURR')
k2400.sense('VOLT')
k2400.write('SYST:RSEN ON') # enables 4point measurement
k2400.compliancev(compliance_voltage)

k2400.write(f'SENS:VOLT:RANG {volt_range}')
k2400.write(f'SOUR:CURR:RANG {current_range}')
k2400.curr(current_out)

k2400.output(True)

# ======================== Title ==================================

general_time_1D(
    time_axis=time_axis, 
    measurement_parameters=[measured_parameter1,
                            measured_parameter2,
                            measured_parameter3], 
    sample_name=SAMPLE_NAME,
    sleep_time=sleep)