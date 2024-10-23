import time
from instruments import k2400, triton2
from experiment_functions.general_functions import general_time_1D
from qcodes.instrument.specialized_parameters import ElapsedTimeParameter

#------user input----------------
SAMPLE_NAME = 'test'
sleep = 1
volt_range = 200e-3
current_range = 10e-6
current_out = 5e-6
compliance_voltage = 1
triton2.T4() # to refresh the temperature

# Fast axis
time_axis = ElapsedTimeParameter('time_axis')

# Parameters definition
measured_parameter1 = k2400.volt

measured_parameter2 = triton2.T4

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
    measurement_parameters=[measured_parameter1,measured_parameter2], 
    sample_name=SAMPLE_NAME,
    sleep_time=sleep)