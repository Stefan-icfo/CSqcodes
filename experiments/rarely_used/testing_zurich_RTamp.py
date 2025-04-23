from instruments import  zurich
import time
from qcodes import Parameter
from qcodes.dataset import Measurement, new_experiment
from tqdm import tqdm

tc=1e-2
vsdac=1 #dummy
measured_parameter1 = zurich.demod0
measured_parameter2 = zurich.demod1
source1_amplitude_param = zurich.output0_amp0
source2_amplitude_param = zurich.sigouts.sigouts1.amplitudes.amplitudes0.value

start_value=1e-7
length=90
instr_power_sweep=[0]+[start_value * (1.1 ** i) for i in range(length)]
drive_mag_param = Parameter('drive_mag', label='drive_mag', unit='Vrms')
print(f"max value {max(instr_power_sweep)}" )

time.sleep(20)

exp_name='instrument_test'
device_name='_test_slow_20dBattn_div40_switchchannels'
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(drive_mag_param)
meas.register_custom_parameter('V_1', unit='V', basis=[], setpoints=[drive_mag_param])
meas.register_custom_parameter('V_2', unit='V', basis=[], setpoints=[drive_mag_param])
with meas.run() as datasaver:
    for power_value in tqdm(instr_power_sweep, leave=False, desc='instr_power_sweep', colour='green'):
        source1_amplitude_param(power_value)
        source2_amplitude_param(power_value)
        time.sleep(3*tc+1)
        _, v_r_calc1, _, _ = zurich.phase_voltage_current_conductance_compensate(vsdac,x_avg=0,y_avg=0)
        time.sleep(3*tc)
        _, v_r_calc2, _, _ = zurich.phase_voltage_current_conductance_compensate(vsdac,measured_value=zurich.demod1(),x_avg=0,y_avg=0)

        datasaver.add_result(('V_1', v_r_calc1/40), ('V_2', v_r_calc2), (drive_mag_param, power_value))

source1_amplitude_param(0)
source2_amplitude_param(0)