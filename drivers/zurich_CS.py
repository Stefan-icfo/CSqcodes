from os import name
import zhinst.qcodes as ziqc
from qcodes.instrument_drivers.zurich_instruments.uhfli import UHFLI
from qcodes.instrument.parameter import DelegateParameter, ScaledParameter
from qcodes.instrument.channel import InstrumentChannel, ChannelList
import numpy as np
import time


class MyZurich(ziqc.UHFLI):
    '''
    [use_cache_conductance]

    This parameter (bool) of source.conductance
    that controls whether to use the last returned value of the 
    "source.voltage" function. e.g.

    measured_parameter = zurich.source.voltage()
    measured_parameter2 = zurich.source.conductance
    zurich.source.conductance.use_cache_conductance = True

    Will not read the voltage from the Zurich twice, whereas the 
    following will:
    
    measured_parameter = zurich.source.voltage()
    measured_parameter2 = zurich.source.conductance
    zurich.source.conductance.use_cache_conductance = False

    Note here that the order is important. Using the cached value 
    relies on it being the right one.

    '''


    def __init__(self, manual_instr, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.manual_instr = manual_instr

        output0_amp0=self.sigouts.sigouts0.amplitudes.amplitudes0.value
        output1_amp1=self.sigouts.sigouts1.amplitudes.amplitudes1.value
        demod0=self.demods.demods0.sample
        demod1=self.demods.demods1.sample
        demod3=self.demods.demods2.sample
        demod2=self.demods.demods3.sample




        '''
        Convention:

        Measurement always at frequency of osc0 on demod3.
        Source frequency set on osc1 and output on amplitudes4.
        Gate output always on amplitudes6 (modulation SB C-M)
        
        '''
    def phase_voltage_current_conductance_compensate(self, vsdac, x_avg=0, y_avg=0,measured_value=None, gain_RT=200, gain_HEMT=5.64, Z_tot=7521):
        """
        This function calculates the compensated phase, voltage, current, and conductance
        based on measured values and calibration parameters.

        Parameters:
            measured_value (dict): Contains 'x' and 'y' measurement components.
            vsdac (float): The AC voltage value used in the measurement.
            x_avg (float): The average value for x to compensate.
            y_avg (float): The average value for y to compensate.
            gain_RT (float): Room temperature gain, default is 200.
            gain_HEMT (float): HEMT gain, default is 5.64.
            Z_tot (float): Total impedance, default is 7521 ohms.

        Returns:
            tuple: Contains (theta, v_r, I, G) - phase angle, voltage, current, conductance.
        """
        if measured_value is None:
            measured_value = self.demods.demods0.sample
        # Compensate x and y with the provided averages
        x = measured_value['x'][0] - x_avg  # Compensated x
        y = measured_value['y'][0] - y_avg  # Compensated y

        # Calculate complex representation of compensated x and y
        xy_complex = complex(x, y)
        v_r = np.absolute(xy_complex)  # Voltage magnitude
        theta = np.angle(xy_complex)  # Phase angle

        # Calculate current (I) and conductance (G)
        I = v_r / (gain_RT * gain_HEMT * Z_tot)
        G = 1 / ((vsdac / I) - Z_tot)

        return theta, v_r, I, G
    
    def x_y_avg(self, measured_parameter, tc=100e-3, avg_nr=100):
        """
        This function calculates the average x and y values over a specified number of measurements.

        Parameters:
            measured_parameter (callable): A callable that returns a measurement containing 'x' and 'y'.
            tc (float): Time in seconds to wait between measurements, default is 100 ms.
            avg_nr (int): Number of measurements to average, default is 100.

        Returns:
            tuple: Contains (x_avg, y_avg) - the averaged x and y values.
        """
        x_sum = 0
        y_sum = 0
        for n in range(avg_nr):
            time.sleep(tc)  # Wait for specified time constant between measurements
            measured_value = measured_parameter()  # Get the current measurement
            x_sum += measured_value['x'][0]  # Add the x component to the sum
            y_sum += measured_value['y'][0]  # Add the y component to the sum
            
        # Calculate the average values for x and y
        x_avg = x_sum / avg_nr
        y_avg = y_sum / avg_nr

        return x_avg, y_avg
    
    def save_config_to_metadata(self, datasaver):
        """
        Adds the dc_constant_V metadata for each channel to the given datasaver.

        Args:
            datasaver: The datasaver object where metadata will be stored.
            prefix (str): Prefix for metadata keys (default is 'qdac').
        """
        keys,values=[],[]
        keys.append("self.sigouts.sigouts0.amplitudes.amplitudes0.value")
        values.append(self.sigouts.sigouts0.amplitudes.amplitudes0.value) 
        keys.append("self.sigouts.sigouts1.amplitudes.amplitudes1.value")
        values.append(self.sigouts.sigouts1.amplitudes.amplitudes1.value)
        for key,value in zip(keys,values):
            datasaver.dataset.add_metadata(key, value)