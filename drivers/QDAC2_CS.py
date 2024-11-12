from drivers.QDAC2 import QDac2
import numpy as np
import time
from tqdm import tqdm
from qcodes.utils import validators

max_dc_V=3
min_dc_V=-3


class QDac2Channel(InstrumentChannel):

    def __init__(self, parent: 'QDac2', name: str, channum: int):
        super().__init__(parent, name)
        self._channum = channum

        # Modify dc_constant_V to use the new _set_limited_voltage function
        self.add_parameter(
            name='dc_constant_V',
            label=f'ch{channum}',
            unit='V',
            set_cmd=self._set_limited_voltage,  # Overwrite with _set_limited_voltage
            get_cmd=f'sour{channum}:volt?',
            get_parser=float,
            vals=validators.Numbers(min_dc_V, max_dc_V)
        )

    def _set_limited_voltage(self, value):
        """
        Custom method to set voltage with additional constraints.
        """
        # Retrieve the current voltage using the ask_channel method
        current_voltage = float(self.ask_channel('sour{0}:volt?'))

        # Check if the change is within the allowed 100 mV limit
        max_change = 0.1  # 100 mV
        if abs(value - current_voltage) > max_change:
            raise ValueError(f"Attempted to change voltage by more than {max_change * 1000:.1f} mV. "
                             f"Current voltage: {current_voltage} V, Requested voltage: {value} V")


        # Set the voltage directly if within limits
        self.write(f'sour{self._channum}:volt {value}')

    def ramp_ch(self, target_voltage: float):
        """
        Sets the voltage with a delay to allow the instrument to ramp at the specified slew rate.

        Args:
            target_voltage (float): The desired voltage to reach.
            slew_rate (float): The slew rate in V/s (default 0.01 V/s).
        
        Raises:
            ValueError: If the target voltage is out of range.
        """
        current_voltage = float(self.dc_constant_V())  # Get current voltage
        slew_rate =  float(self.dc_dc_slew_rate_V_per_s())
        max_voltage = 3.0
        min_voltage = -3.0

        if not (min_voltage <= target_voltage <= max_voltage):
            raise ValueError(f"Target voltage {target_voltage} V is out of bounds. Allowed range is {min_voltage}V to {max_voltage}V.")

        # Calculate the total ramping time based on the slew rate
        voltage_difference = abs(target_voltage - current_voltage)
        wait_time = voltage_difference / slew_rate  # Time required for ramping

        # Set the target voltage
        self._set_limited_voltage(target_voltage)

        # Wait for the ramping to complete
        time.sleep(wait_time)


class QDac2_CS(QDac2):
    def __init__(self, name: str, address: str, **kwargs):
        super().__init__(name, address, **kwargs)
        # Additional initialization here if needed

    def ramp_multi_ch_slowly(self, channels, final_vgs, step_size: float = 10e-3, ramp_speed: float = 1e-3):
        # Calculate the wait time per step based on ramp speed
        wait_time = step_size / ramp_speed
        
        # Get initial voltages for each channel and calculate the required steps
        start_points = [self.channel(ch).dc_constant_V() for ch in channels]
        step_counts = [int(abs(final_vgs[i] - start_points[i]) / step_size) for i in range(len(channels))]
        max_steps = max(step_counts)

        # Generate voltage sweeps for each channel to reach the final voltages
        V_sweeps = [
            np.linspace(start_points[i], final_vgs[i], num=max_steps) for i in range(len(channels))
        ]

        # Apply each step of the voltage ramp to all channels
        for step in tqdm(range(max_steps)):
            for j, ch in enumerate(channels):
                self.channel(ch).dc_constant_V(V_sweeps[j][step])
            time.sleep(wait_time)

    def ramp_multi_ch_fast(self, qdac_channels, final_vgs):
        
        # Get initial voltages for each channel and calculate the required steps
        start_points = [qdac_channels.dc_constant_V() for ch in qdac_channels]
        slew_rates   = [qdac_channels.dc_dc_slew_rate_V_per_s() for ch in qdac_channels]
        wait_times   = [(end_point-start_point)/slew_rate for end_point,start_point,slew_rate in zip(start_points,final_vgs,slew_rates)]
        wait_time=max(wait_times)
        for ch,final_vg in zip(qdac_channels,final_vgs):
            ch.dc_constant_V(final_vg)
        
        
        time.sleep(wait_time)

    def read_channels(self, chan_nr: int = 7):
        """
        Reads the voltage, slew rate, and filter settings of the specified number of channels
        and prints them in a formatted output.
        
        Args:
            chan_nr (int): The number of channels to read (default is 7).
        """
        # Collect information about each channel in a loop
        for i in range(1, chan_nr + 1):
            ch_voltage = self.channel(i).dc_constant_V()  # Voltage of channel
            ch_slew_rate = self.channel(i).dc_slew_rate_V_per_s()  # Slew rate of channel
            ch_filter = self.channel(i).output_filter()  # Filter setting of channel
            
            # Print formatted information for each channel
            print(f"QDac channel {i}: Voltage = {ch_voltage:.3f} V, "
                  f"Slew Rate = {ch_slew_rate:.3f} V/s, "
                  f"Filter = {ch_filter}")
    
    def set_all_slewrates(self, slew_rate: float = 0.01):
        """
        Sets the slew rate for all channels to the specified value.

        Args:
            slew_rate (float): The slew rate in V/s to set for each channel.
        """
        for ch in range(1, 8):  # Assuming you want to set this for channels 1 through 7
            self.channel(ch).dc_slew_rate_V_per_s(slew_rate)