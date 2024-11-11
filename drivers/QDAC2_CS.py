from drivers.QDAC2 import QDac2
import numpy as np
import time

class QDac2_CS(QDac2):
    def __init__(self, name: str, address: str, **kwargs):
        super().__init__(name, address, **kwargs)
        # Additional initialization here if needed

    def ramp_QDAC_multi_ch_slowly(self, channels, final_vgs, step_size: float = 10e-3, ramp_speed: float = 1e-3):
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
        for step in range(max_steps):
            for j, ch in enumerate(channels):
                self.channel(ch).dc_constant_V(V_sweeps[j][step])
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