from drivers.QDAC2 import QDac2
import numpy as np
import time
import json
from tqdm import tqdm
from qcodes.utils import validators
from utils.CS_utils import get_gate_Vs_from_metadata



class QDac2_CS(QDac2):
    def __init__(self, name: str, address: str, **kwargs):
        super().__init__(name, address, **kwargs)
        # Additional initialization here if needed

    def ramp_multi_ch_slowly(self, channels, final_vgs, step_size: float = 10e-3, ramp_speed: float = 1e-3):
   
        # Calculate the wait time per step based on ramp speed
        wait_time = step_size / ramp_speed

        # Get initial voltages for each channel
        if all(isinstance(ch, int) for ch in channels):
            start_points = [self.channel(ch).dc_constant_V() for ch in channels]
        else:
            start_points = [ch.dc_constant_V() for ch in channels]

        # Calculate the required steps and ensure at least one step for very small changes
        step_counts = [
            max(1, int(abs(final_vgs[i] - start_points[i]) / step_size))
            for i in range(len(channels))
        ]
        max_steps = max(step_counts)

        # Generate voltage sweeps for each channel to reach the final voltages
        V_sweeps = [
            np.linspace(start_points[i], final_vgs[i], num=max_steps)
            for i in range(len(channels))
        ]

        # Apply each step of the voltage ramp to all channels
        for step in tqdm(range(max_steps)):
            if all(isinstance(ch, int) for ch in channels):
                for j, ch in enumerate(channels):
                    self.channel(ch).dc_constant_V(V_sweeps[j][step])
            else:
                for j, ch in enumerate(channels):
                    ch.dc_constant_V(V_sweeps[j][step])  # Handles channel objects
            time.sleep(wait_time)

        # Ensure the final voltages are set accurately
        if all(isinstance(ch, int) for ch in channels):
            for j, ch in enumerate(channels):
                self.channel(ch).dc_constant_V(final_vgs[j])
        else:
            for j, ch in enumerate(channels):
                ch.dc_constant_V(final_vgs[j])  # Handles channel objects



    def ramp_multi_ch_fast(self, qdac_channels, final_vgs):     
        # Get initial voltages for each channel and calculate the required steps
        start_points = [ch.dc_constant_V() for ch in qdac_channels]
        slew_rates   = [ch.dc_slew_rate_V_per_s() for ch in qdac_channels]
        wait_times   = [(end_point-start_point)/slew_rate for end_point,start_point,slew_rate in zip(start_points,final_vgs,slew_rates)]
        wait_time=max(wait_times)
        for ch,final_vg in zip(qdac_channels,final_vgs):
            ch.dc_constant_V(final_vg)
        
        
        time.sleep(abs(wait_time))

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
            print(f"QDac channel {i}: Voltage = {ch_voltage:.4g} V, "
                  f"Slew Rate = {ch_slew_rate:.2f} V/s, "
                  f"Filter = {ch_filter}")
    
    def set_all_slewrates(self, slew_rate: float = 0.01):
        """
        Sets the slew rate for all channels to the specified value.

        Args:
            slew_rate (float): The slew rate in V/s to set for each channel.
        """
        for ch in range(1, 8):  # Assuming you want to set this for channels 1 through 7
            self.channel(ch).dc_slew_rate_V_per_s(slew_rate)

    def add_dc_voltages_to_metadata(self, datasaver, ch_num = 7, prefix: str = 'qdac'):
        """
        Adds the dc_constant_V metadata for each channel to the given datasaver.

        Args:
            datasaver: The datasaver object where metadata will be stored.
            prefix (str): Prefix for metadata keys (default is 'qdac').
        """
        for i in range(1, ch_num+1):
            channel = getattr(self, f'ch{i:02}')
            voltage = channel.dc_constant_V()  # Retrieve the current voltage
            metadata_key = f"{prefix}_ch{i:02}_dc_constant_V"
            datasaver.dataset.add_metadata(metadata_key, voltage)

    def set_gates_to_metadata_config(self,meas_id,pre_str='qdac_ch0',post_str='_dc_constant_V',gate_nrs=[1,2,3,4,5,6]):
    
        gates_dict=get_gate_Vs_from_metadata(meas_id,pre_str=pre_str,post_str=post_str,gate_nrs=gate_nrs)
        gate_nrs,target_Vs=[],[]
        for key, value in gates_dict.items():
            gate_nrs.append(int(key[-1]))
            target_Vs.append(value)
        self.ramp_multi_ch_slowly(gate_nrs,target_Vs)

    def set_gates_to_json_config(self, unique_name, file_path="parameters.json",gate_nrs=None):
        """
        Set gates to the voltage configuration saved in a JSON file.

        Args:
            file_path (str): Path to the JSON file containing voltage configurations.
            unique_name (str): The unique name of the saved configuration to apply.
            gate_nrs (list, optional): List of gate numbers to consider; if None, all gates in the configuration are used.
        """
        try:
            # Load the JSON file
            with open(file_path, 'r') as file:
                data = json.load(file)

            # Retrieve the configuration for the unique name
            gates_dict = data.get(unique_name)
            if not gates_dict:
                print(f"No configuration found for '{unique_name}' in {file_path}")
                return

            # Filter gates based on the provided gate numbers, if any
            gate_nrs_list = []
            target_Vs = []
            for key, value in gates_dict.items():
                gate_nr = int(key.lstrip("qdac_ch"))  # Extract the gate number from the key
                if gate_nrs is None or gate_nr in gate_nrs:
                    gate_nrs_list.append(gate_nr)
                    target_Vs.append(value)

            # Apply the voltages to the gates
            print(gate_nrs_list)
            print(target_Vs)
            self.ramp_multi_ch_slowly(gate_nrs_list, target_Vs)
            print(f"Gates set to configuration '{unique_name}' from {file_path}")

        except (FileNotFoundError, json.JSONDecodeError):
            print(f"Error loading configurations from {file_path}")

    def save_current_config_to_json(self, unique_name,file_path="parameters.json"):
        """
        Save the current QDAC configuration to a JSON file under a unique name.

        Args:
            file_path (str): Path to the JSON file where the configuration will be saved.
            unique_name (str): The unique name for this configuration.
        """
        try:
            # Load the existing file or create a new dictionary if the file doesn't exist
            try:
                with open(file_path, 'r') as file:
                    data = json.load(file)
            except (FileNotFoundError, json.JSONDecodeError):
                data = {}

            # Get the current configuration (voltages on each channel)
            current_config = {}
            for ch in range(1, 8):  # Assuming channels are 1 through 7
                channel = getattr(self, f'ch{ch:02}')
                current_config[f"qdac_ch{ch}"] = channel.dc_constant_V()

            # Add or update the configuration under the unique name
            data[unique_name] = current_config

            # Save the updated data back to the file
            with open(file_path, 'w') as file:
                json.dump(data, file, indent=4)
            print(f"Current configuration saved as '{unique_name}' in {file_path}")

        except Exception as e:
            print(f"An error occurred while saving the configuration: {e}")