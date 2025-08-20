import numpy as np
import time
from instruments import station, qdac, Triton, zurich, exp
from experiments.zurich_experiments.spectrum_0825 import run_thermomech_temp_meas

# Set oscillator 0 frequency (using the driver shortcut)
#osc0_target_freq = 136.74e6
#zurich.freq0(osc0_target_freq)
#print(f"Set oscillator 0 to {osc0_target_freq / 1e6:.5f} MHz")

# Start amplitude sweep from 2 mV to 75 mV with +20% steps
amplitude_V = 2e-3
amplitude_max_V = 75e-3
scale_factor = 1.20
step_index = 0

while amplitude_V <= amplitude_max_V:
    print(f"\n==============================")
    print(f" Step {step_index}: Measuring at amplitude = {amplitude_V * 1e3:.2f} mV")
    print(f"==============================\n")

    # Set amplitude using your custom Zurich driver
    zurich.output0_amp0(amplitude_V)

    # Optional: short delay
    time.sleep(0.5)

    # Run measurement
    run_thermomech_temp_meas()

    amplitude_V *= scale_factor
    step_index += 1
