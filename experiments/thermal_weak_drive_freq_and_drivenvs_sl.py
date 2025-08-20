import numpy as np
import time
from instruments import station, qdac, Triton, zurich, exp
from experiments.zurich_experiments.spectrum_0825 import run_thermomech_temp_meas
import copy

# Set oscillator 0 frequency (using the driver shortcut)
#osc0_target_freq = 136.74e6
#zurich.freq0(osc0_target_freq)
#print(f"Set oscillator 0 to {osc0_target_freq / 1e6:.5f} MHz")

# Start amplitude sweep from 2 mV to 75 mV with +20% steps
zurich.set_mixdown(137987000,side="+")
frequency_start = 137.95e6
frequency_stop = 138.05e6#137.99e6
frequency_increment = 5e3
frequency=copy.copy(frequency_start)
frequency=137.99e6


while frequency <= frequency_stop:
    print(f"\n==============================")
    print(f" Measuring at frequency = {frequency * 1e3:.2f} mV")
    print(f"==============================\n")

    # Set amplitude using your custom Zurich driver
    zurich.freq1(frequency)

    # Optional: short delay
    time.sleep(200)#wait

    # Run measurement
    run_thermomech_temp_meas(exp_name=f"thermal spectrum+mix_and_weakly_driven_at_{frequency}")

    frequency+=frequency_increment

#going to mix-

zurich.set_mixdown(137987000)
frequency=copy.copy(frequency_start)   
while frequency <= frequency_stop:
    print(f"\n==============================")
    print(f" Measuring at frequency = {frequency * 1e3:.2f} mV")
    print(f"==============================\n")

    # Set amplitude using your custom Zurich driver
    zurich.freq1(frequency)

    # Optional: short delay
    time.sleep(200)#wait

    # Run measurement
    run_thermomech_temp_meas(exp_name=f"thermal spectrum-mix_and_weakly_driven_at_{frequency}")

    frequency+=frequency_increment

zurich.set_mixdown(137987000)
#driven_vs_sl
amplitude_V = 2e-3
drive_V=100e-6
amplitude_max_V = 75e-3
scale_factor = 1.5
step_index = 0

zurich.output1_amp1(drive_V)

while amplitude_V <= amplitude_max_V:
    print(f"\n==============================")
    print(f" Step {step_index}: Measuring at amplitude = {amplitude_V * 1e3:.2f} mV")
    print(f"==============================\n")

    # Set amplitude using your custom Zurich driver
    zurich.output0_amp0(amplitude_V)

    # Optional: short delay
    time.sleep(200)#wait

    # Run measurement
    run_thermomech_temp_meas(exp_name="thermal_and_drive_vs_sl")

    amplitude_V *= scale_factor
    step_index += 1
  