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
frequency_start = 137.95e6
frequency_stop = 138.05e6
frequency_increment = 5e3
frequency=copy.copy(frequency_start)



while frequency <= frequency_stop:
    print(f"\n==============================")
    print(f" Measuring at frequency = {frequency * 1e3:.2f} mV")
    print(f"==============================\n")

    # Set amplitude using your custom Zurich driver
    zurich.freq1(frequency)

    # Optional: short delay
    time.sleep(0.5)

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
    time.sleep(0.5)

    # Run measurement
    run_thermomech_temp_meas(exp_name=f"thermal spectrum-mix_and_weakly_driven_at_{frequency}")

    frequency+=frequency_increment


    