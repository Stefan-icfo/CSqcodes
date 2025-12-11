import numpy as np
import time
from instruments import station, qdac, Triton, zurich, exp
from experiments.zurich_experiments.spectrum_0925B import run_thermomech_temp_meas

# Set oscillator 0 frequency (using the driver shortcut)
#osc0_target_freq = 136.74e6
#zurich.freq0(osc0_target_freq)
#print(f"Set oscillator 0 to {osc0_target_freq / 1e6:.5f} MHz")

# Start amplitude sweep from 2 mV to 75 mV with +20% steps
amplitude_V = 5e-3
amplitude_max_V = 75e-3
scale_factor = 1.20
step_index = 0
mixdown_f=zurich.freq1()
while amplitude_V <= amplitude_max_V:
    print(f"\n==============================")
    print(f" Step {step_index}: Measuring at amplitude = {amplitude_V * 1e3:.2f} mV")
    print(f"==============================\n")
    exp.set_params(source_amplitude_instrumentlevel_GVg=amplitude_V)
    exp.sit_at_max_Isens(side="right")
    # Set amplitude using your custom Zurich driver
    zurich.output0_amp0(amplitude_V)
    zurich.set_mixdown(mixdown_f)
    time.sleep(100)

    # Optional: short delay
    time.sleep(0.5)

    # Run measurement
    run_thermomech_temp_meas(exp_name=f'sl_amp_sweep{amplitude_V*1e3:.3g}mV',reps_nodrive=50,background_id=354)

    amplitude_V *= scale_factor
    step_index += 1
exp.set_params(source_amplitude_instrumentlevel_GVg=20e-3)