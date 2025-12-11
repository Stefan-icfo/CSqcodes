import numpy as np
import time
from instruments import station, qdac, Triton, zurich, exp
from experiments.zurich_experiments.spectrum_0925B import run_thermomech_temp_meas

# Set oscillator 0 frequency (using the driver shortcut)
#osc0_target_freq = 136.74e6
#zurich.freq0(osc0_target_freq)
#print(f"Set oscillator 0 to {osc0_target_freq / 1e6:.5f} MHz")
mixdown_f=zurich.freq1()
saved_original_bias=qdac.ch07.dc_constant_V()
bias = 0
bias_max=300e-6
step=20e-6
step_index = 0

while bias <= bias_max:
    print(f"\n==============================")
    print(f" Step {step_index}: Measuring at amplitude = {bias * 1e6:.2f} uV")
    print(f"==============================\n")

    # Set amplitude using your custom Zurich driver
    qdac.ch07.dc_constant_V(bias)

    # Optional: short delay
    time.sleep(0.5)
    exp.sit_at_max_Isens(side="right")
    zurich.set_mixdown(mixdown_f)
    time.sleep(100)
    # Run measurement
    run_thermomech_temp_meas(exp_name=f'sl_amp_sweep{bias*1e6:.3g} uV',reps_nodrive=30,background_id=472)

    bias +=step
    step_index += 1
qdac.ch07.dc_constant_V(saved_original_bias)