import numpy as np
import time
from instruments import station, qdac, Triton, zurich, exp


amplitude_V =75e-3
amplitude_min_V = 500e-6#0.1e-3
scale_factor = 0.9#<1!
#scale_substract=20e-6
#reps=10
#initial_reps=5
step_index = 0
exp.sit_at_max_Isens(side="left")

while amplitude_V >= amplitude_min_V:
    print(f"\n==============================")
    print(f" Step {step_index}: Measuring at amplitude = {amplitude_V * 1e3:.2f} mV")
    print(f"==============================\n")

    # Set amplitude using your custom Zurich driver
    zurich.output1_amp1(amplitude_V)

    # Optional: short delay
    time.sleep(0.5)

    # Run measurement
    
    
    #for n in range(reps):
    exp.mech_simple_fun_db(costum_prefix='g2_dot_g2_drive')
    first_step=False
    amplitude_V *= scale_factor
    #amplitude_V -= scale_substract
    step_index += 1
exp.sit_at_max_Isens(side="left")