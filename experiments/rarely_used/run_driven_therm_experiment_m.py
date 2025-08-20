import numpy as np
import time
from experiments.zurich_experiments.spectrum_0825 import run_thermomech_temp_meas
from instruments import zurich

# Parameters
start_drive = 50e-6       # 50 μV
end_drive = 5e-3          # 5 mV
step_drive = 100e-6       # 100 μV
Temp = 0.035              # Temperature in K

# Loop over drive amplitudes
drive = start_drive
while drive <= end_drive:
    # Set drive amplitude
    zurich.output1_amp1(drive)
    time.sleep(5)  # Let the system settle

    # Create experiment name with drive info
    drive_uV = int(drive * 1e6)
    exp_name = f"Spectrum_{Temp:.3f}_drive{drive_uV}u"

    print(f"\n▶️ Running measurement at {drive_uV} μV drive...")

    # Run measurement
    run_thermomech_temp_meas(exp_name=exp_name)

    # Increment drive
    drive += step_drive


