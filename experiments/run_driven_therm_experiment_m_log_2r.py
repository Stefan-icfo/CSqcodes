import numpy as np
import time
from experiments.zurich_experiments.spectrum_0825 import run_thermomech_temp_meas
from instruments import zurich

# Spectrum acquisition parameters
filter_bw = 50e3
rbw = 3.353
BURST_DURATION = 298.262e-3
SAMPLING_RATE = 219.72656250e3
nr_bursts = 7
reps_nodrive = 50
demod_ch = 3
avg_num = 21
Temp = 0.035

# Drive loop parameters
drive = 10e-6       # Start at 10 μV
end_drive = 30e-3   # End at 30 mV
factor = 1.2        # +20% increment per step

while drive <= end_drive:
    zurich.output1_amp1(drive)
    time.sleep(5)  # Let system settle

    drive_uV = int(round(drive * 1e6))
    exp_name = f"Spectrum_{Temp:.3f}_notmaxsensitivity_{drive_uV}u"

    print(f"\n▶️ Running measurement at {drive_uV} μV drive...")

    run_thermomech_temp_meas(
        reps_nodrive=reps_nodrive,
        exp_name=exp_name,
        fit_lorentzian=False
    )
    
    drive *= factor






