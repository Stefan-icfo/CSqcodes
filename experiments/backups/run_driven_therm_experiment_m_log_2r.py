import numpy as np
import time
from experiments.zurich_experiments.spectrum_0825 import run_thermomech_temp_meas
from instruments import zurich,qdac,exp

# Spectrum acquisition parameters
#filter_bw = 50e3
#rbw = 3.353
#BURST_DURATION = 298.262e-3
#SAMPLING_RATE = 219.72656250e3
Temp=0.13
nr_bursts = 7
reps_nodrive = 10
demod_ch = 3
#avg_num = 14035

# Drive loop parameters
drive = 30e-6       # Start at 10 μV
end_drive = 20e-3   # End at 30 mV
factor = 1.2      # +100% increment per step

exp.sit_at_max_Isens()
#qdac.ch06.dc_constant_V(0.966)#setting to the other side
time.sleep(50)
zurich.set_mixdown(151.14983000e6)
while drive <= end_drive:
    zurich.output1_amp1(drive)
    time.sleep(5)  # Let system settle

    drive_uV = int(round(drive * 1e6))
    exp_name = f"Spectrum_{Temp:.3f}_after_disaster_{drive_uV}u"

    print(f"\n▶️ Running measurement at {drive_uV} μV drive...")

    run_thermomech_temp_meas(
        reps_nodrive=reps_nodrive,
        exp_name=exp_name,
        fit_lorentzian=False
    )
    
    drive *= factor

exp_name = f"Spectrum_{Temp:.3f}_thermal"
zurich.output1_amp1(0)
run_thermomech_temp_meas(
        reps_nodrive=reps_nodrive,
        exp_name=exp_name,
        fit_lorentzian=False
    )
zurich.move_mixdown(-8e6)
exp_name = f"Spectrum_{Temp:.3f}_background"
zurich.output1_amp1(0)
run_thermomech_temp_meas(
        reps_nodrive=reps_nodrive,
        exp_name=exp_name,
        fit_lorentzian=False
    )




