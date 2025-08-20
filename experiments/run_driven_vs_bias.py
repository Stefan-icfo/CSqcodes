#from instruments import *

#measure nonlinearity 310725
import time
import copy
from instruments import *
#Vg=0#init
from experiments.zurich_experiments.spectrum_0825 import *

start_bias=0
end_bias=300e-6#starting value drive amp
bias_step=20e-6
drive_V=300e-6
zurich.output1_amp1(drive_V)
saved_original_bias=qdac.ch07.dc_constant_V()
print(f"saved_original_bias {saved_original_bias}")

current_bias=copy.copy(start_bias)
qdac.ch07.dc_constant_V(start_bias)
time.sleep(10)#waiting for cs gate

mech_freq=137.589e6
zurich.set_mixdown(mech_freq)

while current_bias<end_bias:
        print(f"setting ch07 bias to {current_bias*1e6} uV")
        if abs(current_bias)>300e-6:
                raise ValueError("bias too high")
        time.sleep(10)
        qdac.ch07.dc_constant_V(current_bias)
        exp.sit_at_max_Isens()
        zurich.set_mixdown(mech_freq)
        time.sleep(100)
        run_thermomech_temp_meas(exp_name=f'driven_vs_bias={current_bias*1e6} uV_test')
        
        current_bias+=bias_step
        
        
        time.sleep(5)
print(f"setting back to saved_original_bias {saved_original_bias}")
qdac.ch07.dc_constant_V(saved_original_bias)