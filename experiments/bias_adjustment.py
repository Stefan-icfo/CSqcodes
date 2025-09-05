#from instruments import *

#measure nonlinearity 310725
import time
import copy
from instruments import *
#Vg=0#init
#from experiments.zurich_experiments.spectrum_0825 import *

start_bias=100e-6#starting value drive amp
bias_step=20e-6
end_bias=300e-6
saved_original_bias=qdac.ch07.dc_constant_V()
print(f"saved_original_bias {saved_original_bias}")

current_bias=copy.copy(start_bias)
qdac.ch07.dc_constant_V(start_bias)
time.sleep(10)#waiting for cs gate



while current_bias<end_bias:
        print(f"setting ch07 bias to {current_bias*1e6} uV")
        if abs(current_bias)>end_bias:
                raise ValueError("bias too high")
        time.sleep(10)
        qdac.ch07.dc_constant_V(current_bias)
        exp.GVG_fun_sensitivity(costum_prefix="adjust_bias")
        current_bias+=bias_step
        
        
        time.sleep(5)
print(f"setting back to saved_original_bias {saved_original_bias}")
qdac.ch07.dc_constant_V(saved_original_bias)