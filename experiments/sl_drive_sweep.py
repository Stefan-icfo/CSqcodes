#from instruments import *

#measure nonlinearity 310725
import time
import copy
from instruments import *
#Vg=0#init
#from experiments.zurich_experiments.spectrum_0825 import *

start_drive=20e-3#100e-6#starting value drive amp
multiplier=1.5
end_drive=750e-3



current_drive=copy.copy(start_drive)
zurich.output0_amp0(current_drive)
time.sleep(10)#waiting for cs gate



while current_drive<=end_drive:
        #zurich.output0_amp0(current_drive)
        exp.set_params(source_amplitude_instrumentlevel_GVg=current_drive)
        print(f"setting drive to {current_drive*1e3} mV")
        time.sleep(10)
        
        exp.GVG_fun_sensitivity(costum_prefix=f"adjust_drive_{current_drive*1e3:.3g}mV")
        current_drive=current_drive*multiplier
        
        time.sleep(5)

zurich.output0_amp0(start_drive)