#from instruments import *

import time
import copy
from instruments import *
#Vg=0#init
start_V=0.927#1.14#starting value Vgcs
current_A=5e-3#starting value drive amp
#max_a=20e-3
max_V=0.929
#A_factor=2
exp.sit_at_max_Isens()#put in proximity

#while current_A<max_a:
#    current_V=copy.copy(start_V)
#    qdac.ch06.dc_constant_V(start_V)
#    time.sleep(10)#waiting for cs gate
current_V=copy.copy(start_V)
qdac.ch06.dc_constant_V(start_V)
zurich.output1_amp1(current_A)
time.sleep(10)
while current_V<max_V:
        qdac.ch06.dc_constant_V(current_V)
        exp.mech_simple_fun_db(costum_prefix=f"driven_vs_sitpos_{current_V*1e3:6g}mV__28holes")
        current_V=qdac.ch06.dc_constant_V()
        current_V+=0.04e-3
        
        print(f"setting ch06  to {current_V:6g} mV")
        time.sleep(5)
    #current_A=current_A*A_factor

exp.sit_at_max_Isens()