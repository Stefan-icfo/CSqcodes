#from instruments import *

import time
import copy
from instruments import *
#Vg=0#init
#start_V=0.9065#1.14#starting value Vgcs
current_A=10e-3#starting value drive amp
min_a=2e-3


#max_V=0.9015#1.148

mech_freq=136730885
avg_num=10

Vg,G,sens=exp.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=True)
peakpos=Vg[np.argmax(G)]
print(f"setting cs to {peakpos*1e3:.5g} mV")
approx_maxpos=peakpos
start_V=0.77775-0.3e-3
max_V=0.77775+0.3e-3
A_factor=2
#exp.sit_at_max_Isens()
f_mech_opt_sitpos=mech_freq
start_f=f_mech_opt_sitpos-100e3
stop_f=f_mech_opt_sitpos+300e3
step_num_f=round((stop_f-start_f)/250)

while current_A>min_a:
    current_V=copy.copy(start_V)
   # qdac.ch06.dc_constant_V(start_V)
   # time.sleep(10)#waiting for cs gate

    qdac.ch06.dc_constant_V(start_V)
    zurich.output1_amp1(current_A)
    time.sleep(10)
    while current_V<max_V:
        qdac.ch06.dc_constant_V(current_V)
        for n in range(avg_num):
            exp.mech_simple_fun_db(start_f=start_f,stop_f=stop_f,step_num_f=step_num_f,costum_prefix=f"driven_vs_sitpos_{current_V*1e3:6g}mV_g3")
        current_V=qdac.ch06.dc_constant_V()
        current_V+=0.2e-4
        
        print(f"setting ch06  to {current_V:6g} mV")
        time.sleep(5)
    current_A=current_A/A_factor

exp.sit_at_max_Isens()