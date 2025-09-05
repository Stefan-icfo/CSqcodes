#from instruments import *

#measure nonlinearity 310725
import time
import copy
from instruments import *
#Vg=0#init
#from experiments.zurich_experiments.spectrum_0825 import *

side_pos=1.18468
top_pos=1.1843
qdac.ch06.dc_constant_V(side_pos)
time.sleep(10)
exp.mech_simple_fun_db(start_f=400e6,stop_f=500e6)
exp.mech_simple_fun_db(start_f=500e6,stop_f=600e6)
exp.mech_simple_fun_db(start_f=600e6,stop_f=700e6)

qdac.ch06.dc_constant_V(top_pos)#delft method
time.sleep(10)
exp.mech_simple_fun_db(start_f=100e6,stop_f=200e6,Delft=True)
exp.mech_simple_fun_db(start_f=200e6,stop_f=300e6,Delft=True)
exp.mech_simple_fun_db(start_f=300e6,stop_f=400e6,Delft=True)
exp.mech_simple_fun_db(start_f=400e6,stop_f=500e6,Delft=True)
exp.mech_simple_fun_db(start_f=500e6,stop_f=600e6,Delft=True)
exp.mech_simple_fun_db(start_f=600e6,stop_f=700e6,Delft=True)
exp.mech_simple_fun_db(start_f=20e6,stop_f=100e6,Delft=True)



