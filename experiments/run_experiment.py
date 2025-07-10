from instruments import *

while exp.stop_f<700e6:
    exp.do_GVg_and_adjust_sitpos()
    zurich.sigout1_amp1_enabled_param.value(1)
    exp.mech_simple_fun_db()
    zurich.sigout1_amp1_enabled_param.value(0)
    exp.start_f+=100e6
    if exp.stop_f<600e6:
        exp.stop_f+=100e6 
    else:
       exp.stop_f=700e6