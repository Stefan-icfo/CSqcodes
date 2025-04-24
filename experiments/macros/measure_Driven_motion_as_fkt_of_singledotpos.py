from instruments import exp,qdac,zurich
import time
#offset=3e-3
half_range=50e-3
slopes=[]
for n in range(10):
    slope,max_pos=exp.do_GVg_and_adjust_sitpos(sitfraction=0.87,sit_side="right")
    exp.start_vg_cs=max_pos-50
    exp.stop_vg_cs=max_pos+50
    slopes.append(slope)
    time.sleep(10)
    #max_pos=qdac.ch06.dc_constant_V()
    #qdac.ch06.dc_constant_V(max_pos+offset)
    time.sleep(10)
    exp.mech_simple_fun_db()
    current_g1=qdac.ch01.dc_constant_V()
    current_g3=qdac.ch03.dc_constant_V()
    time.sleep(1)
    qdac.ch01.dc_constant_V(current_g1+0.02)
    qdac.ch03.dc_constant_V(current_g3-0.02)
    time.sleep(1)
    print(f"start and stop vgcs{exp.start_vg_cs:4g} and {exp.stop_vg_cs:4g}, max_pos {max_pos:4g}")
    
    
