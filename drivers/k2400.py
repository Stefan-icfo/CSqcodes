
import numpy as np
from time import sleep
from tqdm import tqdm


def ramp_k2400(ramp_param, final_vg:float = 0, step_size:float = 1e-3,  ramp_speed:float = 100e-6, *args, **kwargs):
    end_point = final_vg
    start_point = ramp_param.volt()
    wait_time = (step_size/ramp_speed) *1e-3 # in seconds
    if start_point>=end_point:
        step_size = -step_size
    ampl_sweep = np.append(np.arange(start_point,end_point,step_size),end_point)
    y0 = end_point
    x0 = start_point
    dY = y0-x0
    
    if (abs(dY)>1e-6):
        for i in ampl_sweep:
            ramp_param.volt(i)
            sleep(wait_time)
            # x = gate.volt()
            # status = (x-x0)/dY
            # tqdm.write(f"{status*100:.0f}% done", end='\r')
        
        # final_read = gate.volt()
        # print(f'Ramped to {final_read} V \n')
