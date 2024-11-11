from instruments import qdac
import numpy as np
from time import sleep
from tqdm import tqdm

def read_channels(chan_nr = 7):
    ch1=qdac.ch01.dc_constant_V()
    ch2=qdac.ch02.dc_constant_V()
    ch3=qdac.ch03.dc_constant_V()
    ch4=qdac.ch04.dc_constant_V()
    ch5=qdac.ch05.dc_constant_V()
    ch6=qdac.ch06.dc_constant_V()
    ch7=qdac.ch07.dc_constant_V()
    sl_ch1=qdac.ch01.dc_slew_rate_V_per_s()
    sl_ch2=qdac.ch02.dc_slew_rate_V_per_s()
    sl_ch3=qdac.ch03.dc_slew_rate_V_per_s()
    sl_ch4=qdac.ch04.dc_slew_rate_V_per_s()
    sl_ch5=qdac.ch05.dc_slew_rate_V_per_s()
    sl_ch6=qdac.ch06.dc_slew_rate_V_per_s()
    sl_ch7=qdac.ch07.dc_slew_rate_V_per_s()
    fl_ch1=qdac.ch01.output_filter()
    fl_ch2=qdac.ch02.output_filter()
    fl_ch3=qdac.ch03.output_filter()
    fl_ch4=qdac.ch04.output_filter()
    fl_ch5=qdac.ch05.output_filter()
    fl_ch6=qdac.ch06.output_filter()
    fl_ch7=qdac.ch07.output_filter()
    channels=[ch1,ch2,ch3,ch4,ch5,ch6,ch7]
    slewrates=[sl_ch1,sl_ch2,sl_ch3,sl_ch4,sl_ch5,sl_ch6,sl_ch7]
    filters=[fl_ch1,fl_ch2,fl_ch3,fl_ch4,fl_ch5,fl_ch6,fl_ch7]
    for i in range(chan_nr):
        print(f"qdac.channel {i+1} is {channels[i]} V,its slew rate is {slewrates[i]} V/s and the filter is {filters[i]}")
        #print(f"qdac.channel {i+1} slew rate is {slewrates[i]} V/s")

def set_all_slewrates(slew_rate = 0.01):
    qdac.ch01.dc_slew_rate_V_per_s(slew_rate)
    qdac.ch02.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch03.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch04.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch05.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch06.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch07.dc_slew_rate_V_per_s(slew_rate) 


def ramp_to_zero(chan_nr = 7,slew_rate = 0.01):
    qdac.ch01.dc_slew_rate_V_per_s(slew_rate)
    qdac.ch02.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch03.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch04.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch05.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch06.dc_slew_rate_V_per_s(slew_rate) 
    qdac.ch07.dc_slew_rate_V_per_s(slew_rate)
    qdac.ch01.dc_constant_V(0)
    qdac.ch02.dc_constant_V(0)
    qdac.ch03.dc_constant_V(0)
    qdac.ch04.dc_constant_V(0)
    qdac.ch05.dc_constant_V(0)
    qdac.ch06.dc_constant_V(0)
    qdac.ch07.dc_constant_V(0)
    


def ramp_QDAC_channel(ramp_param, slew_rate = 1e-2,final_vg:float = 0, step_size:float = 10e-3,  ramp_speed:float = 1e-3, *args, **kwargs):
    #ramp_param.dc_slew_rate_per_s(slew_rate)
    
    end_point = final_vg
    start_point = ramp_param.dc_constant_V()
    print(f"sleep time={(final_vg-start_point)/ramp_speed}")
    wait_time = (step_size/ramp_speed)  # in seconds
    if start_point>=end_point:
        step_size = -step_size
    V_sweep = np.append(np.arange(start_point,end_point,step_size),end_point)
    y0 = end_point
    x0 = start_point
    dY = y0-x0
    
    if (abs(dY)<=1e-6):
            print(f'No change in voltage \r')
    else:
        for i in V_sweep:
            ramp_param.dc_constant_V(i)
            time.sleep(wait_time)
            

def ramp_QDAC_multi_channel(ramp_params, final_vgs, slew_rate = 1e-2, step_size:float = 10e-3,  ramp_speed:float = 1e-3, *args, **kwargs):
    #ramp_param.dc_slew_rate_per_s(slew_rate)
    end_points = final_vgs
    
    wait_time = (step_size/ramp_speed)  # in seconds
    
    V_sweeps=[]
    step_nums=[]
    start_points=[]
    j=0
    for ramp_param in ramp_params:
        start_points=start_points+[ramp_param.dc_constant_V()]
        step_nums=step_nums+[round(abs(start_points[j]-end_points[j])/step_size)]
        j=j+1
    step_num=max(step_nums)
    print(f"step_num={step_num}")
    j=0
    for ramp_param in ramp_params:
        V_sweeps = V_sweeps+[np.linspace(start_points[j],end_points[j],num=step_num)]
        j=j+1
        
    
    for i in tqdm(range(step_num)):
            
        for j in range(len(ramp_params)):
                ramp_params[j].dc_constant_V(V_sweeps[j][i])
        sleep(wait_time)

#def combined_gate_scan(gates = [qdac.ch02,qdac.ch04], starting?point):

def ramp_time(gates,setpoints,slew_rate=0.01):#in construction
     single_gate_ramp_times=[]
     for gate,setpoint in zip(gates,setpoints):
         single_gate_ramp_times.append(abs(gate.dc_constant_V()-setpoint)/slew_rate)
     return max(single_gate_ramp_times)