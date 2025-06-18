import json

#general parameters

#device
device_name = 'CD13_E3_C2'

#saving parameters
costum_prefix='_'
costum_prefix='RT'
#zurich 

tc = 30e-3   # in seconds. Doesn't get overwritten by ZI called value.
tg=5e-3
attn_dB_source = 42.3+20 # attenuation at the source in dB
attn_dB_gate = 46+20
source_amplitude_instrumentlevel_GVg = 20e-3
mix_down_f = 1.25e6 # RLC frequency
slew_rate=0.01

max_ramp_speed=5e-3#for long ramps
ramp_step_size=1e-1#for long ramps



#compensation
x_avg=+1.35e-5  #+1.51e-5@75#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-1.4e-5  #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100




#used params
start_vg_cs =-2
stop_vg_cs =3
step_num_cs=5000*10

#saved old
#for squeesedsingledot2
fit_type='data'
data_avg_num=11#make it ODD!!
sitfraction=0.5#"l_max_slope"

min_acceptable_peak=50e-9


RLC_frequency=1.25e6

idt_point1_x=-1.51742
idt_point1_y=-2.25909
idt_point2_x=-1.50758
idt_point2_y=-2.25254

start_f = 340e6 #Hz unit
stop_f =  440e6 #Hz unit
step_num_f =100000
freq_sweep_avg_num=11


#linesweep
start_vgo_ls=-0.62
stop_vgo_ls=-0.59
step_vgo_num_ls=30*10
start_vgi_ls= -0.477
stop_vgi_ls= -0.473
step_vgi_num_ls=4*50
start_vgi_scan_ls=-0.475
scan_range_ls=20e-3
increments_ls=0

"""

MAPI:


#used params
start_vg_cs = -1.655
stop_vg_cs = -1.645
step_num_cs=10*100


#saved old
#for squeesedsingledot2
fit_type='data'
data_avg_num=15#over 50uV
sitfraction=0.5#"l_max_slope"

min_acceptable_peak=50e-9

RLC_frequency=1.25e6

idt_point1_x=-1.52262
idt_point1_y=-2.26213
idt_point2_x=-1.51551
idt_point2_y=-2.25797
"""