import json

#general parameters

#device
device_name = 'CD12_B5_F4'

#saving parameters
costum_prefix='_'
costum_prefix=''
#zurich 

tc =  30e-3   # in seconds. Doesn't get overwritten by ZI called value.
tg=5e-3
attn_dB_source = 42.3+20 # attenuation at the source in dB
attn_dB_gate = 46+20
source_amplitude_instrumentlevel_GVg =20e-3
mix_down_f = 1.25e6 # RLC frequency
slew_rate=0.01

max_ramp_speed=9e-3#for long ramps
ramp_step_size=0.9e-1#for long ramps

pre_ramping_required=True

#compensation
x_avg=+1.24465881e-06#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-1.07161223e-06 #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100




#used params
start_vg_cs=0.794#1.12 #0.960
stop_vg_cs =0.8#1.17#0.970
step_num_cs=6*50#1000*5#10*100
#saved old
#for squeesedsingledot2
fit_type='data'#'tunnel_broadened'#'thermal'#'data'#'tunnel_broadened'
data_avg_num=15#make it ODD!!
sitfraction="l_max_slope"

min_acceptable_peak=50e-9

RLC_frequency=1.25e6

idt_point1_x=-1.51742
idt_point1_y=-2.25909
idt_point2_x=-1.50758
idt_point2_y=-2.25254

start_f =135e6#159.1e6 #Hz unit
stop_f =  145e6#159.3e6 #Hz unit
step_num_f =10*1000
freq_sweep_avg_num=21


#linesweep
start_vgo_ls=0.9#-1.88
stop_vgo_ls=1.6#-1.68
step_vgo_num_ls=700
start_vgi_ls= 0.72#0.920
stop_vgi_ls= 0.8#0.970
step_vgi_num_ls=100*10
start_vgi_scan_ls=0.771#966e-3
scan_range_ls=10e-3
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