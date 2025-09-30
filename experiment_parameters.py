import json

#general parameters

#device
device_name = 'CD12_B5_F4'

#saving parameters
costum_prefix='_'
costum_prefix=''
#zurich 

tc =  30e-3#100e-3   # in seconds. Doesn't get overwritten by ZI called value.
tg=5e-3
attn_dB_source = 42.3+20 # attenuation at the source in dB
attn_dB_gate = 46+20
source_amplitude_instrumentlevel_GVg =20e-3
mix_down_f = 1.25e6 # RLC frequency
slew_rate=0.01

max_ramp_speed=5e-3#for long ramps
ramp_step_size=5e-2#for long ramps

pre_ramping_required=True

#compensation
x_avg=3.72811964e-06#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-3.01203684e-06 #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100




#used params
start_vg_cs=0.85#0.955#0.905#0.870#0.783#0.803#1.12 #0.960
stop_vg_cs =0.95#0.975#0.910#0.875#11#1.17#0.970
step_num_cs=100*50#1000*5#10*100
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


f_mech_opt_sitpos=146.75e6
start_f=f_mech_opt_sitpos-500e3
stop_f=f_mech_opt_sitpos+500e3
step_num_f=round((stop_f-start_f)/250)
start_f =120e6#141.4e6#141.6e6#150.5e6#150e6#159.1e6 #Hz unit #438e6#
stop_f =  170e6#142.2e6#142.1e6#151.5e6#160e6#159.3e6 #Hz unit #442e6#
step_num_f =5000#160*4#300*4#500*4#80*100#10000*4 #4000#
#step_num_f=round((stop_f-start_f)/250)
freq_sweep_avg_num=21


#linesweep
start_vgo_ls=-4#-1.88
stop_vgo_ls=4#-1.68
step_vgo_num_ls=800#800*2#5002mV
start_vgi_ls= 0.700#0.920
stop_vgi_ls= 1#0.970
step_vgi_num_ls=300*10
start_vgi_scan_ls=0.89#966e-3
scan_range_ls=6e-3
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