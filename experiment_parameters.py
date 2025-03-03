import json

#general parameters

#device
device_name = 'CD11_D7_C1'

#saving parameters

#zurich 

tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
tg=5e-3
attn_dB_source = 42.3+20 # attenuation at the source in dB
attn_dB_gate = 46+20
source_amplitude_instrumentlevel_GVg = 20e-3
mix_down_f = 1.25e6 # RLC frequency
slew_rate=0.01

#compensation
x_avg=+4.6e-6  #+1.51e-5@75#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-6.6e-6  #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100


#used params
start_vg_cs = -1.647
stop_vg_cs = -1.644
step_num_cs=3*100


#saved old
#for squeesedsingledot2
fit_type='data'
data_avg_num=15#over 50uV
sitfraction=0.5#"l_max_slope"

min_acceptable_peak=50e-9

RLC_frequency=1.25e6

idt_point1_x=-1.51742
idt_point1_y=-2.25909
idt_point2_x=-1.50758
idt_point2_y=-2.25254


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