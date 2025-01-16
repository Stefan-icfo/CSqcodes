import json

#general parameters

#device
device_name = 'CD11_D7_C1'

#saving parameters

#zurich 

tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
attn_dB_source = 39 # attenuation at the source in dB
attn_dB_gate = 46+20
source_amplitude_instrumentlevel_GVg = 20e-3
mix_down_f = 1.25e6 # RLC frequency

#compensation
x_avg=+3.4e-6  #+1.51e-5@75#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-5.4e-6  #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100


#used params
start_vg_cs = -1.513
stop_vg_cs = -1.511
step_num_cs=2*10000

#saved old
#for squeesedsingledot2
fit_type='data'
data_avg_num=301
sitfraction=0.5#"l_max_slope"

