import json

#######################general parameters#############################

#device
device_name = 'CD12_B5_F4_test'

#saving parameters
costum_prefix='_'


tc =  30e-3#100e-3   # in seconds. Doesn't get overwritten by ZI called value.
tg=5e-3
attn_dB_source = 42.3+20 # attenuation at the source in dB
attn_dB_gate = 46+20
source_amplitude_instrumentlevel_GVg =20e-3
mix_down_f = 1.25e6 # RLC frequency
slew_rate=0.01

max_ramp_speed=5e-3#for long ramps
ramp_step_size=5e-2#for long ramps

RLC_frequency=1.25e6

#compensation
x_avg=3.72811964e-06#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-3.01203684e-06 #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100


##################device specific params: in construction###########3333

gapcenters_5g=[0.8,0.3,0.3,0.3,0.3]
gapcenter_cs=0.5 #check

crosscap_5g_onto_g1=[0.4,0.115,0,0]
crosscap_5g_onto_g2=[0,0,0,0.25]
crosscap_5g_onto_g3=[0,0,0,0]
crosscap_5g_onto_g4=[0,0.25,0,0]
crosscap_5g_onto_g5=[0,0,0,0]

crosscap_5g_onto_cs=[0,0.019,0.012,0]

crosscapg2g4=0.25 #for first try
#################run-specific parameters SQD###############################################
pre_ramping_required=True

#GVg params
start_vg_cs=0.760#0.94#0.954#0.955#0.905#0.870#0.783#0.803#1.12 #0.960
stop_vg_cs =0.860#1.04#0.958#0.975#0.910#0.875#11#1.17#0.970
step_num_cs=40*50#4*50#1000*5#10*100

sitside="left"


#GVg fitting
fit_type='data'#'tunnel_broadened'#'thermal'#'data'#'tunnel_broadened'
data_avg_num=15#make it ODD!!
sitfraction="l_max_slope"

min_acceptable_peak=50e-9


###########mechanics sweep######################



f_mech_opt_sitpos=153.54e6
start_f=f_mech_opt_sitpos-100e3
stop_f=f_mech_opt_sitpos+100e3

start_f =152e6
stop_f =155e6# 1136.8e6#
step_num_f=15*1000
#step_num_f =8*500#160*4#300*4#500*4#80*100#10000*4 #4000#
#step_num_f=round((stop_f-start_f)/250)

#sstep_num_f=round((stop_f-start_f)/1000)
#sf_mech_opt_sitpos=153.54e6
#start_f=f_mech_opt_sitpos-5e6
#stop_f=f_mech_opt_sitpos+500e3

freq_sweep_avg_num=5


#linesweep
start_vgo_ls= -0.5#0.1127 
stop_vgo_ls=4#-1.68
step_vgo_num_ls=450*5

#start_vgo_ls= 2.53707 
#stop_vgo_ls=0.3
#step_vgo_num_ls=220
#800
#5002mV
start_vgi_ls= 0.75
stop_vgi_ls= 0.9
step_vgi_num_ls=150*10#300*10
start_vgi_scan_ls=0.8615
#966e-3
scan_range_ls=6e-3#10e-3
increments_ls=0

#####for cs_meta
increment_meta=-0.4
startpos_gate_meta=4
startpos_auxgate_meta=-0.67


##########for finding mechanical mode, not in use yet#################
findM_start_drive=75e-3
findM_end_drive=200e-6
#freq_range=None,#this uses the generalmech_freuqency range
findM_found_range=1e6
findM_start_step_pitch=0.25e3#0.25e6 for first few e on 1946 ls
findM_div_factor=4
findM_div_f=2
findM_min_sig_I=1.5e-12
findM_min_initial_sig_I=1.5e-12
findM_avg_num=1
#freq_bands=[[135e6,144e6],[150e6,154e6]]#full span
#freq_bands=[[135e6,139e6]]#reduced span for first few e
freq_bands=[[152e6,154e6]]#for first few e on l1946 at 150M

###for meta###

therm_reps=40
temp_meas_counts=3



########################DQD params######################
idt_point1_x=-1.51742
idt_point1_y=-2.25909
idt_point2_x=-1.50758
idt_point2_y=-2.25254