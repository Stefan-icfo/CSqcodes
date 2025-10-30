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

###Julie's cross-capactiance-matrix
cc_Julie=[[1,         0.61643226, 0.28104037, 0.14392452, 0.08360364, 0.10567371],
[0.73769445, 1,         0.66117359, 0.320361,   0.16656599, 0.07562857],
[0.34698729, 0.6877848,  1,         0.68721441, 0.34055084, 0.05488676],
[0.16999987, 0.32066261, 0.66124586, 1,         0.72614275, 0.03970722],
[0.08558164, 0.14441613, 0.28221459, 0.62090503, 1,         0.02646617],
[0.09056547, 0.06373444, 0.0424606,  0.02988223, 0.0222321,  1        ]]

#################run-specific parameters SQD###############################################
pre_ramping_required=True

#GVg params
start_vg_cs=0.83#0.94#0.954#0.955#0.905#0.870#0.783#0.803#1.12 #0.960
stop_vg_cs =0.86#1.04#0.958#0.975#0.910#0.875#11#1.17#0.970
step_num_cs=30*50#4*50#1000*5#10*100

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

start_f =145e6
stop_f =155e6# 1136.8e6#
step_num_f=10*1000*4
#step_num_f =8*500#160*4#300*4#500*4#80*100#10000*4 #4000#
#step_num_f=round((stop_f-start_f)/250)

#sstep_num_f=round((stop_f-start_f)/1000)
#sf_mech_opt_sitpos=153.54e6
#start_f=f_mech_opt_sitpos-5e6
#stop_f=f_mech_opt_sitpos+500e3

freq_sweep_avg_num=5


#linesweep
start_vgo_ls= 0.3#0.1127 
stop_vgo_ls=-1.3#-1.68
step_vgo_num_ls=100*2

#start_vgo_ls= 2.53707 
#stop_vgo_ls=0.3
#step_vgo_num_ls=220
#800
#5002mV
start_vgi_ls= 0.8
stop_vgi_ls= 0.9
step_vgi_num_ls=100*10#300*10
start_vgi_scan_ls=0.8615
#966e-3
scan_range_ls=6e-3#10e-3
increments_ls=0

#####for cs_meta
increment_meta=-0.4
startpos_gate_meta=4#not used
startpos_auxgate_meta=-0.67#not used


##########for finding mechanical mode, not in use yet#################
findM_start_drive=75e-3
findM_end_drive=200e-6
#freq_range=None,#this uses the generalmech_freuqency range
findM_found_range=1e6
findM_start_step_pitch=0.25e3#0.25e6 for first few e on 1946 ls
findM_div_factor=4
findM_div_f=2
findM_min_sig_I=1.5e-12
findM_min_initial_sig_I=1.8e-12
findM_avg_num=1
#freq_bands=[[135e6,144e6],[150e6,154e6]]#full span
freq_bands=[[141e6,144e6]]#reduced span for first few e
#freq_bands=[[152e6,154e6]]#for first few e on l1946 at 150M

###for meta###

therm_reps=40
temp_meas_counts=3
softening_pitch=0.05e-3
softening_reps=20
background_reps=80

#pos_list = [1.12736, 1.27197,1.41379, 1.5584, 1.703, 1.8476, 1.98943, 2.12847, 2.2703, 2.41212, 2.54838, 2.68186, 2.82091, 2.96273, 3.09621, 3.22691, 3.36317, 3.49666, 3.62458, 3.7525, 3.88598]#for the ever-repeated 1946 ls, 271025
#pos_list=[0.471717, 0.693939, 0.885859, 1.06768]#tensioned config in broken database 21
pos_list2=[0.470854, 0.694472, 0.885427, 1.06633]#same tensioned config in new database 22 - run 33
pos_list1=[0.453266, 0.661809, 0.850251, 1.03869]#35
pos_list_for_ac=[2.52836]

pos_list=pos_list1

########################DQD params######################
idt_point1_x=-1.51742
idt_point1_y=-2.25909
idt_point2_x=-1.50758
idt_point2_y=-2.25254