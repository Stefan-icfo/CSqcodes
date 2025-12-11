import json

#######################general parameters#############################

#device
device_name = 'CD12_B5_F4'

#saving parameters
costum_prefix='_'
tc = 0.03
tg=5e-3
attn_dB_source = 42.3+20 # attenuation at the source in dB
attn_dB_gate = 46+20
source_amplitude_instrumentlevel_GVg = 0.02
mix_down_f = 1.25e6 # RLC frequency
slew_rate=0.01

max_ramp_speed=5e-3#for long ramps
ramp_step_size=5e-2#for long ramps

RLC_frequency=1.25e6

#compensation
x_avg=3.72811964e-06#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-3.01203684e-06 #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100


##################measured cross-capacitances###########

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
start_vg_cs = 0.7
stop_vg_cs = 1.2
step_num_cs = 12500
sitside = 'left'


#GVg fitting
fit_type='data'#'tunnel_broadened'#'thermal'#'data'#'tunnel_broadened'
data_avg_num=15#make it ODD!!
sitfraction="l_max_slope"
min_acceptable_peak=50e-9


###########mechanics sweep######################



#f_mech_opt_sitpos=153.54e6
#start_f = 150000000.0
#stop_f = 155000000.0
start_f = 150800000.0
stop_f = 154800000.0
step_num_f =20*1000*2
#step_num_f =8*500#160*4#300*4#500*4#80*100#10000*4 #4000#
#step_num_f=round((stop_f-start_f)/250)


freq_sweep_avg_num=5


#linesweep
start_vgo_ls=-2#0.1127 
stop_vgo_ls=2
step_vgo_num_ls=80*2
scan_range_ls = 10e-3
start_vgi_ls = 1.7
stop_vgi_ls = 1.8
step_vgi_num_ls = 100*5
start_vgi_scan_ls=0.87
increments_ls=0#not used I think

#####for cs_meta
increment_meta=-0.4
startpos_gate_meta=0.3#not used
startpos_auxgate_meta=0.8#not used


##########for finding mechanical mode, not in use yet#################
findM_start_drive=75e-3
findM_end_drive=5e-3
#freq_range=None,#this uses the generalmech_freuqency range
findM_found_range=1e6
findM_start_step_pitch = 500
findM_div_factor=2
findM_div_f=2
findM_min_sig_I=1.5e-12
findM_min_initial_sig_I=1.9e-12
findM_avg_num=1
#freq_bands=[[135e6,144e6],[150e6,154e6]]#full span
freq_bands = [[149000000.0, 156000000.0], [140000000.0, 144000000.0]]
#freq_bands=[[152e6,154e6]]#for first few e on l1946 at 150M

###for meta###
manual_thermomech_frequency=None
manual_background_set=None
update_therm_freq=False
therm_reps = 50
temp_meas_counts = 3
softening_pitch=50e-6
softening_reps=30
background_reps = 50
autocorr_reps = 1
therm_autocorr_pitch=50e-6
autocorr_Vg_pitch=50e-6#doubled, I think


driven_avg_num_meta=1
driven_range_meta=100e3
driven_pitch_meta=0.1e3
driven_amp_meta=1e-3




##############for moving dot##############, shape [[g3_firststep,g2_firststep,g1_firststep],[g3_secondstep,g2_secondstep,g1_secondstep],[...],...]




mech_freq_list = [153.315e6]

pos_listg3h2g1_for_softening= [
    [0.0,    0.585,   0.7205],
    [0.1,    0.535,   0.769]]
    #finegrained for softening
pos_listg3h2g1_fine = [
    [0.0,    0.585,   0.7205],
    [0.025,  0.5725,  0.732625],
    [0.05,   0.56,    0.74475],
    [0.075,  0.5475,  0.756875],
    [0.1,    0.535,   0.769],
    [0.125,  0.51875, 0.77263],
    [0.15,   0.5025,  0.77625]]


pos_listg3h2g1_holes= [
    [0.3,    0.58,   0.3],
    [0.3,    0.276,   0.3],
    [0.3,    -0.223,   0.3],
    [0.3,    -0.648,   0.3],
    [0.3,    -1.62,   0.3],
    [0.3,    -2.03,   0.3],
    [0.3,    -2.49,   0.3],
    [0.3,    -2.555,   0.3]]

pos_list=pos_listg3h2g1_fine

pos_listg3h2g1=pos_listg3h2g1_fine

pos_list_5g_freq=[
    [2.0, 2.0, 2.0, 2.0, 2.0],
    [1.8, 1.8, 1.8, 1.8, 1.8],
    [1.6, 1.6, 1.6, 1.6, 1.6],
    [1.4, 1.4, 1.4, 1.4, 1.4],
    [1.2, 1.2, 1.2, 1.2, 1.2],
    [1.0, 1.0, 1.0, 1.0, 1.0],
    [0.8, 0.8, 0.8, 0.8, 0.8],
    [0.6, 0.6, 0.6, 0.6, 0.6],
    [0.4, 0.4, 0.4, 0.4, 0.4],
    [0.2, 0.2, 0.2, 0.2, 0.2],
    [0.0, 0.0, 0.0, 0.0, 0.0],
    [-0.2, -0.2, -0.2, -0.2, -0.2],
    [-0.4, -0.4, -0.4, -0.4, -0.4],
    [-0.6, -0.6, -0.6, -0.6, -0.6],
    [-0.8, -0.8, -0.8, -0.8, -0.8],
    [-1.0, -1.0, -1.0, -1.0, -1.0],
    [-1.2, -1.2, -1.2, -1.2, -1.2],
    [-1.4, -1.4, -1.4, -1.4, -1.4],
    [-1.6, -1.6, -1.6, -1.6, -1.6],
    [-1.8, -1.8, -1.8, -1.8, -1.8],
    [-2.0, -2.0, -2.0, -2.0, -2.0]
]




pos_list_5g_alt1=[
    [1.0, 1.0, 1.0, 1.0, 1.0],
    [0.8, 0.8, 0.8, 0.8, 0.8],
    [0.6, 0.6, 0.6, 0.6, 0.6],
    [0.4, 0.4, 0.4, 0.4, 0.4],
    [0.2, 0.2, 0.2, 0.2, 0.2],
    [0.0, 0.0, 0.0, 0.0, 0.0],
    [-0.2, -0.2, -0.2, -0.2, -0.2],
    [-0.4, -0.4, -0.4, -0.4, -0.4],
    [-0.6, -0.6, -0.6, -0.6, -0.6],
    [-0.8, -0.8, -0.8, -0.8, -0.8],
    [-1.0, -1.0, -1.0, -1.0, -1.0],
    [-1.2, -1.2, -1.2, -1.2, -1.2],
    [-1.4, -1.4, -1.4, -1.4, -1.4],
    [-1.6, -1.6, -1.6, -1.6, -1.6],
    [-1.8, -1.8, -1.8, -1.8, -1.8],
    [-2.0, -2.0, -2.0, -2.0, -2.0]
]


pos_list_5g_alt2=[
    [2.0, 2.0, 2.0, 2.0, 2.0],
    [1.8, 1.8, 1.8, 1.8, 1.8],
    [1.6, 1.6, 1.6, 1.6, 1.6],
    [1.4, 1.4, 1.4, 1.4, 1.4],
    [1.2, 1.2, 1.2, 1.2, 1.2],
    [1.0, 1.0, 1.0, 1.0, 1.0]
]
cs_ranges=[[0.8,0.85],[1.44,1.49],[1.75,1.8],[1.88,1.93]]



########################DQD params######################
idt_point1_x=-1.51742
idt_point1_y=-2.25909
idt_point2_x=-1.50758
idt_point2_y=-2.25254
DQD_stability_start_vg1 = 0.28
DQD_stability_start_vg2 = 0.23
