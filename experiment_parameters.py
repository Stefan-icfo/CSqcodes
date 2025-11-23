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
source_amplitude_instrumentlevel_GVg = 0.04
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
start_vg_cs = 0.75
stop_vg_cs = 0.9
step_num_cs = 7500
sitside = 'left'


#GVg fitting
fit_type='data'#'tunnel_broadened'#'thermal'#'data'#'tunnel_broadened'
data_avg_num=15#make it ODD!!
sitfraction="l_max_slope"
min_acceptable_peak=50e-9


###########mechanics sweep######################



f_mech_opt_sitpos=153.54e6
#start_f = 150000000.0
#stop_f = 155000000.0
start_f = 151000000.0
stop_f = 155000000.0
step_num_f = 70*1000
#step_num_f =8*500#160*4#300*4#500*4#80*100#10000*4 #4000#
#step_num_f=round((stop_f-start_f)/250)


freq_sweep_avg_num=5


#linesweep
start_vgo_ls= 0#0.1127 
stop_vgo_ls=1.5
step_vgo_num_ls=150*5
scan_range_ls=6e-3#10e-3

start_vgi_ls= 0.750
stop_vgi_ls= 0.9
step_vgi_num_ls=150*10#300*10
start_vgi_scan_ls=0.845
increments_ls=0#not used I think

#####for cs_meta
increment_meta=-0.4
startpos_gate_meta=0.3#not used
startpos_auxgate_meta=0.8#not used


##########for finding mechanical mode, not in use yet#################
findM_start_drive=50e-3
findM_end_drive=200e-6
#freq_range=None,#this uses the generalmech_freuqency range
findM_found_range=1e6
findM_start_step_pitch = 250
findM_div_factor=2
findM_div_f=2
findM_min_sig_I=1.5e-12
findM_min_initial_sig_I=1.9e-12
findM_avg_num=1
#freq_bands=[[135e6,144e6],[150e6,154e6]]#full span
freq_bands = [[151000000.0, 155000000.0]]
#freq_bands=[[152e6,154e6]]#for first few e on l1946 at 150M

###for meta###
manual_thermomech_frequency=None
manual_background_set=None
update_therm_freq=False
therm_reps = 100
temp_meas_counts = 2
softening_pitch=40e-6
softening_reps=30
background_reps = 100
autocorr_reps = 2
therm_autocorr_pitch=50e-6
autocorr_Vg_pitch=50e-6#doubled, I think


driven_avg_num_meta=1
driven_range_meta=100e3
driven_pitch_meta=0.1e3
driven_amp_meta=1e-3

#pos_list = [1.12736, 1.27197,1.41379, 1.5584, 1.703, 1.8476, 1.98943, 2.12847, 2.2703, 2.41212, 2.54838, 2.68186, 2.82091, 2.96273, 3.09621, 3.22691, 3.36317, 3.49666, 3.62458, 3.7525, 3.88598]#for the ever-repeated 1946 ls, 271025
#pos_list=[0.471717, 0.693939, 0.885859, 1.06768]#tensioned config in broken database 21
pos_list_tensioned_attempt1_08=[0.470854, 0.694472, 0.885427, 1.06633]#same tensioned config in new database 22 - run 33
pos_list1=[0.453266, 0.661809, 0.850251, 1.03869]#35
#pos_list_for_ac=[ 1.26453, 1.40481, 1.5501, 2.4018, 2.53707, 2.67234, 2.81263, 3.61924, 3.74449, 3.87976,0.412826, 0.563126, 0.693387]
pos_list_for_ac=[ 2.53707, 2.67234, 2.81263, 1.26453, 1.40481, 1.5501, 2.4018, 3.61924, 3.74449, 3.87976,0.412826, 0.563126, 0.693387]#301025

pos_list_rep1946_on311025=[0.410505, 0.561978, 0.693995, 0.832961, 0.977486, 1.12201, 1.26931, 1.41106, 1.55419, 1.69872, 1.84324, 1.98499, 2.12256, 2.26431, 2.40605, 2.54224, 2.67704, 2.81739, 2.95775, 3.09116, 3.22178, 3.35658, 3.49277, 3.62062, 3.74707, 3.88048]
#3110
pos_15e_ls1946_rep=[2.40605]
pos_list_tensioned_0411=[0.608363, 0.754535, 0.919729, 1.08192, 1.2371]
pos_list_tensioned_attempt1_08_rep=[0.689778, 0.881161, 1.06352,0.453266]
pos_list_tensioned_attempt1_08_rep_extra_pt=[1.24852,1.43352]

pos_list_rep1946_on311025_softening=[   2.40605, 1.69872,3.09116,  3.74707 ]

pos_list_078_tensioned=[0.453266, 0.661809, 0.850251, 1.03869]

pos_list_tensioned_attempt1_08_rep_step3and4and5=[1.06352,1.24852,1.43352]
last_point_tensioned_attempt1_08_rep_step0=[0.453266]
pos_list=pos_list_tensioned_attempt1_08_rep_step3and4and5
#pos_list=pos_list_rep1946_on311025_softening

##############for moving dot##############, shape [[g3_firststep,g2_firststep,g1_firststep],[g3_secondstep,g2_secondstep,g1_secondstep],[...],...]
#pos_listg3h2g1=[[0,0.755,0.6525],[0.1,0.6925,0.666]]#these values are for 2 electrons. forget them for now

#pos_listg3h2g1=[[0,0.755,0.6525],[0.1,0.6925,0.666]]
pos_listg3h2g1_base = [
    [0,    0.585,   0.07205],
    [0.1,  0.535,   0.729],
    [0.2,  0.47,    0.7435],
    [0.3,  0.4025,  0.759],
    [0.4,  0.31,    0.7845],
    [0.45, 0.25,    0.80275],
    [0.5,  0.1875,  0.822],
    [0.55, 0.1325,  0.83825],
]
pos_listg3h2g1 = [
    [0,    0.585,   0.07605],
    [0.1,  0.535,   0.769],
    [0.2,  0.47,    0.7835],
    [0.3,  0.4025,  0.799],
    [0.4,  0.31,    0.8245],
    [0.45, 0.25,    0.84275],
    [0.5,  0.1875,  0.862],
    [0.55, 0.1325,  0.87825],
]#40mV higher; might try different angle eventually

mech_freq_list = [153.43e6, 153.27e6, 153.28e6, 153.30e6, 153.31e6, 153.3e6, 153.3e6, 153.3e6]#first six from previous measurements, last two are just guesses as of 221125

########################DQD params######################
idt_point1_x=-1.51742
idt_point1_y=-2.25909
idt_point2_x=-1.50758
idt_point2_y=-2.25254
DQD_stability_start_vg1 = 0.28
DQD_stability_start_vg2 = 0.23
