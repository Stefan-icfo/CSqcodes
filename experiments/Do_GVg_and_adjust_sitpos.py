

import numpy as np


from instruments import station, zurich, qdac,Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
from utils.CS_utils import *
from experiments.GVg_qdac_zurich_general import GVG_fun
import scipy as scp
import time
from tqdm import tqdm
from experiment_functions.CS_functions import *

#------User input----------------

#adjustable hardware params

tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd_dB = 39 # attenuation at the source in dB
source_amplitude_instrumentlevel_GVg = 20e-3
mix_down_f = 1.25e6 # RLC frequency

#channel assignment
gate=qdac.ch06
freq = zurich.oscs.oscs0.freq
source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
measured_parameter = zurich.demods.demods0.sample 

#compensation
x_avg=+3.4e-6  #+1.51e-5@75#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=-5.4e-6  #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100


#gate sweep params
start_vg = -2.5
stop_vg = 0
step_num= 2500*10

pre_ramping_required=True

#costum name
device_name = 'CD11_D7_C1'
prefix_name = 'Conductance_rf_'
exp_name="_"

#fixed hardware params
#####################
gain_RT = 200       
gain_HEMT = 5.64   
Z_tot = 7521       
###################

def do_GVg_and_adjust_sitpos(
            fit_type='tunnel_broadened',
            initial_guess=None, 
            sitfraction="l_max_slope",
            start_vg=start_vg,
            stop_vg=stop_vg,
            tc=tc,
            source_amplitude_instrumentlevel_GVg=source_amplitude_instrumentlevel_GVg,
            x_avg=x_avg,
            y_avg=y_avg,
            device_name=device_name,
            prefix_name=prefix_name,
            exp_name=exp_name,
            pre_ramping_required=pre_ramping_required,
            save_in_database=True,
            return_full_fit_data=False,
            reverse=False,
            data_avg_num=3,
            data_fit_side='left'
            ):

    Vg,G_vals=GVG_fun(start_vg=start_vg,
                stop_vg=stop_vg,
                tc=tc,
                source_amplitude_instrumentlevel_GVg=source_amplitude_instrumentlevel_GVg,
                x_avg=x_avg,
                y_avg=y_avg,
                device_name=device_name,
                prefix_name=prefix_name,
                exp_name=exp_name,
                pre_ramping_required=pre_ramping_required,
                save_in_database=True,
                return_data=True,
                return_only_G=True,
                reverse=False
                )
    print("GVg done")
    if fit_type=='tunnel_broadened':
        slope,sitpos=fit_and_find_sitpos_singlepeak_tunnel(Vg,G_vals,initial_guess=initial_guess, sitfraction=sitfraction,return_full_fit_data=return_full_fit_data)
     
    if fit_type=='thermal':
        slope,sitpos=fit_and_find_sitpos_singlepeak_thermal(Vg,G_vals,initial_guess=initial_guess, sitfraction=sitfraction,return_full_fit_data=return_full_fit_data)
    
    if fit_type=='data':
        avg_G=centered_moving_average(G_vals,data_avg_num)
        max_avg=max(avg_G)
        if data_fit_side=='left':
            left_idx = np.argmax(avg_G > max_avg*sitfraction)
            sitpos=Vg[left_idx]
            x=[Vg[right_idx-2:right_idx+2]]
            y=[G_vals[right_idx-2:right_idx+2]]
            slope=scp.stats.linregress(x,y)

        if data_fit_side=='right':
            right_idx = len(avg_G) - 1 - np.argmax((avg_G[::-1] > x))
            sitpos=Vg[right_idx]
            x=[Vg[right_idx-data_avg_num:right_idx+data_avg_num]]
            y=[G_vals[right_idx-data_avg_num:right_idx+data_avg_num]]
            slope=scp.stats.linregress(x,y)
            
    
    gate.ramp_ch(sitpos) 


   
    return sitpos