

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
import copy
import time
from tqdm import tqdm
from experiment_functions.CS_functions import *
import matplotlib.pyplot as plt

#------User input----------------
run=False
#run=True
#adjustable hardware params

tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
attn_dB_source = 39 # attenuation at the source in dB
attn_dB_gate = 46+20
source_amplitude_instrumentlevel_GVg = 20e-3
mix_down_f = 1.25e6 # RLC frequency

#channel assignment
gate=qdac.ch06
freq = zurich.oscs.oscs0.freq
source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
measured_parameter = zurich.demods.demods0.sample 

#compensation - now defined globally in experiment_parameters
#x_avg=+3.4e-6  #+1.51e-5@75#+4.38e-6#@20mVpk -2.41e-5@100
#y_avg=-5.4e-6  #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100


#gate sweep params
start_vg = -1.87
stop_vg = -1.86
step_num= 100*100



pre_ramping_required=True

#costum name
device_name = 'CD11_D7_C1'
prefix_name = 'Conductance_rf_fittingandsitpos'
exp_name=prefix_name+device_name

#params
fit_type='data'
sitfraction='l_max_slope'
data_avg_num=5
min_acceptable_peak=50e-9
#fit_type='thermal'
#sitfraction=0.2



#fixed hardware params
#####################
gain_RT = 200       
gain_HEMT = 5.64   
Z_tot = 7521       
###################

def do_GVg_and_adjust_sitpos(
            fit_type=fit_type,
            initial_guess=None, 
            sitfraction=sitfraction,
            start_vg=start_vg,
            stop_vg=stop_vg,
            step_num=step_num,
            tc=tc,
            source_amplitude_instrumentlevel_GVg=source_amplitude_instrumentlevel_GVg,
            #x_avg=x_avg,
            #y_avg=y_avg,
            device_name=device_name,
            prefix_name=prefix_name,
            exp_name=exp_name,
            pre_ramping_required=pre_ramping_required,
            save_in_database=True,
            return_full_data=False,
            return_data=False,
            reverse=False,
            data_avg_num=data_avg_num,
            #data_fit_side='left',
            gate=gate
            ):

    Vg,G_vals=GVG_fun(start_vg=start_vg,
                stop_vg=stop_vg,
                step_num=step_num,
                tc=tc,
                source_amplitude_instrumentlevel_GVg=source_amplitude_instrumentlevel_GVg,
                #x_avg=x_avg,
                #y_avg=y_avg,
                device_name=device_name,
                prefix_name=prefix_name,
                exp_name=exp_name,
                pre_ramping_required=pre_ramping_required,
                save_in_database=False,
                return_data=True,
                return_only_Vg_and_G=True,
                reverse=False,
                gate=gate
                )
    #metadata
    vars_to_save=[tc,attn_dB_source,source_amplitude_instrumentlevel_GVg,x_avg,y_avg,start_vg,stop_vg,step_num,fit_type,sitfraction,data_avg_num,min_acceptable_peak]
    names_of_vars_to_save="tc,attn_dB_source,source_amplitude_instrumentlevel_GVg,x_avg,y_avg,start_vg,stop_vg,step_num,fit_type,sitfraction,data_avg_num,min_acceptable_peak"

    print("GVg done")
    if max(G_vals)<min_acceptable_peak:
        raise ValueError(f"maximum conductance lower than {min_acceptable_peak} S")

    if fit_type=='tunnel_broadened':
        popt, pcov,slope,sitpos=fit_and_find_sitpos_singlepeak_tunnel(Vg,G_vals,initial_guess=initial_guess, sitfraction=sitfraction,return_full_fit_data=True)
        fit_vals=breit_wigner_fkt(Vg,popt[0],popt[1],popt[2],popt[3])
    if fit_type=='thermal':
        popt, pcov,slope,sitpos=fit_and_find_sitpos_singlepeak_thermal(Vg,G_vals,initial_guess=initial_guess, sitfraction=sitfraction,return_full_fit_data=True)
        fit_vals=thermal_CB_peak(Vg,popt[0],popt[1],popt[2])
    if fit_type=='data':
        popt,pcov=None,None
        avg_G=centered_moving_average(G_vals,data_avg_num)
        fit_vals=avg_G
        max_avg=max(avg_G)
        deriv_avg=avg_G[:-1] - avg_G[1:]

        if isinstance(sitfraction, (int, float)):
            left_idx = np.argmax(avg_G > max_avg*sitfraction)
            sitpos=Vg[left_idx]
            x=[Vg[left_idx-data_avg_num:left_idx+data_avg_num]]
            y=[G_vals[left_idx-data_avg_num:left_idx+data_avg_num]]
            result=scp.stats.linregress(x,y)
            slope=result.slope

        elif sitfraction=="r_max_slope":
            rmax_id=np.argmax(deriv_avg)
            sitpos=(Vg[rmax_id]+Vg[rmax_id+1])/2
            slope=deriv_avg[rmax_id]/(Vg[rmax_id-1]-Vg[rmax_id])
        elif sitfraction=="l_max_slope":
            lmax_id=np.argmin(deriv_avg)
            sitpos=(Vg[lmax_id]+Vg[lmax_id+1])/2
            slope=deriv_avg[lmax_id]/(Vg[lmax_id-1]-Vg[lmax_id])
        elif sitfraction=="max":
            max_id=np.argmax(avg_G)
            sitpos=Vg[max_id]
            slope=0
        else:
            raise ValueError("sitpos must be a string or a number")
            

    gate.ramp_ch(sitpos) 

    if save_in_database:
        vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_num)
        if save_in_database:
            experiment = new_experiment(name=exp_name, sample_name=device_name)
            meas = Measurement(exp=experiment)
            meas.register_parameter(vgdc_sweep.parameter)
            meas.register_custom_parameter('G', unit='S', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('fit', unit='S', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('sitpos', unit='S', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('slope', unit='S', setpoints=[vgdc_sweep.parameter])
            with meas.run() as datasaver:
                qdac.add_dc_voltages_to_metadata(datasaver=datasaver)
                zurich.save_config_to_metadata(datasaver=datasaver)
                var_names=names_of_vars_to_save.split(',')
                for varname,var in zip(var_names,vars_to_save):
                    datasaver.dataset.add_metadata(varname,var)
                # Find the index of the value in Vg closest to sitpos
    
                approx_sitpos_index = np.argmin(np.abs(Vg - sitpos))

                if approx_sitpos_index in {0, len(Vg)-1}:
                    raise ValueError("sitpos is at beginning or end of sweep")


                # Define the approx_sitpos_array
                approx_sitpos_array = np.zeros_like(Vg)
                approx_sitpos_array[approx_sitpos_index] = G_vals[approx_sitpos_index]

                # Define the slope_array
                slope_array = np.zeros_like(Vg)
                slope_array[approx_sitpos_index] = fit_vals[approx_sitpos_index]
                #print(f"slope{slope}")
                #print(f"Vg[approx_sitpos_index]{Vg[approx_sitpos_index]}")
                #print(f"Vg[approx_sitpos_index-1]{Vg[approx_sitpos_index-1]}")
                
                slope_array[approx_sitpos_index-1] = fit_vals[approx_sitpos_index] - slope * (Vg[approx_sitpos_index] - Vg[approx_sitpos_index-1])
                slope_array[approx_sitpos_index+1] = fit_vals[approx_sitpos_index] + slope * (Vg[approx_sitpos_index+1] - Vg[approx_sitpos_index])
                #print(f"test:{slope_array[approx_sitpos_index]}")
                #plt.plot(Vg,slope_array)
                #plt.show()
                datasaver.add_result(('G', G_vals), ('fit',fit_vals),('sitpos',approx_sitpos_array),('slope',slope_array), (vgdc_sweep.parameter, Vg))
                datasaver.dataset.add_metadata("slope",slope)
                datasaver.dataset.add_metadata("sitpos",sitpos)
            if return_full_data:
                return Vg,G_vals,popt, pcov,slope,sitpos
            else:
                return slope,sitpos
if run:
    slope,sitpos=do_GVg_and_adjust_sitpos()