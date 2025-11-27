

import numpy as np

from instruments import station, zurich,qdac,Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from instruments import exp
#from experiments.Do_GVg_and_adjust_sitpos import do_GVg_and_adjust_sitpos

from utils.d2v import d2v
from utils.v2d import v2d

import time
from tqdm import tqdm
import matplotlib.pyplot as plt
from utils.CS_utils import *
import copy


#------User input----------------


#channel assignment
source_amplitude_param = zurich.output0_amp0
gate_amplitude_param = zurich.output1_amp1
freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq
gate=qdac.ch06
measured_parameter = zurich.demods.demods2.sample #for mechanics












def cs_mechanics_simple_setpoint(start_f,stop_f,step_num_f,
                                start_vg,stop_vg,step_num,
                                fit_type='data',data_avg_num=3,sitfraction="l_max_slope",freq_sweep_avg_nr=10,
                                check_at_end=False,return_GVgs=False,return_all_fit_data=False,switch_off_gate_drive_for_GVg=False,tc=100e-3,vsdac=15e-6):
    if switch_off_gate_drive_for_GVg:
        zurich.sigout1_amp1_enabled_param.value(0)
    #Vg,G_vals,popt, pcov,slope,sitpos
    Vg_before,slope_before=exp.sit_at_max_Isens(avg_num=3,return_sitpos_and_sens=True,side=None,start_vg=None,stop_vg=None,step_num=None)

    if switch_off_gate_drive_for_GVg:
        zurich.sigout1_amp1_enabled_param.value(1)
    Vg_before_freq_sweep=qdac.ch06.dc_constant_V()
    #print(f"I've just set the gate to {Vg_before_freq_sweep}")

    freq_sweep = freq_rf.sweep(start=start_f, stop=stop_f, num = step_num_f)
    freq_sweep_list=list(freq_sweep)

    Vr_list,I_list,Phaselist=[],[],[]
    for f_value in tqdm(freq_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
        freq_rf(f_value-freq_rlc())
        freq_mech(f_value)
        time.sleep(1.1*tc) # Wait 1.1 times the time contanst of the lock-in
        measured_value=measured_parameter()#now measuring demod2 ie non-standard
        theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac=vsdac,measured_value=measured_value)

        Vr_list.append(v_r_calc)        
        I_list.append(I)
        Phaselist.append(theta_calc)

    I_avg=centered_moving_average(I_list,n=freq_sweep_avg_nr)
    freq_np=np.array(freq_sweep_list)
    Vr_np,I_np,Phase_np,I_avg_np=np.array(Vr_list),np.array(I_list),np.array(Phaselist),np.array(I_avg)

    values_to_return={"freq" : freq_np,"V": Vr_np, "I": I_np, "Phase": Phase_np, "I_avg" : I_avg_np}
    

    if return_GVgs:
        values_to_return.update({"Vg_before" : Vg_before, "G_vals_before" :G_vals_before}) 
        if return_all_fit_data:
            values_to_return.update({"popt_before": popt_before,"pcov_before": pcov_before,"slope_before": slope_before,"sitpos_before": sitpos_before})

    if check_at_end:
        if switch_off_gate_drive_for_GVg:
            zurich.sigout1_amp1_enabled_param.value(0)
        Vg_after,G_vals_after,popt_after, pcov_after,slope_after,sitpos_after=exp.do_GVg_and_adjust_sitpos(start_vg=start_vg,
                                                                                                    stop_vg=stop_vg,
                                                                                                    step_num=step_num,
                                                                                                    fit_type=fit_type,
                                                                                                    sitfraction=sitfraction,
                                                                                                    data_avg_num=data_avg_num,
                                                                                                    gate=gate,
                                                                                                    pre_ramping_required=False
                                                                                                    )
        Vg_after_freq_sweep=qdac.ch06.dc_constant_V()
        #print(f"I've just readjusted the gate by {Vg_after_freq_sweep-Vg_before_freq_sweep}")

        if return_GVgs:
            values_to_return.update({"Vg_after": Vg_after, "G_vals_after": G_vals_after})
            if return_all_fit_data:
                    values_to_return.update({"popt_after": popt_after, "pcov_after": pcov_after, "slope_after": slope_after, "sitpos_after": sitpos_after})
        slope_estimate=(slope_before+slope_after)/2
    else:
        slope_estimate=slope_before
    
    values_to_return.update({"I_slope_normalized" : I_np/slope_estimate})

    return values_to_return