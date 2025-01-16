

import numpy as np

from instruments import station, zurich,qdac,Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from experiments.Do_GVg_and_adjust_sitpos import do_GVg_and_adjust_sitpos

from utils.d2v import d2v
from utils.v2d import v2d

import time
from tqdm import tqdm
import matplotlib.pyplot as plt
from utils.CS_utils import *
import copy


#------User input----------------
#costum name
device_name = 'CD11_D7_c1'
prefix_name = 'chargesensing_mechanics1dotsqueeze'



#adjustable hardware params
tc = 100e-3   # in seconds. Doesn't get overwritten by ZI called value.
vsd_dB = 39 # attenuation at the source in dB
mix_down_f = 1.25e6 # RLC frequency
source_amplitude_instrumentlevel_GVg = 20e-3

#channel assignment
source_amplitude_param = zurich.output0_amp0
gate_amplitude_param = zurich.output1_amp1
freq_mech = zurich.oscs.oscs1.freq
freq_rf = zurich.oscs.oscs0.freq
freq_rlc = zurich.oscs.oscs2.freq
gate=qdac.ch06
measured_parameter = zurich.demods.demods2.sample #for mechanics

#frequency sweep params
start_f = 159e6#162.62e6 #Hz unit
stop_f =  162e6 #Hz unit
step_num_f = 1000*3#

freq_sweep_avg_nr=9

#####################

#gate sweep params
#gate sweep params
start_vg = -1.1555
stop_vg = -1.1525
step_num= 3*30



#GVg fit params
fit_type='data'
sitfraction="r_max_slope"
data_avg_num=7

switch_off_gate_drive_for_GVg=True




#calculate derived quantities
vsdac=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_GVg)-vsd_dB)/10 #rf amplitude at source

#create postfix, labels, and other names
postfix = f"_{round(gate_amplitude_param()*1000,3)}mV on gate@inst,_{round(source_amplitude_param()*1000,3)}mV on source@inst, g1={round(qdac.ch01.dc_constant_V(),2)},g2={round(qdac.ch02.dc_constant_V(),5)},g3={round(qdac.ch03.dc_constant_V(),2)},g4={round(qdac.ch04.dc_constant_V(),5)},g5={round(qdac.ch05.dc_constant_V(),2)},gcs={round(qdac.ch06.dc_constant_V(),5)}"
#exp_dict = dict(vsdac = vsdac)
exp_name = prefix_name+device_name+postfix #sample_name(prefix_name,exp_dict,postfix)
freq_rlc(mix_down_f)
freq_mech(start_f)
freq_rf(start_f-mix_down_f)
time.sleep(1)

#define frequency sweep
freq_sweep = freq_rf.sweep(start=start_f, stop=stop_f, num = step_num_f)
freq_sweep_list=list(freq_sweep)
 
print("preramping")
qdac.ramp_multi_ch_slowly(channels=[gate], final_vgs=[start_vg])


# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(freq_sweep.parameter)  # 
meas.register_custom_parameter('I_rf', 'current', unit='I', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('I_rf_avg', 'current_avg', unit='I', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[freq_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[freq_sweep.parameter])

# # -----------------Start the Measurement-----------------------

with meas.run() as datasaver:

    qdac.add_dc_voltages_to_metadata(datasaver=datasaver)
    zurich.save_config_to_metadata(datasaver=datasaver)
    
    if switch_off_gate_drive_for_GVg:
        zurich.sigout1_amp1_enabled_param.value(0)
    slope_first,sitpos_first=do_GVg_and_adjust_sitpos(start_vg=start_vg,
                             stop_vg=stop_vg,
                             step_num=step_num,
                             fit_type=fit_type,
                             sitfraction=sitfraction,
                             data_avg_num=data_avg_num,
                             gate=gate,
                             save_in_database=True,
                             pre_ramping_required=False
                             )
    if switch_off_gate_drive_for_GVg:
        zurich.sigout1_amp1_enabled_param.value(1)

    print(f"I've just set the gate to {qdac.ch06.dc_constant_V()}")
    theta_calc, v_r_calc, I, G = theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
    first_sit_G=copy.copy(G)
    print(f"initial conductance is {first_sit_G:.4g}")
    I_list=[]
    for f_value in tqdm(freq_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
        freq_rf(f_value-freq_rlc())
        freq_mech(f_value)
        time.sleep(1.1*tc) # Wait 1.1 times the time contanst of the lock-in
        measured_value=measured_parameter()#now measuring demod2 ie non-standard
        theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac=vsdac,measured_value=measured_value)
                
        #G calculation
        I_list.append(I)
    
        datasaver.add_result(('I_rf', I),
                            ('V_r', v_r_calc),
                            ('Phase', theta_calc),
                            (freq_sweep.parameter,f_value))
    #averaged values
    I_avg=centered_moving_average(I_list,n=freq_sweep_avg_nr)

    datasaver.add_result(('I_rf_avg', I_avg),(freq_sweep.parameter,freq_sweep_list))#try this first
    
    #final check:
    #measured_value=measured_parameter()
    freq_rf(mix_down_f)
    time.sleep(0.1)
    theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
    end_sit_G=copy.copy(G)
    print(f"final conductance is {end_sit_G:.4g}")
    if switch_off_gate_drive_for_GVg:
        zurich.sigout1_amp1_enabled_param.value(0)
    slope_last,sitpos_last=do_GVg_and_adjust_sitpos(start_vg=start_vg,
                             stop_vg=stop_vg,
                             step_num=step_num,
                             fit_type=fit_type,
                             sitfraction=sitfraction,
                             data_avg_num=data_avg_num,
                             gate=gate,
                             save_in_database=True,
                             pre_ramping_required=False
                             )
    print(f"I've in the end set the gate to {qdac.ch06.dc_constant_V()}")
    G_delta=end_sit_G-first_sit_G
    sitpos_delta=sitpos_last-sitpos_first
    slope_delta=slope_last-slope_first
    G_delta_on_G=G_delta/end_sit_G
    slope_delta_on_slope=slope_delta/slope_last
    print(f"G_delta={G_delta}")
    print(f"sitpos_delta={sitpos_delta}")
    print(f"slope_delta={slope_delta}")
    print(f"G_delta_on_G={G_delta_on_G*100} %")
    print(f"slope_delta_on_slope={slope_delta_on_slope*100} %")

    vars_to_save_postrun=[first_sit_G,end_sit_G,sitpos_last,sitpos_first,slope_last,slope_first,G_delta,sitpos_delta,slope_delta,slope_delta_on_slope,G_delta_on_G]
    names_of_vars_to_save_postrun = "first_sit_G,end_sit_G,sitpos_last,sitpos_first,slope_last,slope_first,G_delta,sitpos_delta,slope_delta,slope_delta_on_slope,G_delta_on_G"
    var_names=names_of_vars_to_save_postrun.split(',')
    for varname,var in zip(var_names,vars_to_save_postrun):
        datasaver.dataset.add_metadata(varname,var)

"""    

def cs_mechanics_simple_setpoint(start_f,stop_f,step_num_f,
                                start_vg,stop_vg,step_num,
                                fit_type='data',data_avg_num=3,sitfraction="l_max_slope",freq_sweep_avg_nr=freq_sweep_avg_nr,
                                check_at_end=False,return_GVgs=False,return_all_fit_data=False,switch_off_gate_drive_for_GVg=False):
    if switch_off_gate_drive_for_GVg:
        zurich.sigout1_amp1_enabled_param.value(0)
    Vg_before,G_vals_before,popt_before, pcov_before,slope_before,sitpos_before=do_GVg_and_adjust_sitpos(start_vg=start_vg,
                                                                                                        stop_vg=stop_vg,
                                                                                                        step_num=step_num,
                                                                                                        fit_type=fit_type,
                                                                                                        sitfraction=sitfraction,
                                                                                                        data_avg_num=data_avg_num,
                                                                                                        gate=gate,
                                                                                                        save_in_database=False,
                                                                                                        return_full_data=True,
                                                                                                        pre_ramping_required=False
                                                                                                        )
    if switch_off_gate_drive_for_GVg:
        zurich.sigout1_amp1_enabled_param.value(1)
    Vg_before_freq_sweep=qdac.ch06.dc_constant_V()
    print(f"I've just set the gate to {Vg_before_freq_sweep}")

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
        Vg_after,G_vals_after,popt_after, pcov_after,slope_after,sitpos_after=do_GVg_and_adjust_sitpos(start_vg=start_vg,
                                                                                                    stop_vg=stop_vg,
                                                                                                    step_num=step_num,
                                                                                                    fit_type=fit_type,
                                                                                                    sitfraction=sitfraction,
                                                                                                    data_avg_num=data_avg_num,
                                                                                                    gate=gate,
                                                                                                    pre_ramping_required=False
                                                                                                    )
        Vg_after_freq_sweep=qdac.ch06.dc_constant_V()
        print(f"I've just readjusted the gate by {Vg_after_freq_sweep-Vg_before_freq_sweep}")

        if return_GVgs:
            values_to_return.update({"Vg_after": Vg_after, "G_vals_after": G_vals_after})
            if return_all_fit_data:
                    values_to_return.update({"popt_after": popt_after, "pcov_after": pcov_after, "slope_after": slope_after, "sitpos_after": sitpos_after})
        slope_estimate=(slope_before+slope_after)/2
    else:
        slope_estimate=slope_before
    
    values_to_return.update({"I_slope_normalized" : I_np/slope_estimate})

    return values_to_return
"""