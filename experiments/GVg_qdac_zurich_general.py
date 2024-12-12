

import numpy as np


from instruments import station, zurich, qdac,Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
from utils.CS_utils import save_metadata_var, get_var_name
import experiment_parameters

import time
from tqdm import tqdm


#------User input----------------
run=False
run=True
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
x_avg=experiment_parameters.x_avg#+3.4e-6  #+1.51e-5@75#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=experiment_parameters.y_avg#-5.4e-6  #-1.75e-5#@75-4.41e-6#@20mVpk -6.14e-5@100


#gate sweep params
start_vg = -0.8345
stop_vg = -0.8355
step_num= 1*100

#for metadata
vars_to_save=[tc,vsd_dB,source_amplitude_instrumentlevel_GVg,x_avg,y_avg,start_vg,stop_vg,step_num]

pre_ramping_required=True

#costum name
device_name = 'CD11_D7_C1'
prefix_name = 'Conductance_rf_'
exp_name=f"vsac@inst={source_amplitude_instrumentlevel_GVg*1e3} mV"

#fixed hardware params
#####################
gain_RT = 200       
gain_HEMT = 5.64   
Z_tot = 7521       
###################


#define function to call later
def GVG_fun(start_vg=start_vg,
            stop_vg=stop_vg,
            step_num=step_num,
            tc=tc,
            source_amplitude_instrumentlevel_GVg=source_amplitude_instrumentlevel_GVg,
            x_avg=x_avg,
            y_avg=y_avg,
            device_name=device_name,
            prefix_name=prefix_name,
            exp_name=exp_name,
            pre_ramping_required=pre_ramping_required,
            save_in_database=True,
            return_data=False,
            return_only_Vg_and_G=True,
            reverse=False,
            gate=gate
            ):
    #calculate derived quantities
    step_vg=np.absolute((start_vg-stop_vg)/step_num) #gate step size
    vsdac=d2v(v2d(np.sqrt(1/2)*source_amplitude_instrumentlevel_GVg)-vsd_dB)/10 #rf amplitude at source
    #vars_to_save.extend([step_vg,vsdac])
    #print(f"source amp at CNT for GVg:{vsdac*1e6} uV")

    #define sweep
    vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_num)
    
    #create postfix, labels, and other names
    #Temp=Triton.MC()
    if save_in_database:
        postfix = f"_g1={round(qdac.ch01.dc_constant_V(),2)},g2={round(qdac.ch02.dc_constant_V(),2)},g3={round(qdac.ch03.dc_constant_V(),2)},g4={round(qdac.ch04.dc_constant_V(),2)},g5={round(qdac.ch05.dc_constant_V(),2)}"
        #postfix = f"_5g={round(qdac.ch05.dc_constant_V(),7)}"
        #postfix = f"_{zurich.sigouts.sigouts0.amplitudes.amplitudes0.value()}mVrf_pk"
        gate.label = 'cs_gate' # Change the label of the gate chaneel
        #instr_dict = dict(gate=[gate])
        exp_dict = dict(mV = vsdac*1000)

        exp_name = sample_name(prefix_name,exp_dict,postfix)

    #select variables to be saved in metadata
    vars_to_save=[tc,vsd_dB,source_amplitude_instrumentlevel_GVg,vsdac,x_avg,y_avg]

    #------------init--------------------
    #print("starting init")
    freq(mix_down_f)
    source_amplitude_param(source_amplitude_instrumentlevel_GVg)

    if pre_ramping_required:
        #print("preramping")
        qdac.ramp_multi_ch_slowly(channels=[gate], final_vgs=[start_vg])

    gate.ramp_ch(start_vg)
   # print(f"init done, gate at {gate.dc_constant_V()}")

    
    Glist, Vlist, Ilist, Phaselist, Rlist = [], [], [], [], []

    if save_in_database:
        experiment = new_experiment(name=exp_name, sample_name=device_name)
        meas = Measurement(exp=experiment)
        meas.register_parameter(vgdc_sweep.parameter)
        meas.register_custom_parameter('G', unit='S', setpoints=[vgdc_sweep.parameter])
        meas.register_custom_parameter('V_r', unit='V', setpoints=[vgdc_sweep.parameter])
        meas.register_custom_parameter('Phase', unit='rad', setpoints=[vgdc_sweep.parameter])
        meas.register_custom_parameter('R', unit='Ohm', setpoints=[vgdc_sweep.parameter])
        meas.register_custom_parameter('I', unit='A', setpoints=[vgdc_sweep.parameter])

        with meas.run() as datasaver: 
            qdac.add_dc_voltages_to_metadata(datasaver=datasaver)
            varnames = [get_var_name(var) for var in vars_to_save]
            save_metadata_var(datasaver.dataset, varnames, vars_to_save)

            

            for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='gate voltage Sweep', colour='green'):
                gate.ramp_ch(vgdc_value)
                time.sleep(1.1 * tc)  
                measured_value = measured_parameter()#fix the following line according to driver. push driver
                theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                R = 1 / G
                
                datasaver.add_result(('R', R), ('G', G), ('V_r', v_r_calc), ('Phase', theta_calc),
                                     ('I', I), (vgdc_sweep.parameter, vgdc_value))

                if return_data:
                    Glist.append(G)
                    Vlist.append(v_r_calc)
                    Ilist.append(I)
                    Phaselist.append(theta_calc)
                    Rlist.append(R)

    if save_in_database==False:
        for vgdc_value in tqdm(vgdc_sweep, leave=False, desc='gate voltage Sweep', colour='green'):
                gate.ramp_ch(vgdc_value)
                time.sleep(1.1 * tc)  
                measured_value = measured_parameter()#fix the following line according to driver. push driver
                theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                R = 1 / G
                
                Glist.append(G)
                Vlist.append(v_r_calc)
                Ilist.append(I)
                Phaselist.append(theta_calc)
                Rlist.append(R)

    if return_data and reverse:
            Glist.reverse()
            Vlist.reverse()
            Ilist.reverse()
            Phaselist.reverse()
            Rlist.reverse()
        
    Vglist =list(vgdc_sweep)

    if return_data:
            if return_only_Vg_and_G:
                return np.array(Vglist),np.array(Glist)
            else:
                return np.array(Vglist),np.array(Glist), np.array(Vlist), np.array(Ilist), np.array(Phaselist), np.array(Rlist)

# Call function
if run:
    GVG_fun()