# standard charge-sensing stability diagram procedure
#suggested name: DQD_charge_stability_viaCS_standard
import numpy as np
import copy

from instruments import  manual, station, qdac,zurich, Triton
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name


from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
import time
from tqdm import tqdm
import scipy as scp
import matplotlib.pyplot as plt
from utils.CS_utils import breit_wigner_fkt,  breit_wigner_detuning, save_metadata_var, get_var_name
#from drivers.Qdac_utils import ramp_QDAC_multi_channel
from qcodes import Parameter
from instruments import exp
#------User input----------------

ramp_speed = 0.01 # V/s for large ramps
step_ramp_speed=0.01 # between steps, V/s
tc = exp.tc#30e-3   # in seconds. 
vsdac =  15.8e-6 # source DC voltage in volt
att_source_dB = 39 # attenuation at the source in dB
att_gate_dB =46 
device_name = exp.device_name
prefix_name = 'charge_stability_'
#postfix = '20mK_'
#offset = -10e-6 #voltage offset of k2400
#offset_i=-44e-12



debug=False
x_avg=exp.x_avg#+1.24465881e-06#+4.38e-6#@20mVpk -2.41e-5@100
y_avg=exp.y_avg#-1.07161223e-06


#outer voltage range (slow axis2)
#####################
start_vg1 = -0.5
stop_vg1 = +0.5
step_vg1_num =50
step_vg1=np.absolute((start_vg1-stop_vg1)/step_vg1_num)

vars_to_save=[ramp_speed,step_ramp_speed,tc,att_source_dB,att_gate_dB,debug,x_avg,y_avg,step_vg1]#more to add later



#inner voltage range (fast axis)
#####################
start_vg2 = -0.5
stop_vg2 =  +0.5
#stop_vg2 =  -1.571#-1.875#delta=10mV
step_vg2_num=1000
step_vg2=np.absolute((start_vg2-stop_vg2)/step_vg2_num)
vars_to_save.append(step_vg2)
time.sleep(10)


step_cs_num=500*1#10uV
delta=5e-3#10mV
vars_to_save.extend([step_cs_num,delta])

sitfraction=0.55# dhow far up the peak
lower_G_bound_fraction=0.6# no big problem if too low
upper_G_bound_fraction=1.3#not too high to make sure we dont fall over peak

upper_noise_bound=20e-9#Siemens, lowest permissible value of measured G that's not considered noise
lower_peak_bound=50e-9#Siemens, lowest value of peak conductance that allows it to be considered a peak
vars_to_save.extend([sitfraction,lower_G_bound_fraction,upper_G_bound_fraction])
######################ramping gates

qdac.ramp_multi_ch_slowly([1,2,3,4,5,6],[0.8,-0.5,0.3,-0.5,0.5,0.7])



V_GVg,G_GVg=exp.GVG_fun(start_vg=0.7,
            stop_vg=0.9,
            step_num=200*10,
            pre_ramping_required=True,
            save_in_database=True,
            return_data=True,
            return_only_Vg_and_G=True,
            reverse=False
            )

start_vgcs=V_GVg[np.argmax(G_GVg)]

print(f"automatically chosen highest peak at {start_vgcs}, max conductance is {max(G_GVg)*1e6} uS")

#start_vgcs=-1.2195 #-0lowerV slope, 140nS

#
crosscap_outer_gate=-0.18
crosscap_inner_gate=-0.03

Run_GVg_for_each_outer_value=True

vars_to_save.extend([start_vgcs,crosscap_outer_gate,crosscap_inner_gate,Run_GVg_for_each_outer_value])

#swept contacts
gate1=qdac.ch02
  # swept outer gate voltage
gate2=qdac.ch04 #swept inner gate voltage
#source = k2400 # source 
csgate=qdac.ch06

def GVg(start,stop,num,csgate=csgate):
    Glistcs=[]
    cs_sweep=csgate.dc_constant_V.sweep(start=start, stop=stop, num = num)
    csgate.dc_constant_V(GVg_startpos-delta_reduced)

                        ###########131025 first step in cleaning up this code: replace the whole next block with a d-_GVG_and....block
    time.sleep(abs(GVg_startpos-delta_reduced-current_csvg)/step_ramp_speed)
    for gatecs_value in (cs_sweep):
            csweeplist=list(cs_sweep)
            cs_sweep.set(gatecs_value)
            time.sleep(1.1*tc+2*delta_reduced/step_num_reduced/step_ramp_speed)
            measured_value = measured_parameter()
            _, _, _, G = zurich.phase_voltage_current_conductance_compensate(vsdac)

            Glistcs.append(G)
    return np.array(csweeplist),np.array(Glistcs)


print("bias")
print(qdac.ch07.dc_constant_V())


#--------Definitions-------------



postfix =f"g1={qdac.ch01.dc_constant_V():.4g},g3={qdac.ch03.dc_constant_V():.4g},g5={qdac.ch05.dc_constant_V():.4g},gcsstart={start_vgcs:.6g}"

#gate1.label = 'gate2' # 
#gate2.label = 'gate4' # 
#instr_dict = dict(gate1=[gate1])
#exp_dict = dict(mV = vsdac*1000)
exp_name = prefix_name+device_name#sample_name(prefix_name,exp_dict,postfix)

#----------- defined values------#----------- defined values------

###################
freq_rlc = zurich.oscs.oscs0.freq
freq_rlc(exp.mix_down_f)

# ------------------define sweep axes-------------------------

gate1_sweep=gate1.dc_constant_V.sweep(start=start_vg1, stop=stop_vg1, num = step_vg1_num)
gate2_sweep=gate2.dc_constant_V.sweep(start=start_vg2, stop=stop_vg2, num = step_vg2_num)
g1sweeplist=list(gate1_sweep)
g2sweeplist=list(gate2_sweep)




current_csvg=start_vgcs



sleeptime=100#max([abs((gate_V_ch1-qdac.ch01.dc_constant_V())/ramp_speed),abs((gate_V_ch3-qdac.ch03.dc_constant_V())/ramp_speed),abs((gate_V_ch5-qdac.ch05.dc_constant_V())/ramp_speed),abs((start_vg1-gate1.dc_constant_V())/ramp_speed),abs((start_vg2-gate2.dc_constant_V())/ramp_speed),abs((start_vgcs-csgate.dc_constant_V())/ramp_speed)])+30  #wait for the time it takes to do both ramps plus one second
print(f"sleeptime={sleeptime}")
time.sleep(sleeptime)
measured_parameter = zurich.demods.demods0.sample
vsdac0 = rms2pk(d2v(v2d(vsdac)+att_source_dB))  

print("gate channels")
print(qdac.ch01.dc_constant_V())
print(qdac.ch02.dc_constant_V())
print(qdac.ch03.dc_constant_V())
print(qdac.ch04.dc_constant_V())
print(qdac.ch05.dc_constant_V())
print(qdac.ch06.dc_constant_V())
print(qdac.ch07.dc_constant_V())
#set fast ramp speeds



# ----------------Create a measurement-------------------------
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(gate1_sweep.parameter)  # 
meas.register_parameter(gate2_sweep.parameter)  # 
meas.register_custom_parameter('signal_shift_Vx_deriv', 'signal_shift_Vx_deriv', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('signal_shift_V', 'V_sig', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('peak_position', 'V_peak', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('peak_fit_position', 'V_peak_fit', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])

meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
#meas.register_custom_parameter('DeltaG', 'DeltaG', unit='S', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])

meas.register_custom_parameter('peak_Value', 'G_peak', unit='S', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
#for visualizing signal

meas.register_custom_parameter('signal_shift_V_deriv', 'V_sig_deriv', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('signal_shift_Vxn_deriv', 'signal_shift_Vxn_deriv', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])

meas.register_custom_parameter('signal_shift_cbrt_V', 'V_sig_cbrt', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('signal_shift_2pt_deriv', 'signal_shift_2pt_deriv', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])

#for tests of code and debugging
meas.register_custom_parameter('V_detune', 'V_detune', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('V_corr', 'V_corr', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('V_corr_imp', 'V_corr_imp', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('V_corr_pos', 'V_corr_pos', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])
meas.register_custom_parameter('fit_detuning', 'fit_detuning', unit='Volt', basis=[], setpoints=[gate1_sweep.parameter,gate2_sweep.parameter])


GVg_num_param = Parameter('GVg_num_param',
                            label='GVg_num',
                           # unit='V',  # Or the appropriate unit
                            set_cmd=None,  # If there is no instrument to communicate with directly
                            get_cmd=None)  # Define get_cmd if you need to read a value

CS_V_param = Parameter('CS_V_param',
                            label='CS_V',
                           # unit='V',  # Or the appropriate unit
                            set_cmd=None,  # If there is no instrument to communicate with directly
                            get_cmd=None)  # Define get_cmd if you need to read a value

GVg_num_param.set(0)
experiment_aux1 = new_experiment(name=exp_name+'_aux1', sample_name=device_name)
meas_aux1 = Measurement(exp=experiment_aux1)
meas_aux1.register_parameter(GVg_num_param)  # 
meas_aux1.register_parameter(CS_V_param)  #
#meas_aux1.register_parameter(gate1_sweep.parameter)  # 
#meas_aux1.register_parameter(gate2_sweep.parameter)  # 
meas_aux1.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[GVg_num_param,CS_V_param])
meas_aux1.register_custom_parameter('Vg1', 'Vg1', unit='V', basis=[], setpoints=[GVg_num_param,CS_V_param])
meas_aux1.register_custom_parameter('Vg2', 'Vg2', unit='V', basis=[], setpoints=[GVg_num_param,CS_V_param])






# # -----------------Start the Measurement-----------------------
 
# inverse_source_sweep=source_sweep.reverse() # or define function

with meas.run() as datasaver:
    with meas_aux1.run() as datasaver_aux1:
        #save barrier gates
        qdac.add_dc_voltages_to_metadata(datasaver=datasaver)
    # datasaver.dataset.add_metadata('qdac_ch06_dc_constant_V',qdac.ch06.dc_constant_V())
        varnames=[]
        for i in range(len(vars_to_save)):
            varnames.append(get_var_name(vars_to_save[i]))
        save_metadata_var(datasaver.dataset,varnames,vars_to_save)

        fast_axis_unreversible_list = list(gate2_sweep) #(to deal with snake)
        reversed_sweep=False#(to deal with snake)
        First_run=True#to not correct for cross-capacitance in first run, and to not refer to previously saved values in first run
        First_outer_run=True#to not calculate derivative in first run
        
        #initalize use values
        Glast=0
        last_GVg_G=0
        peak_G_fit=0
        overall_implemented_correction=0#total adjustment of CSgate voltage implented between gvgs
        G=0
        time_spent_on_simple_measurement=0#init
        time_spent_on_GVgs=0#init
        #for test/debugging
        do_GVg_for_test=False
        do_GVg_anyway=False
        GinRange=True
        nr_of_topslips=0
        nr_of_bottomslips=0
        nr_GVgs=0
        nr_simple_meas=0
        last_sitpos=start_vgcs
        peak_fit=copy.copy(start_vgcs)
        sweeplists=[]#cs sweeps voltages
        Glistscs=[]#cs sweeps conductances
        overall_start_time=time.time_ns()
        for gate1_value in tqdm(gate1_sweep, leave=False, desc='outer gate sweep', colour = 'green'): #slow axis loop (gate)

            #make one aux_meas for each outer run,maybe; or define one new variable for each, maybe - try that
            #if gate1_value>-1.44:
            #    do_GVg_for_test=True

            gate1_sweep.set(gate1_value)
            if First_run==False:
                current_csvg+=step_vg1*crosscap_outer_gate
                csgate.dc_constant_V(current_csvg)
            time.sleep(abs(step_vg1/step_ramp_speed))#assuming crosscapactivance<1 
            #test
        
            Glist=[]
            Vlist=[]
            Phaselist=[]
            #peakpos_list=[]
            peakpos_fit_list=[]
            peakG_list=[]
            #detune_V_list=[]
            Signal_V_list=[]
            Signal_V_deriv_list=[]
            Signal_V_cbrt_list=[]
            Vcorr_list=[]
            V_corr_imp_list=[]
            V_corr_pos_list=[]
            #sitposlist=[]
            GVg_peakposlist=[]
            GVg_sitposlist=[]
            fit_detuning_list=[]

            
    
            
            First_inner_run=True
            for gate2_value in tqdm(gate2_sweep, leave=False, desc='inner gate sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
                gate2_sweep.set(gate2_value)
                correctionV=(gate1_value-start_vg1)*crosscap_outer_gate+(gate2_value-start_vg2)*crosscap_inner_gate
                # Wait 3 times the time constant of the lock-in, plus the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
                if First_inner_run==False:
                    do_GVg_anyway=False
                    if reversed_sweep:
                        current_csvg-=step_vg2*crosscap_inner_gate
                        peak_fit-=step_vg2*crosscap_inner_gate
                        if start_vgcs-100e-3 < current_csvg < start_vgcs+100e-3:
                            csgate.dc_constant_V(current_csvg)
                    else:
                        current_csvg+=step_vg2*crosscap_inner_gate
                        peak_fit+=step_vg2*crosscap_inner_gate
                        if start_vgcs-100e-3 < current_csvg < start_vgcs+100e-3:
                            csgate.dc_constant_V(current_csvg)
                else:#ie it's the first inner run
                    if Run_GVg_for_each_outer_value:
                        do_GVg_anyway=True       
                start_time=time.time_ns()           
                Vcorr_list.append(correctionV)
                time.sleep(1.1*tc+step_vg2/step_ramp_speed)#crosscapacitance corrections are smaller
                measured_value = measured_parameter()
                theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                #print(f"G={G}")

                G_min=peak_G_fit*sitfraction*lower_G_bound_fraction
                G_max=peak_G_fit*sitfraction*upper_G_bound_fraction
                
                if First_run==False:
                    
                    if G<=peak_G_fit:
                        fit_detuning=-breit_wigner_detuning(G,peak_G_fit,hgamma_fit)#-for lower voltage side
                        GinRange=True
                    else:
                        #print("G>peak_G_fit")
                        GinRange=False

                    peak_fit=current_csvg-fit_detuning
                    fit_calculated_sitpos=peak_fit-breit_wigner_detuning(peak_G_fit*sitfraction,peak_G_fit,hgamma_fit)#estimated peak position-normal detuning
                    fit_position_correction=fit_calculated_sitpos-current_csvg

               
                    
                    
                if G_min<G<G_max and First_run==False and do_GVg_anyway==False and GinRange and G>upper_noise_bound: #no GVg
                    
                    nr_simple_meas+=1
                    deltaG=G-last_GVg_G
                    #peakpos=last_measured_peakpos-correctionV-deltaG/CS_peak_slope
                    #determine position correction, locked to lower V side
                    corrected_peakpos_fit=last_measured_corrected_peakpos+fit_position_correction+overall_implemented_correction
                    

                    if abs(fit_position_correction)>5e-6:#nearly always
                        #print("correcting position")
                        if start_vgcs-100e-3 < fit_calculated_sitpos < start_vgcs+100e-3:
                            csgate.dc_constant_V(fit_calculated_sitpos)
                        if First_run:
                            print(f"sleeptime first step{abs(fit_calculated_sitpos-current_csvg)/step_ramp_speed}")
                        time.sleep(abs(fit_calculated_sitpos-current_csvg)/step_ramp_speed)
                        current_csvg=fit_calculated_sitpos
                        overall_implemented_correction+=fit_position_correction#overall correction since last GVg
                        last_sitpos=current_csvg

                    #save data, since now it's clear there wont be a GVG
                    fit_detuning_list.append(fit_detuning)   
                    V_corr_imp_list.append(overall_implemented_correction)#
                    
                    Glist.append(G)#test new syntax
                    Vlist.append(v_r_calc)
                    Phaselist.append(theta_calc)
                    peakG_list.append(0)

                    Signal_V_list.append(corrected_peakpos_fit)
                    Signal_V_cbrt_list.append(np.cbrt(corrected_peakpos_fit))
                    peakpos_fit_list.append(peak_fit)
                    V_corr_pos_list.append(fit_position_correction)

                    First_run=False
                    
                    Glast=copy.copy(G)
                    
                    time_spent_on_simple_measurement=time_spent_on_simple_measurement+(time.time_ns()-start_time)/1e9
                else:#do GVg
                    

                    if First_run==False:
                        last_sitpos=copy.copy(current_csvg)
                    nr_GVgs+=1
                    if G>=G_max:
                        nr_of_topslips+=1
                    elif G<=G_min:
                        nr_of_bottomslips+=1
                    

                    start_time=time.time_ns()
                    
                    overall_implemented_correction=0
                    if start_vgcs-100e-3 < peak_fit < start_vgcs+100e-3:
                        GVg_startpos=copy.copy(peak_fit) 
                    for m in range(4):
                       
                        Glistcs=[]
                        p=3-m
                        delta_reduced=delta/2**p
                        step_num_reduced=step_cs_num/2**p
                        cs_sweep=csgate.dc_constant_V.sweep(start=GVg_startpos-delta_reduced, stop=GVg_startpos+delta_reduced, num = round(step_num_reduced))
                        csgate.dc_constant_V(GVg_startpos-delta_reduced)

                        ###########131025 first step in cleaning up this code: replace the whole next block with a d-_GVG_and....block
                        time.sleep(abs(GVg_startpos-delta_reduced-current_csvg)/step_ramp_speed)
                        for gatecs_value in (cs_sweep):
                            csweeplist=list(cs_sweep)
                            cs_sweep.set(gatecs_value)
                            time.sleep(1.1*tc+2*delta_reduced/step_num_reduced/step_ramp_speed)
                            measured_value = measured_parameter()
                            theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)

                            Glistcs.append(G)
                            #save all data in big array
                            #AllData3d_conductance[g1sweeplist.index(gate1_value),g2sweeplist.index(gate2_value),csweeplist.index(gatecs_value)]=G
                            #AllData3d_Vg[g1sweeplist.index(gate1_value),g2sweeplist.index(gate2_value),csweeplist.index(gatecs_value)]=gatecs_value
                        current_csvg=csweeplist[-1]
                        peak_G=max(Glistcs)
                        maxindex=Glistcs.index(max(Glistcs))
                        peakpos=csweeplist[maxindex]
                        if peak_G>lower_peak_bound and maxindex<(len(csweeplist)-20) and maxindex>20 and Glistcs[0]<sitfraction*peak_G: #peak found
                            break
                    
                    #for now just check if error message
                    initial_guess=[peak_fit,5e-3,max(Glistcs)]
                    print(f"initial guess {max(Glistcs)*1e6:3g} uS")
                    
                    popt, pcov = scp.optimize.curve_fit(breit_wigner_fkt, csweeplist, Glistcs, p0=initial_guess)
                    

                    peak_fit, hgamma_fit, peak_G_fit=popt
                    fit_detuning=-breit_wigner_detuning(peak_G_fit*sitfraction,peak_G_fit,hgamma_fit)
                    fit_calculated_sitpos=peak_fit+fit_detuning

                    ###########131025 - unitl here -##########
                    peakpos_fit_list.append(peak_fit)
                    fit_detuning_list.append(fit_detuning)

                    if First_run:
                        initial_peakpos=copy.copy(peak_fit)
                    corrected_peakpos_fit=peak_fit-initial_peakpos-correctionV
                    last_measured_corrected_peakpos=copy.copy(corrected_peakpos_fit)
                    last_measured_peakpos=copy.copy(peak_fit)

                    n=0
                    while Glistcs[n]<peak_G*sitfraction: #sit at approx half max of cs peak
                        n=n+1
                    
                        #n=n-1
                    sit_index=copy.copy(n)
                    Glist=Glist+[Glistcs[sit_index]]

                    ####################################
                    gate1_list=[gate1_value]*len(Glistcs)
                    gate2_list=[gate2_value]*len(Glistcs)
                    datasaver_aux1.add_result(('G', Glistcs),('Vg1', gate1_list),('Vg2', gate2_list),(GVg_num_param,nr_GVgs),
                                (CS_V_param,csweeplist))
                    
                    ####################################
                    #print(f"Glistcs[sit_index]={Glistcs[sit_index]}")
                    last_GVg_G=copy.copy(Glistcs[sit_index])
                    Vlist.append(0)
                    Phaselist.append(0)
                    Signal_V_list.append(corrected_peakpos_fit)
                    Signal_V_cbrt_list.append(np.cbrt(corrected_peakpos_fit))
                    V_corr_imp_list.append(0)
                    V_corr_pos_list.append(0)
                    peakG_list.append(peak_G_fit)
                    
                    if start_vgcs-100e-3 < fit_calculated_sitpos < start_vgcs+100e-3:
                        csgate.dc_constant_V(fit_calculated_sitpos)
                    time.sleep(abs(fit_calculated_sitpos-current_csvg)/step_ramp_speed)
                    current_csvg=copy.copy(fit_calculated_sitpos)
                
           
                    if (peak_G_fit<30e-9):
                    #plt.figure(4)
                        print("no peak found")
                        #plt.plot(csweeplist,Glistcs)
                        #plt.plot(last_sitpos, 10e-9, 'bo')
                        #plt.plot(calculated_sitpos, 10e-9, 'go')
                        #plt.show()
                    #
                    Glistscs.append(Glistcs)
                    sweeplists.append(csweeplist)
                    if First_run==False:
                        GVg_sitposlist.append(last_sitpos)
                        GVg_peakposlist.append(peak_fit)
                    
                        #print(f"{current_csvg-last_sitpos}")
                
                    First_run=False
                    First_inner_run=False
                    time_spent_on_GVgs=time_spent_on_GVgs+(time.time_ns()-start_time)/1e9
                
            #Rlist.reverse
            if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
                Glist.reverse()
                Vlist.reverse()
                #Rlist.reverse()
                Phaselist.reverse()
                peakG_list.reverse()
                Vcorr_list.reverse()
                Signal_V_list.reverse()
                #peakpos_list.reverse()
                #detune_V_list.reverse()
                Signal_V_cbrt_list.reverse()
                V_corr_imp_list.reverse()
                V_corr_pos_list.reverse()
                peakpos_fit_list.reverse()
                fit_detuning_list.reverse()
            
            if First_outer_run:
                divergence_list=list(np.zeros(len(Glist)))
                divergence_x_list=list(np.zeros(len(Glist)))
                divergence_2nd_list=list(np.zeros(len(Glist)))
                divergence_x_n_array=list(np.zeros(len(Glist)))
            else:

                #calculate derivative
                divergence_list=[0]#add fist value to account for missing last value
                divergence_x_list=[0]
                divergence_2nd_list=[0,0]
                deriv_ct=5
                divergence_x_n_array=list(np.zeros(deriv_ct+1))
                
                for i in range(len(Signal_V_list)):
                    if i>0:
                        divergence_x_list.append((Signal_V_list[i]-Signal_V_list[i-1]))
                        divergence_list.append((Signal_V_list[i]-Signal_V_list[i-1])+(Signal_V_list[i]-last_Signal_V_list[i]))
                        if i>1 and 'secondlast_Signal_V_list' in locals():
                            divergence_2nd_list.append((Signal_V_list[i]-Signal_V_list[i-2])+(Signal_V_list[i]-secondlast_Signal_V_list[i]))
                        if i>deriv_ct:
                            divergence_x_n_array.append(Signal_V_list[i]-Signal_V_list[i-2])
            if 'last_Signal_V_list' in locals():
                secondlast_Signal_V_list=copy.copy(last_Signal_V_list)
            if len(divergence_2nd_list)==2:
                divergence_2nd_list=list(np.zeros(len(Glist)))
            last_Signal_V_list=copy.copy(Signal_V_list)
            
            First_outer_run=False

            datasaver.add_result(('G', Glist),
                                ('V_r', Vlist),
                                ('Phase', Phaselist),
                                ('V_corr',Vcorr_list),
                                ('signal_shift_V',Signal_V_list),
                                ('signal_shift_V_deriv',divergence_list),
                                ('signal_shift_Vx_deriv',divergence_x_list),
                                ('signal_shift_Vxn_deriv',divergence_x_n_array),
                                ('signal_shift_2pt_deriv',divergence_2nd_list),
                                ('signal_shift_cbrt_V',Signal_V_cbrt_list),
                                ('V_corr_imp',V_corr_imp_list),
                                ('V_corr_pos',V_corr_pos_list),
                                ('peak_Value',peakG_list),
                            #  ('peak_position',peakpos_list),
                                ('peak_fit_position',peakpos_fit_list),
                                #('V_detune',detune_V_list),
                                ('fit_detuning',fit_detuning_list),
                                (gate1_sweep.parameter,gate1_value),
                                (gate2_sweep.parameter,fast_axis_unreversible_list))
            
           
            gate2_sweep.reverse()
            reversed_sweep= not reversed_sweep 
        total_measurement_time=(time.time_ns()-overall_start_time)/1e9
        final_vars_to_save=[time_spent_on_GVgs,time_spent_on_simple_measurement,nr_of_topslips,nr_of_bottomslips,nr_GVgs,nr_simple_meas,total_measurement_time]   
        varnames=[]
        for i in range(len(final_vars_to_save)):
            varnames.append(get_var_name(final_vars_to_save[i]))
        #save_metadata_var(datasaver.dataset,varnames,vars_to_save)


run_number = experiment.last_counter
file_path_y = f".\Data\Raw_data\CD11_D7_C1_\meas{run_number}_conductance.npy"
file_path_x = f".\Data\Raw_data\CD11_D7_C1_\meas{run_number}_vgcs.npy"

# Save the array to the file
#np.save(file_path_y, AllData3d_conductance)
#np.save(file_path_x, AllData3d_Vg)
# Ramp down everything
gate1.dc_slew_rate_V_per_s(ramp_speed)
gate2.dc_slew_rate_V_per_s(ramp_speed)
print(f"time_spent_on_GVgs:{time_spent_on_GVgs}")
print(f"time_spent_on_simple_measurement:{time_spent_on_simple_measurement}")
print(f"topslips:{nr_of_topslips}")
print(f"bottomslips:{nr_of_bottomslips}")
print(f"GVgs:{nr_GVgs}")
print(f"simple measurements:{nr_simple_meas}")
print(f"total measuremnt time(w/o initial ramp):{total_measurement_time}")
#print(f"final_peak_pos:{peakpos_list[-1]}")
#plt.figure(20)
#for i in range(len(sweeplists)):
#    plt.plot(sweeplists[i],Glistcs[i])
#    plt.plot(GVg_sitposlist[i], 10e-9, 'bo')
#    plt.plot(GVg_peakposlist[i], 10e-9, 'ro')
#plt.show()
#now try save metadata - test
#datasets=meas.experiment.data_sets()
#datasets[-1].metadata['new_key'] = 'new_test_value'

#datasets[-1].metadata['new_key'] = 'new_value'

