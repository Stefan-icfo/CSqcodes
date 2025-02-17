# cs_experiment.py

import numpy as np
import time
from tqdm import tqdm
import os

import qcodes as qc
# Import your parameters
import experiment_parameters as params

# Example placeholders for instruments
from instruments import station, zurich, qdac, Triton
from qcodes.dataset import Measurement, new_experiment

# Utility functions
from utils.sample_name import sample_name
from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
from utils.CS_utils import *
from experiment_functions.CS_functions import *
import database

class thermomech_measurement:
    def __init__(self):
        self.area_values_scaled = []
        self.area_values_unscaled = []
        self.area_values_scaled_by_area = []
        self.area_values_scaled_by_slope = []
        self.slopes = []



class CSExperiment:
    """
    A class that copies defaults from experiment_parameters.py
    and provides measurement methods (GVG_fun, other_fun, etc.).
    """

    def __init__(self):
        # Copy from experiment_parameters
        self.device_name = params.device_name
        self.tc = params.tc
        self.tg=params.tg
        self.attn_dB_source = params.attn_dB_source
        self.source_amplitude_instrumentlevel_GVg = params.source_amplitude_instrumentlevel_GVg
        self.mix_down_f = params.mix_down_f
        self.x_avg = params.x_avg
        self.y_avg = params.y_avg
        self.start_vg_cs = params.start_vg_cs
        self.stop_vg_cs = params.stop_vg_cs
        self.step_num_cs = params.step_num_cs
        self.slew_rate=params.slew_rate

        self.sitfraction = params.sitfraction
        self.GVg_data_avg_num=params.data_avg_num
        self.fit_type=params.fit_type
        self.device_name=params.device_name

        self.min_acceptable_peak=params.min_acceptable_peak
        self.cs_gate=qdac.ch06
        self.freq_RLC=1.25e6

        self.idt_point1_x=-1.7467
        self.idt_point1_y=-2.6343
        self.idt_point2_x=-1.74596
        self.idt_point2_y=-2.63382

        self.source_amplitude_CNT = d2v(v2d(np.sqrt(1/2) * self.source_amplitude_instrumentlevel_GVg) - self.attn_dB_source) / 10
       # self.area_values_scaled=[]
       # self.area_values_unscaled=[]
       # self.area_values_scaled_by_area=[]
       # self.area_values_scaled_by_slope=[]
       # self.slopes=[]

    def do_testruns(self):
        qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'sometestruns.db'

    def do_real_measurements(self):
        qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1_part3.db'

    def GVG_fun(
        self,
        start_vg=None,
        stop_vg=None,
        step_num=None,
        device_name=None,
        #run=False,
        save_in_database=True,
        return_data=False,
        return_only_Vg_and_G=True,
        reverse=False,
        pre_ramping_required=False,
        costum_prefix='_'
    ):
        """
        Example measurement that uses self.xxx (from experiment_parameters).
        Ad-hoc overrides can be done by directly changing self.xxx in the run file.
        """
        #if not run:
        #    print("GVG_fun: run=False, skipping measurement.")
        #    return
        gate=self.cs_gate
        tc = self.tc
        vsd_dB = self.attn_dB_source
        amp_lvl = self.source_amplitude_instrumentlevel_GVg
        f_mix = self.mix_down_f
        x_avg = self.x_avg
        y_avg = self.y_avg
        if start_vg==None:
            start_vg = self.start_vg_cs
        if stop_vg==None:
            stop_vg = self.stop_vg_cs
        if step_num==None:
            step_num = self.step_num_cs
        if device_name==None:
            device_name = self.device_name

        # Instrument references
        gate = qdac.ch06
        freq = zurich.oscs.oscs0.freq
        source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
        measured_parameter = zurich.demods.demods0.sample

        # Initialize hardware
        freq(f_mix)
        source_amplitude_param(amp_lvl)
        if pre_ramping_required:
            print(f"Pre-ramping gate to {start_vg}")
            qdac.ramp_multi_ch_slowly(channels=[gate], final_vgs=[start_vg])
        gate.ramp_ch(start_vg)

        vsdac = d2v(v2d(np.sqrt(1/2) * amp_lvl) - vsd_dB) / 10
        vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num=step_num)
        
        prefix_name = 'Conductance_rf_'+costum_prefix
        
        postfix = (f"vsac@inst={amp_lvl*1e3:4g} mV",
            f"_g1={qdac.ch01.dc_constant_V():4g},"
            f"g2={qdac.ch01.dc_constant_V():4g},"
            f"g3={qdac.ch01.dc_constant_V():4g},"
            f"g4={qdac.ch01.dc_constant_V():4g},"
            f"g5={qdac.ch01.dc_constant_V():4g}"
        )
        postfix_str = "".join(postfix)
        gate.label = 'cs_gate'
        exp_name = exp_name=prefix_name+device_name+postfix_str

        # Prepare data arrays
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
                #varnames = [str(name) for name in [tc, vsd_dB, amp_lvl, x_avg, y_avg]]
                #save_metadata_var(datasaver.dataset, varnames, [tc, vsd_dB, amp_lvl, x_avg, y_avg])

                for vgdc_value in tqdm(vgdc_sweep, desc='Gate voltage Sweep'):
                    gate.ramp_ch(vgdc_value)
                    time.sleep(1.1 * tc)

                    # Some measurement
                    #_ = measured_parameter()
                    theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                    R = 1 / G

                    datasaver.add_result(
                        ('R', R),
                        ('G', G),
                        ('V_r', v_r_calc),
                        ('Phase', theta_calc),
                        ('I', I),
                        (vgdc_sweep.parameter, vgdc_value)
                    )

                    if return_data:
                        Glist.append(G)
                        Vlist.append(v_r_calc)
                        Ilist.append(I)
                        Phaselist.append(theta_calc)
                        Rlist.append(R)
        else:
            # Not saving in DB
            for vgdc_value in tqdm(vgdc_sweep, desc='Gate voltage Sweep'):
                gate.ramp_ch(vgdc_value)
                time.sleep(1.1 * tc)

                #_ = measured_parameter()
                theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                R = 1 / G

                if return_data:
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

        if return_data:
            Vglist = list(vgdc_sweep)
            if return_only_Vg_and_G:
                return np.array(Vglist), np.array(Glist)
            else:
                return (
                    np.array(Vglist),
                    np.array(Glist),
                    np.array(Vlist),
                    np.array(Ilist),
                    np.array(Phaselist),
                    np.array(Rlist),
                )

   

    def do_GVg_and_adjust_sitpos(
            self,
            fit_type=None,
            initial_guess=None, 
            sitfraction=None,
            start_vg=None,
            stop_vg=None,
            step_num=None,
            exp_name=None,
            pre_ramping_required=False,
            device_name=None,
            save_in_database=True,
            return_full_data=False,
            return_data=False,
            reverse=False,
            data_avg_num=None,
            sit_side="left",
            costum_prefix='_',
            testplot=False
            ):
        tc = self.tc
        vsd_dB = self.attn_dB_source
        amp_lvl = self.source_amplitude_instrumentlevel_GVg
        f_mix = self.mix_down_f
        gate=self.cs_gate
        #x_avg = self.x_avg
        #y_avg = self.y_avg
        min_acceptable_peak=self.min_acceptable_peak
        if device_name==None:
            device_name = self.device_name
        if start_vg==None:
            start_vg = self.start_vg_cs
        if stop_vg==None:
            stop_vg = self.stop_vg_cs
        if step_num==None:
            step_num = self.step_num_cs
        if fit_type==None:
            fit_type=self.fit_type
        if data_avg_num==None:
            data_avg_num=self.GVg_data_avg_num
        if sitfraction==None:
            sitfraction=self.sitfraction
        
        prefix_name = 'Conductance_rf_'+costum_prefix
        postfix = (f"vsac@inst={amp_lvl*1e3:4g} mV",
            f"_g1={qdac.ch01.dc_constant_V():4g},"
            f"g2={qdac.ch01.dc_constant_V():4g},"
            f"g3={qdac.ch01.dc_constant_V():4g},"
            f"g4={qdac.ch01.dc_constant_V():4g},"
            f"g5={qdac.ch01.dc_constant_V():4g}"
        )
        postfix_str = "".join(postfix)
        exp_name = exp_name=prefix_name+device_name+postfix_str
    
        
        
        Vg,G_vals=self.GVG_fun(start_vg=start_vg,
                stop_vg=stop_vg,
                step_num=step_num,
                pre_ramping_required=pre_ramping_required,
                save_in_database=False,
                return_data=True,
                return_only_Vg_and_G=True,
                reverse=False,
                )
        
        #print("GVg done")
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

            """

            avg_G=centered_moving_average(G_vals,data_avg_num)
            fit_vals=avg_G
            max_avg=max(avg_G)
            deriv_avg=avg_G[:-1] - avg_G[1:]

        if isinstance(sitfraction, (int, float)):
            print("sitfraction is a number")
            left_idx = np.argmax(avg_G > max_avg*sitfraction)
            sitpos=Vg[left_idx]
            x=[Vg[left_idx-data_avg_num:left_idx+data_avg_num]]
            y=[G_vals[left_idx-data_avg_num:left_idx+data_avg_num]]
            if left_idx == 0:
                raise ValueError("left_idx=0. probably no peak found")
            #print(f"left_idx: {left_idx}, max_avg: {max_avg}, sitfraction: {sitfraction}")
            #print(f"x slice indices: {left_idx-data_avg_num} to {left_idx+data_avg_num}")
            #print(f"x values: {Vg[left_idx-data_avg_num:left_idx+data_avg_num]}")
            #print(f"y values: {G_vals[left_idx-data_avg_num:left_idx+data_avg_num]}")
            #y_div = [y_item * 1e7 for y_item in y]
            result=scp.stats.linregress(x,y)#try y*1e7, result/1e7
            slope=result.slope#

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


        """
        fit_vals,slope,sitpos,pos_idx=find_sitpos_from_avg_data(Vg,G_vals,sitfraction=sitfraction,data_avg_num=data_avg_num,sit_side=sit_side,return_avg_data=True)
        gate.ramp_ch(sitpos) 

        if save_in_database:
            vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num = step_num)
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
                #var_names=names_of_vars_to_save.split(',')
                #for varname,var in zip(var_names,vars_to_save):
                #    datasaver.dataset.add_metadata(varname,var)
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
                
                slope_array[approx_sitpos_index-round(data_avg_num/2)] = fit_vals[approx_sitpos_index] - slope * (Vg[approx_sitpos_index] - Vg[approx_sitpos_index-round(data_avg_num/2)])
                slope_array[approx_sitpos_index+round(data_avg_num/2)] = fit_vals[approx_sitpos_index] + slope * (Vg[approx_sitpos_index+round(data_avg_num/2)] - Vg[approx_sitpos_index])
                #print(f"test:{slope_array[approx_sitpos_index]}")
                #plt.plot(Vg,slope_array)
                #plt.show()
                datasaver.add_result(('G', G_vals), ('fit',fit_vals),('sitpos',approx_sitpos_array),('slope',slope_array), (vgdc_sweep.parameter, Vg))
                datasaver.dataset.add_metadata("slope",slope)
                datasaver.dataset.add_metadata("sitpos",sitpos)

                if testplot:
                    run_id = datasaver.run_id
                    foldername=f'C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_\\meas{run_id}_'
                    if not os.path.exists(foldername):
                        os.makedirs(foldername) 
                    filename=f'meas{run_id}_slopeplot.png'
                    path = os.path.join(foldername, filename)
                
                    plt.plot(Vg,G_vals)
                    plt.plot([Vg[approx_sitpos_index-round(data_avg_num/2)],Vg[approx_sitpos_index+round(data_avg_num/2)]],[slope_array[approx_sitpos_index-round(data_avg_num/2)],slope_array[approx_sitpos_index+round(data_avg_num/2)]])
                    #plt.show()
                    plt.savefig(path)
                    plt.close()

        if return_full_data:
            return Vg,G_vals,popt, pcov,slope,sitpos
        else:
            return slope,sitpos


    def GVG_fun_sensitivity(
        self,
        start_vg=None,
        stop_vg=None,
        step_num=None,
        device_name=None,
        #run=False,
        save_in_database=True,
        return_data=False,
        #return_only_Vg_and_G=True,
        return_only_Vg_G_and_Isens=True,
        reverse=False,
        pre_ramping_required=False,
        costum_prefix='_',
        sens_demod=zurich.demod2,
        RF_sens_osc=zurich.freq2,
        mod_gate=None,
        mod_amplitude=100e-6,
        mod_frequency=5e3,
        RF_meas_osc=zurich.freq0,
        RF_drive_osc=zurich.freq1,
        drive_type="LF"#LF for qdac sine wave, RF for zurich
        ):
        """
        Example measurement that uses self.xxx (from experiment_parameters).
        Ad-hoc overrides can be done by directly changing self.xxx in the run file.
        """
        #if not run:
        #    print("GVG_fun: run=False, skipping measurement.")
        #    return
        gate=self.cs_gate
        tc = self.tc
        vsd_dB = self.attn_dB_source
        amp_lvl = self.source_amplitude_instrumentlevel_GVg
        f_mix = self.mix_down_f
        x_avg = self.x_avg
        y_avg = self.y_avg
        if start_vg==None:
            start_vg = self.start_vg_cs
        if stop_vg==None:
            stop_vg = self.stop_vg_cs
        if step_num==None:
            step_num = self.step_num_cs
        if device_name==None:
            device_name = self.device_name

        #######
        if mod_gate==None:
            mod_gate=self.cs_gate

        # Instrument references
        gate = self.cs_gate
        freq = zurich.oscs.oscs0.freq
        source_amplitude_param = zurich.sigouts.sigouts0.amplitudes.amplitudes0.value
        measured_parameter = zurich.demods.demods0.sample

        # Initialize hardware
        freq(f_mix)
        source_amplitude_param(amp_lvl)
        if pre_ramping_required:
            print(f"Pre-ramping gate to {start_vg}")
            qdac.ramp_multi_ch_slowly(channels=[gate], final_vgs=[start_vg])
        gate.ramp_ch(start_vg)

        


        if drive_type=="LF":
            RF_sens_osc(self.freq_RLC+mod_frequency)#set sensing channel frequency
            ################init sinewave generator#########################
            # Define amplitude and frequency for the sine wave
            # Configure the sine wave 
            sine_wave_context = mod_gate.sine_wave(
            frequency_Hz=mod_frequency,  # Set frequency to 5 kHz
            span_V=mod_amplitude,        # Set amplitude span to 10 mV
            offset_V=0.0,            # No offset, centered around 0V
            repetitions=-1,          # Run indefinitely (-1 for infinite repetitions)
            )
            # Start the sine wave
            sine_wave_context.start()

        if drive_type=="RF":
            #set up mixdown
            RF_sens_osc(self.freq_RLC)
            RF_meas_osc(mod_frequency-self.freq_RLC)
            RF_drive_osc(mod_frequency)
            zurich.output1_amp1(mod_amplitude)
            zurich.sigout1_amp1_enabled_param.value(1)




        vsdac = d2v(v2d(np.sqrt(1/2) * amp_lvl) - vsd_dB) / 10
        vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num=step_num)
        if reverse:
            vgdc_sweep.reverse()
        

        prefix_name = 'Conductance_LFsens_'+costum_prefix
        
        postfix = (f"vsac@inst={amp_lvl*1e3:4g} mV",
            f"_g1={qdac.ch01.dc_constant_V():4g},"
            f"g2={qdac.ch01.dc_constant_V():4g},"
            f"g3={qdac.ch01.dc_constant_V():4g},"
            f"g4={qdac.ch01.dc_constant_V():4g},"
            f"g5={qdac.ch01.dc_constant_V():4g},"
            f"freq={mod_frequency:4g},"
            f"amp={mod_amplitude:4g},"
        )
        postfix_str = "".join(postfix)
        gate.label = 'cs_gate'
        exp_name = exp_name=prefix_name+device_name+postfix_str

        # Prepare data arrays
        Glist, Vlist, Ilist, Phaselist, Rlist, V_sens_list, I_sens_list, Phase_sens_list = [], [], [], [], [], [], [], []

        

        if save_in_database:
            experiment = new_experiment(name=exp_name, sample_name=device_name)
            meas = Measurement(exp=experiment)
            meas.register_parameter(vgdc_sweep.parameter)
            meas.register_custom_parameter('G', unit='S', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('V_r', unit='V', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('Phase', unit='rad', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('R', unit='Ohm', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('I', unit='A', setpoints=[vgdc_sweep.parameter])

            meas.register_custom_parameter('V_r_sens', unit='V', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('Phase_sens', unit='rad', setpoints=[vgdc_sweep.parameter])
            meas.register_custom_parameter('I_sens', unit='A', setpoints=[vgdc_sweep.parameter])

            with meas.run() as datasaver:
                #varnames = [str(name) for name in [tc, vsd_dB, amp_lvl, x_avg, y_avg]]
                #save_metadata_var(datasaver.dataset, varnames, [tc, vsd_dB, amp_lvl, x_avg, y_avg])

                

                for vgdc_value in tqdm(vgdc_sweep, desc='Gate voltage Sweep'):
                    gate.ramp_ch(vgdc_value)
                    time.sleep(1.1 * tc)

                    # Some measurement
                    #_ = measured_parameter()
                    theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                    R = 1 / G

                    sens_value=sens_demod()
                    theta_sens, v_r_sens, I_sens, _ = zurich.phase_voltage_current_conductance_compensate(vsdac,x_avg=0,y_avg=0,measured_value=sens_value)

                    datasaver.add_result(
                        ('R', R),
                        ('G', G),
                        ('V_r', v_r_calc),
                        ('V_r_sens', v_r_sens),
                        ('Phase', theta_calc),
                        ('Phase_sens', theta_sens),
                        ('I', I),
                        ('I_sens', I_sens),
                        (vgdc_sweep.parameter, vgdc_value)
                    )

                    if return_data:
                        Glist.append(G)
                        Vlist.append(v_r_calc)
                        Ilist.append(I)
                        Phaselist.append(theta_calc)
                        Rlist.append(R)
        else:
            # Not saving in DB
            for vgdc_value in tqdm(vgdc_sweep, desc='Gate voltage Sweep'):
                gate.ramp_ch(vgdc_value)
                time.sleep(1.1 * tc)

                #_ = measured_parameter()
                theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                R = 1 / G
                sens_value=sens_demod()
                theta_sens, v_r_sens, I_sens, _ = zurich.phase_voltage_current_conductance_compensate(vsdac,x_avg=0,y_avg=0,measured_value=sens_value)

                if return_data:
                    Glist.append(G)
                    Vlist.append(v_r_calc)
                    Ilist.append(I)
                    Phaselist.append(theta_calc)
                    Rlist.append(R)

                    V_sens_list.append(v_r_sens)
                    I_sens_list.append(I_sens)
                    Phase_sens_list.append(theta_sens)

        if return_data and reverse:
            Glist.reverse()
            Vlist.reverse()
            Ilist.reverse()
            Phaselist.reverse()
            Rlist.reverse()

            V_sens_list.reverse()
            I_sens_list.reverse()
            Phase_sens_list.reverse()

        if return_data:
            Vglist = list(vgdc_sweep)
            if return_only_Vg_G_and_Isens:
                return np.array(Vglist), np.array(Glist), np.array(I_sens_list)
            else:
                return (
                    np.array(Vglist),
                    np.array(Glist),
                    np.array(Vlist),
                    np.array(Ilist),
                    np.array(Phaselist),
                    np.array(Rlist),
                    np.array(I_sens_list),
                    np.array(V_sens_list),
                    np.array(Phase_sens_list)
                )
            
        if drive_type=="LF":
            sine_wave_context.abort()

        if drive_type=="RF":
            zurich.sigout1_amp1_enabled_param.value(0)