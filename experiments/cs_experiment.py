# cs_experiment.py

import numpy as np
import time
from tqdm import tqdm
import os

import qcodes as qc
# Import your parameters
import experiment_parameters as params

# Example placeholders for instruments
from instruments import station, qdac, Triton
from qcodes.dataset import Measurement, new_experiment

# Utility functions
from utils.sample_name import sample_name
#from utils.zi_uhfli_GVg_setup import zi_uhfli_GVg_setup
from utils.d2v import d2v
from utils.v2d import v2d
from utils.rms2pk import rms2pk
from utils.CS_utils import *
from experiment_functions.CS_functions import *
import database

import sys,signal
import copy

class CSExperiment:
    def __init__(self):
        # Set any constants that don't come from params
        self.cs_gate = qdac.ch06
        self.max_thermomech_freq = 160e6
        
        # Load all parameters from params module
        self.load_parameters()

    
        
    def load_parameters(self):
        import importlib
        importlib.reload(params)  
    
        # Update all parameters from the params module
        self.device_name = params.device_name
        self.tc = params.tc
        self.tg = params.tg
        self.attn_dB_source = params.attn_dB_source
        self.source_amplitude_instrumentlevel_GVg = params.source_amplitude_instrumentlevel_GVg
        self.mix_down_f = params.mix_down_f
        self.x_avg = params.x_avg
        self.y_avg = params.y_avg
        self.start_vg_cs = params.start_vg_cs
        self.stop_vg_cs = params.stop_vg_cs
        self.step_num_cs = params.step_num_cs
        self.slew_rate = params.slew_rate
        self.sitfraction = params.sitfraction
        self.GVg_data_avg_num = params.data_avg_num
        self.fit_type = params.fit_type
        self.device_name = params.device_name  # This line is duplicated in the original
        self.min_acceptable_peak = params.min_acceptable_peak
        self.freq_RLC = params.RLC_frequency
        self.idt_point1_x = params.idt_point1_x
        self.idt_point1_y = params.idt_point1_y
        self.idt_point2_x = params.idt_point2_x
        self.idt_point2_y = params.idt_point2_y
        self.start_f = params.start_f
        self.stop_f = params.stop_f
        self.step_num_f = params.step_num_f
        self.freq_sweep_avg_num = params.freq_sweep_avg_num
        self.max_ramp_speed=params.max_ramp_speed
        self.ramp_step_size=params.ramp_step_size
        self.costum_prefix=params.costum_prefix

        self.pre_ramping_required=params.pre_ramping_required
        
        # linesweep parameters
        self.start_vgo_ls = params.start_vgo_ls
        self.stop_vgo_ls = params.stop_vgo_ls
        self.step_vgo_num_ls = params.step_vgo_num_ls
        self.start_vgi_ls = params.start_vgi_ls
        self.stop_vgi_ls = params.stop_vgi_ls
        self.step_vgi_num_ls = params.step_vgi_num_ls
        self.start_vgi_scan_ls = params.start_vgi_scan_ls
        self.scan_range_ls = params.scan_range_ls
        self.increments_ls = params.increments_ls
    
        # Recalculate derived parameters
        self.source_amplitude_CNT = d2v(v2d(np.sqrt(1/2) * self.source_amplitude_instrumentlevel_GVg) - self.attn_dB_source)
        try:
            self.tempMC=Triton.MC()
        except:
            print("cant measure MC")
            self.tempMC="temp measurement error"
        
        return self



    def save_all_parameters_to_metadata(self, datasaver):
        """
        Saves all instance attributes to dataset metadata,
        skipping non-serializable types (only saves strings, ints, and floats).
        """
        for key, value in self.__dict__.items():
            if isinstance(value, (int, float, str)):
                datasaver.dataset.add_metadata(key, value)
            elif hasattr(value, 'item'):  # Handle numpy scalar types
                try:
                    datasaver.dataset.add_metadata(key, value.item())
                except Exception as e:
                    print(f"Warning: Could not save {key} (numpy scalar) due to: {e}")
            #else:
            #    print(f"Skipping metadata key '{key}' (unsupported type: {type(value)})")

    def print_parameters(self):
        """
        Prints all parameters and their values defined in the __init__ method.
        Excludes commented out parameters and other methods.
        """
        print(f"Parameters for {self.__class__.__name__}:")
        print("-" * 40)
        
        # Get all instance attributes defined in __init__
        for attr_name, attr_value in self.__dict__.items():
            # Format the output based on the type of value
            if isinstance(attr_value, (int, float)) and abs(attr_value) > 1000:
                # Scientific notation for large numbers
                print(f"{attr_name}: {attr_value:.4e}")
            elif isinstance(attr_value, float):
                print(f"{attr_name}: {attr_value:.6f}")
            else:
                print(f"{attr_name}: {attr_value}")
        
        print("-" * 40)


    def do_testruns(self):
        qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'sometestruns.db'

    def do_real_measurements(self):
        qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD13_E3_C2.db'


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
        costum_prefix='_',
        load_params=True
    ):
        """
        Example measurement that uses self.xxx (from experiment_parameters).
        Ad-hoc overrides can be done by directly changing self.xxx in the run file.
        """
        #if not run:
        #    print("GVG_fun: run=False, skipping measurement.")
        #    return
        if load_params:
            self.load_parameters()
        gate=self.cs_gate
        tc = self.tc
        vsd_dB = self.attn_dB_source
        amp_lvl = self.source_amplitude_instrumentlevel_GVg
        f_mix = self.mix_down_f
        if start_vg==None:
            start_vg = self.start_vg_cs
        if stop_vg==None:
            stop_vg = self.stop_vg_cs
        if step_num==None:
            step_num = self.step_num_cs
        if device_name==None:
            device_name = self.device_name
        if pre_ramping_required==None:
            pre_ramping_required = self.pre_ramping_required
        


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
            qdac.ramp_multi_ch_slowly(channels=[gate], final_vgs=[start_vg],step_size=self.ramp_step_size,ramp_speed=self.max_ramp_speed)
        gate.ramp_ch(start_vg)
        

        vsdac = d2v(v2d(np.sqrt(1/2) * amp_lvl) - vsd_dB) 
        vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num=step_num)
        
        prefix_name = 'Conductance_rf_'+costum_prefix
        
        postfix = (f"vsac@inst={amp_lvl*1e3:.4g} mV",
            f"_g1={qdac.ch01.dc_constant_V():.4g},"
            f"g2={qdac.ch02.dc_constant_V():.4g},"
            f"g3={qdac.ch03.dc_constant_V():.4g},"
            f"g4={qdac.ch04.dc_constant_V():.4g},"
            f"g5={qdac.ch05.dc_constant_V():.4g}"
            f"g7={qdac.ch07.dc_constant_V():.4g}"
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
                qdac.add_dc_voltages_to_metadata(datasaver=datasaver)
                zurich.save_config_to_metadata(datasaver=datasaver)
                self.save_all_parameters_to_metadata(datasaver=datasaver)
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
            data_avg_num=None,
            sit_side="left",
            costum_prefix=None,
            testplot=False,
            load_params=True
            ):
        if load_params:
            self.load_parameters()
        
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
            print("fit type: "+fit_type)
        if data_avg_num==None:
            data_avg_num=self.GVg_data_avg_num
        if sitfraction==None:
            sitfraction=self.sitfraction
        if costum_prefix==None:
            costum_prefix=self.costum_prefix
        
        prefix_name = 'Conductance_rf_'+costum_prefix
        postfix = (f"vsac@inst={amp_lvl*1e3:.4g} mV",
            f"_g1={qdac.ch01.dc_constant_V():.4g},"
            f"g2={qdac.ch02.dc_constant_V():.4g},"
            f"g3={qdac.ch03.dc_constant_V():.4g},"
            f"g4={qdac.ch04.dc_constant_V():.4g},"
            f"g5={qdac.ch05.dc_constant_V():.4g}"
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
            print("fitted to tunnel peak")
        if fit_type=='thermal':
            popt, pcov,slope,sitpos=fit_and_find_sitpos_singlepeak_thermal(Vg,G_vals,initial_guess=initial_guess, sitfraction=sitfraction,return_full_fit_data=True)
            fit_vals=thermal_CB_peak(Vg,popt[0],popt[1],popt[2])
            print("fitted thermal peak")
        if fit_type=='data':
            popt,pcov=None,None
            print("fitted to averaged data")

 
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
                self.save_all_parameters_to_metadata(datasaver=datasaver)
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
        pre_ramping_required=True,
        costum_prefix='_',
        sens_demod=zurich.demod2,
        RF_sens_osc=zurich.freq2,
        mod_gate=None,
        mod_amplitude=0.1e-3,
        mod_frequency=1e3,
        RF_meas_osc=zurich.freq0,
        RF_drive_osc=zurich.freq1,
        drive_type="LF",#LF for qdac sine wave, RF for zurich
        load_params=True):
        """
        Example measurement that uses self.xxx (from experiment_parameters).
        Ad-hoc overrides can be done by directly changing self.xxx in the run file.
        """
        if load_params:
            self.load_parameters()
        
        #if not run:
        #    print("GVG_fun: run=False, skipping measurement.")
        #    return
        gate=self.cs_gate
        tc = self.tc
        vsd_dB = self.attn_dB_source
        amp_lvl = self.source_amplitude_instrumentlevel_GVg
        f_mix = self.mix_down_f
       # x_avg = self.x_avg
       # y_avg = self.y_avg
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
            qdac.ramp_multi_ch_slowly(channels=[gate], final_vgs=[start_vg],step_size=self.ramp_step_size,ramp_speed=self.max_ramp_speed)
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




        vsdac = d2v(v2d(np.sqrt(1/2) * amp_lvl) - vsd_dB) 
        vgdc_sweep = gate.dc_constant_V.sweep(start=start_vg, stop=stop_vg, num=step_num)
        if reverse:
            vgdc_sweep.reverse()
        

        prefix_name = 'Conductance_sens_'+costum_prefix
        
        postfix = (f"vsac@inst={amp_lvl*1e3:.4g} mV",
            f"_g1={qdac.ch01.dc_constant_V():.4g},"
            f"g2={qdac.ch02.dc_constant_V():.4g},"
            f"g3={qdac.ch03.dc_constant_V():.4g},"
            f"g4={qdac.ch04.dc_constant_V():.4g},"
            f"g5={qdac.ch05.dc_constant_V():.4g},"
            f"g7={qdac.ch07.dc_constant_V():.4g},"
            f"freq={mod_frequency:.4g},"
            f"amp={mod_amplitude:.4g},"
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
                qdac.add_dc_voltages_to_metadata(datasaver=datasaver)
                zurich.save_config_to_metadata(datasaver=datasaver)
                self.save_all_parameters_to_metadata(datasaver=datasaver)

                

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
                        V_sens_list.append(v_r_sens)
                        I_sens_list.append(I_sens)
                        Phase_sens_list.append(theta_sens)
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
        if drive_type=="LF":
            sine_wave_context.abort()
            RF_sens_osc(self.freq_RLC)#set back sensing channel frequency

        if drive_type=="RF":
            zurich.sigout1_amp1_enabled_param.value(0)

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
            

    def sit_at_max_Isens(self,avg_num=3,return_sitpos=True,side=None,start_vg=None,stop_vg=None,step_num=None):
        if start_vg==None:
            start_vg = self.start_vg_cs
        if stop_vg==None:
            stop_vg = self.stop_vg_cs
        if step_num==None:
            step_num = self.step_num_cs
        Vg,G,Isens=self.GVG_fun_sensitivity(start_vg=start_vg,stop_vg=stop_vg,step_num=step_num,
        save_in_database=True,
        return_data=True,
        #return_only_Vg_and_G=True,
        return_only_Vg_G_and_Isens=True,
        costum_prefix='sens_sitpos')
        Gmax_id=np.argmax(centered_moving_average(G,n=avg_num))

        I_sens_avg=centered_moving_average(Isens,n=avg_num)
        if side==None:
            Imax_id=np.argmax(I_sens_avg)
            VmaxI=Vg[Imax_id]
            self.cs_gate.ramp_ch(VmaxI)
            print(f"V_max_v {VmaxI}")
        if side=="left":
            I_sens_avg=I_sens_avg[1:Gmax_id]
            Imax_id=np.argmax(I_sens_avg)
            VmaxI=Vg[Imax_id]
            self.cs_gate.ramp_ch(VmaxI)
            print(f"V_max_v_left {VmaxI}")
        if side=="right":
            I_sens_avg=I_sens_avg[Gmax_id:-1]
            Imax_id=np.argmax(I_sens_avg)
            VmaxI=Vg[Gmax_id + Imax_id]
            self.cs_gate.ramp_ch(VmaxI)
            print(f"V_max_v_right {VmaxI}")
        time.sleep(5)
        if return_sitpos:
            return VmaxI


    def mech_simple_fun_db(
            self,
            device_name=None,
            costum_prefix=None,
            source_amplitude_param=zurich.output0_amp0,
            gate_amplitude_param=zurich.output1_amp1,
            gate_rf_enabled_param = getattr(zurich.sigouts.sigouts1.enables, f'enables{1}'),
            freq_mech = zurich.oscs.oscs1.freq,
            freq_rf = zurich.freq0,
            freq_rlc = zurich.freq2,
            drive_amp_at_instr = None,
            start_f=None,
            stop_f=None,
            step_num_f=None,
            measured_parameter=zurich.demod2,
            switch_onoff_g2=True,
            load_params=True,
            Delft=False,
            return_I_and_f=False
        ):
        if load_params:
            self.load_parameters()
        if device_name==None:
            device_name = self.device_name
        if start_f==None:
            start_f = self.start_f
        if stop_f==None:
            stop_f = self.stop_f
        if step_num_f==None:
            step_num_f= self.step_num_f
        if costum_prefix==None:
            costum_prefix=self.costum_prefix
        if not (drive_amp_at_instr==None):
            gate_amplitude_param(drive_amp_at_instr)#sets the drive ampitude if none is given. not yet tested
        #zurich.freq2(freq_rlc) 

        if start_f<1e6 or stop_f<1e6:
            print("FREQUENCY LOW!! MAYBE FOERGOT THE e6??")
        tc=self.tc
        vsdac=self.source_amplitude_CNT
        freq_sweep_avg_nr=self.freq_sweep_avg_num

        postfix = f"_{round(gate_amplitude_param()*1000,3)}mV on gate@inst,_{round(source_amplitude_param()*1000,3)}mV on source@inst, g1={round(qdac.ch01.dc_constant_V(),2)},g2={round(qdac.ch02.dc_constant_V(),5)},g3={round(qdac.ch03.dc_constant_V(),2)},g4={round(qdac.ch04.dc_constant_V(),5)},g5={round(qdac.ch05.dc_constant_V(),2)},gcs={round(qdac.ch06.dc_constant_V(),5)}"
        exp_name = sample_name(costum_prefix+"_simple_mech_"+postfix)
        freq_sweep = freq_rf.sweep(start=start_f, stop=stop_f, num = step_num_f)
        # ----------------Create a measurement-------------------------
        experiment = new_experiment(name=exp_name, sample_name=device_name)
        meas = Measurement(exp=experiment)
        meas.register_parameter(freq_sweep.parameter)  # 
        meas.register_custom_parameter('I_rf', 'current', unit='I', basis=[], setpoints=[freq_sweep.parameter])
        meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[freq_sweep.parameter])
        meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[freq_sweep.parameter])
        meas.register_custom_parameter('I_rf_avg', 'current_avg', unit='I_avg', basis=[], setpoints=[freq_sweep.parameter])

        with meas.run() as datasaver:
            qdac.add_dc_voltages_to_metadata(datasaver=datasaver)
            zurich.save_config_to_metadata(datasaver=datasaver)
            self.save_all_parameters_to_metadata(datasaver=datasaver)
            datasaver.dataset.add_metadata('gate_rf_enabled_param__',gate_rf_enabled_param.value())
            if switch_onoff_g2:
                zurich.sigout1_amp1_enabled_param.value(1)
            print(f"gate 2 on? {gate_rf_enabled_param.value()}")
            if gate_rf_enabled_param.value()==0:
                print("GATE 2 IS OFF!!")
            # for i in range(2):
            I_list=[]
            if Delft:
                freq_rf(freq_rlc())

            for f_value in tqdm(freq_sweep, leave=False, desc='Frequency Sweep', colour = 'green'):
                if not Delft:
                    freq_rf(f_value-freq_rlc())
                freq_mech(f_value)
                #zurich.oscs.oscs4.freq(f_value+freq_rlc())
                time.sleep(1.1*tc) # Wait 1.1 times the time contanst of the lock-in
                measured_value=measured_parameter()
                theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac=vsdac,measured_value=zurich.demods.demods2.sample())
                        
                #G calculation
                I_list.append(I)
            
                datasaver.add_result(('I_rf', I),
                                    ('V_r', v_r_calc),
                                    ('Phase', theta_calc),
                                    (freq_sweep.parameter,f_value))
            I_avg=centered_moving_average(I_list,n=freq_sweep_avg_nr)

            datasaver.add_result(('I_rf_avg', I_avg),(freq_sweep.parameter,list(freq_sweep)))#try this first
            if switch_onoff_g2:
                zurich.sigout1_amp1_enabled_param.value(0)
            if return_I_and_f:
                return np.array(I_list), np.array(list(freq_sweep))
                
        

    def linesweep_parallel_LFsens(self,#in construction
                           device_name=None,
                           costum_prefix='_',
                           start_vgo =  None,#
                           stop_vgo =   None,#
                            step_vgo_num = None,
                            start_vgi = None,#-0.788
                            stop_vgi = None,#-0.776
                            step_vgi_num = None,
                            start_vgi_scan=None,#first guess for peak
                            scan_range=None,
                            mod_amplitude=0.1e-3,
                            mod_frequency=1e3,
                            increments=[0,0,0,0],
                            main_gate=qdac.ch01.dc_constant_V,
                            aux_gates=[],
                            pre_ramping_required=True,
                            return_max=True,
                            load_params=True,
                            unconditional_end_ramp_Vgo=None
             ):
        if load_params:
            self.load_parameters()
        if device_name==None:
            device_name = self.device_name
        if start_vgo == None:
            start_vgo = self.start_vgo_ls
        if stop_vgo == None:
            stop_vgo = self.stop_vgo_ls
        if step_vgo_num == None:
            step_vgo_num = self.step_vgo_num_ls
        if start_vgi == None:
            start_vgi = self.start_vgi_ls
        if stop_vgi == None:
            stop_vgi = self.stop_vgi_ls
        if step_vgi_num == None:
            step_vgi_num = self.step_vgi_num_ls
        if start_vgi_scan == None:
            start_vgi_scan = self.start_vgi_scan_ls
        if scan_range == None:
            scan_range = self.scan_range_ls
        if increments == None:
            increments = self.increments_ls

        vsdac=self.source_amplitude_CNT
        mod_gate=self.cs_gate
        sens_demod=zurich.demod2
        RF_sens_osc=zurich.freq2
        RF_meas_osc=zurich.freq0
        RF_drive_osc=zurich.freq1
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

        tc=self.tc
        #tg=self.tg
        slew_rate=self.slew_rate
        postfix="_linesweep_p_sens"
        inner_gate=self.cs_gate.dc_constant_V

        step_vgo=np.absolute((start_vgo-stop_vgo)/step_vgo_num)
        step_vgi=np.absolute((start_vgi-stop_vgi)/step_vgi_num)
        lower_boundary=start_vgi_scan-scan_range/2
        upper_boundary=start_vgi_scan+scan_range/2
        print(f'Scanning over {step_vgi_num*scan_range/(stop_vgi-start_vgi)} points in vgi')
       
        if pre_ramping_required:
            print(f"Pre-ramping main gate to {start_vgo} and csgate to initial lower boundary {lower_boundary}")
            qdac.ramp_multi_ch_slowly(channels=[main_gate.instrument,inner_gate.instrument], final_vgs=[start_vgo,lower_boundary],step_size=self.ramp_step_size,ramp_speed=self.max_ramp_speed)

        main_gate(start_vgo)
        inner_gate(lower_boundary)
        #for auxgate,increment in zip(aux_gates,increments):
        #    auxgate(start_vgo)
        time.sleep(10)
        #main_gate.label = 'main_gate' # Change the label of the gate chanel
        inner_gate.label = 'CS(inner)' # Change the label of the source chaneel

        exp_name = costum_prefix+postfix
        outer_gate_sweep=main_gate.sweep(start=start_vgo, stop=stop_vgo, num = step_vgo_num)
        inner_gate_sweep=inner_gate.sweep(start=start_vgi, stop=stop_vgi, num = step_vgi_num)
        experiment = new_experiment(name=exp_name, sample_name=device_name)
        meas = Measurement(exp=experiment)
        meas.register_parameter(outer_gate_sweep.parameter)  # 
        meas.register_parameter(inner_gate_sweep.parameter)  # 
        meas.register_custom_parameter('G', 'G', unit='S', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
        meas.register_custom_parameter('V_r', 'Amplitude', unit='V', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
        meas.register_custom_parameter('Phase', 'Phase', unit='rad', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
        #meas.register_custom_parameter('temperature', 'T', unit='K', basis=[], setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
        #meas.register_custom_parameter('V_r_sens', unit='V', setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
        #meas.register_custom_parameter('Phase_sens', unit='rad', setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])
        meas.register_custom_parameter('I_sens', unit='A', setpoints=[outer_gate_sweep.parameter,inner_gate_sweep.parameter])

        # # -----------------Start the Measurement-----------------------
        
        # inverse_source_sweep=source_sweep.reverse() # or define function
        
        with meas.run() as datasaver:

            def cleanup():
                sine_wave_context.abort()
                RF_sens_osc(self.freq_RLC)
                datasaver.dataset.add_metadata(f'endVg1',qdac.ch01.dc_constant_V())
                datasaver.dataset.add_metadata(f'endVg2',qdac.ch02.dc_constant_V())
                datasaver.dataset.add_metadata(f'endVg3',qdac.ch03.dc_constant_V())
                datasaver.dataset.add_metadata(f'endVg4',qdac.ch04.dc_constant_V())
                datasaver.dataset.add_metadata(f'endVg5',qdac.ch05.dc_constant_V())
                datasaver.dataset.add_metadata(f'endVg6',qdac.ch06.dc_constant_V())
                if unconditional_end_ramp_Vgo is not None:
                    qdac.ramp_multi_ch_slowly([main_gate.instrument],[unconditional_end_ramp_Vgo])
            signal.signal(signal.SIGINT, lambda sig, frame: [cleanup(), sys.exit(0)])

            qdac.add_dc_voltages_to_metadata(datasaver=datasaver)
            zurich.save_config_to_metadata(datasaver=datasaver)
            #i=0
            #for increment in zip(increments):
            #    i+=1
            #    datasaver.dataset.add_metadata(f'increments{i}',increment)
            zurich.freq0(self.freq_RLC)
            fast_axis_unreversible_list = list(inner_gate_sweep) #(to deal with snake)
            reversed_sweep=False
            i=0
            #n=0#outer sweep count

            max_sens_list=[]
            max_sens_Vcs_list=[]
            for outer_gate_value in tqdm(outer_gate_sweep, leave=False, desc='outer Gate Sweep', colour = 'green'): #slow axis loop (gate)
                i=i+1#outergatesweepcounter
                #print('temperature')
                #Triton.MC()
                outer_gate_sweep.set(outer_gate_value)
                for auxgate,increment in zip(aux_gates,increments):
                    current_outer_gate_V=auxgate()
                    #print(f"current_outer_gate_V={current_outer_gate_V}")
                    time.sleep(0.5)
                    auxgate(current_outer_gate_V+increment*step_vgo)


                time.sleep(abs(step_vgo/slew_rate)) # Wait  the time it takes for the voltage to settle - doesn't quite work! #SF FIX SLEEP TIMES!
                Glist=[]
                Vlist=[]
                Rlist=[]
                Phaselist=[]
                IsensList=[]
                
                #print(f"lb={lower_boundary},ub={upper_boundary}")
                for inner_gate_value in tqdm(inner_gate_sweep, leave=False, desc='inner gate Sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
                    if (inner_gate_value >= lower_boundary and inner_gate_value <= upper_boundary):
                        inner_gate_sweep.set(inner_gate_value)
                        time.sleep(1.1*tc+step_vgi/slew_rate) # Wait 3 times the time contanst of the lock-in plus gate ramp speed

                        theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                        sens_value=sens_demod()
                        theta_sens, v_r_sens, I_sens, _ = zurich.phase_voltage_current_conductance_compensate(vsdac,x_avg=0,y_avg=0,measured_value=sens_value)
        
                        Glist=Glist+[G]#empirical correction
                        Vlist=Vlist+[v_r_calc]
                        Phaselist=Phaselist+[theta_calc]
                        IsensList=IsensList+[I_sens]
                    else:
                        Glist=Glist+[-1e-15]
                        Vlist=Vlist+[-1e-15]
                        Phaselist=Phaselist+[-1e-15]
                        IsensList=IsensList+[-1e-15]
                #temp_fast_axis_list.reverse()
                Glist_np=np.array(Glist)
                maxid=np.argmax(Glist_np)
                IsensList_np=np.array(IsensList)
                max_sens=max(IsensList_np)
                max_sens_V_id=np.argmax(IsensList_np)
                V_of_max=list(inner_gate_sweep)[maxid]
                Vcs_of_max_sens=list(inner_gate_sweep)[max_sens_V_id]
                max_sens_list.append(max_sens)
                max_sens_Vcs_list.append(Vcs_of_max_sens)
                #print(f"maxid={maxid}")
                #print(f"V_of_max{V_of_max}")
                lower_boundary=V_of_max-scan_range/2
                upper_boundary=V_of_max+scan_range/2
                if reversed_sweep: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
                    Glist.reverse()
                    Vlist.reverse()
                    Rlist.reverse()
                    Phaselist.reverse()
                    IsensList.reverse()
                    #GIVlist.reverse()
                    #VRlist.reverse()
                    #PHASElist.reverse()
                datasaver.add_result(('G', Glist),
                                    ('V_r', Vlist),
                                    ('Phase', Phaselist),
                                    ('I_sens', IsensList),
                                    (outer_gate_sweep.parameter,outer_gate_value),
                                    (inner_gate_sweep.parameter,fast_axis_unreversible_list))
                
                
                inner_gate_sweep.reverse() 
                reversed_sweep= not reversed_sweep
        
        sine_wave_context.abort()
        RF_sens_osc(self.freq_RLC)
        datasaver.dataset.add_metadata(f'endVg1',qdac.ch01.dc_constant_V())
        datasaver.dataset.add_metadata(f'endVg2',qdac.ch02.dc_constant_V())
        datasaver.dataset.add_metadata(f'endVg3',qdac.ch03.dc_constant_V())
        datasaver.dataset.add_metadata(f'endVg4',qdac.ch04.dc_constant_V())
        datasaver.dataset.add_metadata(f'endVg5',qdac.ch05.dc_constant_V())
        datasaver.dataset.add_metadata(f'endVg6',qdac.ch06.dc_constant_V())

        maxmax_sens=max(max_sens_list)
        maxmax_sens_id=np.argmax(max_sens_list)
        print(f"maxmax_sens={maxmax_sens}")
        print(maxmax_sens_id)

        maxmax_sens_vgo=list(outer_gate_sweep)[maxmax_sens_id]
        maxmax_sens_Vgcs=max_sens_Vcs_list[maxmax_sens_id]
        #time.sleep(abs(stop_vg)/ramp_speed/1000 + 10)
        print("wake up, main gate is")
        print(main_gate())

        print(inner_gate())
        if return_max:
            return maxmax_sens_vgo,maxmax_sens_Vgcs,maxmax_sens

        ###continue
    

    def find_mech_mode(self,start_drive=75e-3,end_drive=50e-6,freq_range=None,found_range=4e6,start_step_pitch=2e3,div_factor=4,div_f=2,min_sig_I=1e-12,avg_num=5):
        zurich.output1_amp1(start_drive)
        if freq_range==None:
            start_f = self.start_f
            stop_f = self.stop_f
        else:
            start_f=freq_range[0]
            stop_f=freq_range[1]
        lowest_effective_drive=copy.copy(start_drive)
        
        step_num_f=round((-start_f+stop_f)/start_step_pitch)
        #if find_sitpos:
        #    self.sit_at_max_Isens(side="right")
        #print(f"setting drive to {start_drive}")
        zurich.output1_amp1(start_drive)
        
        I,f=self.mech_simple_fun_db(costum_prefix="find_mech_start",start_f=start_f,stop_f=stop_f,step_num_f=step_num_f,return_I_and_f=True)
        maxI_id=np.argmax(centered_moving_average(I,n=avg_num))
        f_of_max=f[maxI_id]
        intermediate_drive=start_drive/2
        intermediate_start_f=f_of_max-found_range/div_f
        intermediate_stop_f=f_of_max+found_range/div_f
        intermediate_step_num_f=round(found_range/start_step_pitch)
        intermediate_range=found_range
        
        while intermediate_drive>end_drive:
  
            
            zurich.output1_amp1(intermediate_drive)
            I,f=self.mech_simple_fun_db(costum_prefix="find_mech_intermediate",start_f=intermediate_start_f,stop_f=intermediate_stop_f,step_num_f=intermediate_step_num_f,return_I_and_f=True)
            if max(centered_moving_average(I,n=avg_num))>min_sig_I:
                maxI_id=np.argmax(centered_moving_average(I,n=avg_num))
                f_of_max=f[maxI_id]
                print(f"found mode at {f_of_max} with drive amplitude {lowest_effective_drive} ")
                lowest_effective_drive=lowest_effective_drive/2
            intermediate_drive=intermediate_drive/div_factor
            intermediate_range=intermediate_range/2
            intermediate_start_f=f_of_max-intermediate_range/2
            intermediate_stop_f=f_of_max+intermediate_range/2
            

        zurich.output1_amp1(end_drive)
        I,f=self.mech_simple_fun_db(costum_prefix="find_mech_final",start_f=intermediate_start_f,stop_f=intermediate_stop_f,step_num_f=intermediate_step_num_f,return_I_and_f=True)
        if max(centered_moving_average(I,n=avg_num))>min_sig_I:
            maxI_id=np.argmax(centered_moving_average(I,n=avg_num))
            f_of_max=f[maxI_id]
            lowest_effective_drive=end_drive

        return f_of_max,end_drive


    def linesweep_parallel_LFsens_extended(self,#in construction
                           device_name=None,
                           costum_prefix='_',
                           start_vgo =  None,#
                           stop_vgo =   None,#
                            step_vgo_num = None,
                            start_vgi = None,#-0.788
                            stop_vgi = None,#-0.776
                            step_vgi_num = None,
                            start_vgi_scan=None,#first guess for peak
                            scan_range=None,
                            mod_amplitude=0.1e-3,
                            mod_frequency=1e3,
                            increments=[0,0,0,0],
                            main_gate=qdac.ch01.dc_constant_V,
                            aux_gates=[],
                            pre_ramping_required=True,
                            load_params=True,
                            find_startpos=True,
                            check_around_current_V=True,
                            check_V_range=[-0.03,0.03],
                            check_pt_pitch=3e-3,
                            set_best_sitpos=True,#works only for single vgo!
                            sitside="right",
                            sitpos_precision_factor=5, #multiplicator for eventual sitpos determination
                            unconditional_end_ramp_Vgo=None
             ):
        if load_params:
            self.load_parameters()
        if start_vgo == None: #defined here just for initial ramp
            start_vgo = self.start_vgo_ls
        if start_vgi == None:
            start_vgi = self.start_vgi_ls
        if stop_vgi == None:
            stop_vgi = self.stop_vgi_ls
        if step_vgi_num == None:
            step_vgi_num = self.step_vgi_num_ls
        #print(check_around_current_V)

        if aux_gates=="all":
            aux_gates=[qdac.ch02.dc_constant_V,qdac.ch03.dc_constant_V,qdac.ch04.dc_constant_V,qdac.ch05.dc_constant_V]
            if increments==[0,0,0,0]:
                increments=[1,1,1,1]
        
        if check_around_current_V:
            print("checking around current V")
            start_vgo=main_gate()+check_V_range[0]
            stop_vgo=main_gate()+check_V_range[1]
            step_vgo_num=int(max([abs(round(stop_vgo-start_vgo)/check_pt_pitch),10]))
            print(f"step_vgo_num {step_vgo_num}")
            print(f"start{start_vgo} stop{stop_vgo} step_nr {step_vgo_num}")
            unconditional_end_ramp_Vgo=start_vgo

        if pre_ramping_required:
            print(f"Pre-ramping main gate to {start_vgo}")
            qdac.ramp_multi_ch_slowly(channels=[main_gate.instrument], final_vgs=[start_vgo],step_size=self.ramp_step_size,ramp_speed=self.max_ramp_speed)

        if find_startpos:
            Vg,G,_=self.GVG_fun_sensitivity(start_vg=start_vgi,stop_vg=stop_vgi,step_num=step_vgi_num, return_only_Vg_G_and_Isens=True,return_data=True)
            peakpos=Vg[np.argmax(G)]
            print(f"following highest peak detected at {peakpos*1e3:.5g} mV")
        else:
          if start_vgi_scan == None:
            start_vgi_scan = self.start_vgi_scan_ls  

        maxmax_sens_vgo,maxmax_sens_Vgcs,maxmax_sens=self.linesweep_parallel_LFsens(
                           device_name=device_name,
                           costum_prefix=costum_prefix,
                           start_vgo =  start_vgo,#
                           stop_vgo =   stop_vgo,#
                            step_vgo_num = step_vgo_num ,
                            start_vgi = start_vgi,#-0.788
                            stop_vgi = stop_vgi,#-0.776
                            step_vgi_num =step_vgi_num ,
                            start_vgi_scan=peakpos,#first guess for peak
                            scan_range=scan_range,
                            mod_amplitude=mod_amplitude,
                            mod_frequency=mod_frequency,
                            increments=increments,
                            main_gate=main_gate,
                            aux_gates=aux_gates,
                            pre_ramping_required=pre_ramping_required,
                            load_params=load_params,
                            unconditional_end_ramp_Vgo=unconditional_end_ramp_Vgo
             )
        
        if set_best_sitpos:
            qdac.ramp_multi_ch_slowly([main_gate.instrument,self.cs_gate],[maxmax_sens_vgo,start_vgi],step_size=self.ramp_step_size,ramp_speed=self.max_ramp_speed)
            sitpos=self.sit_at_max_Isens(side=sitside,start_vg=start_vgi,stop_vg=stop_vgi,step_num=step_vgi_num*sitpos_precision_factor)
            time.sleep(10)
        

        return maxmax_sens_vgo,maxmax_sens_Vgcs,maxmax_sens