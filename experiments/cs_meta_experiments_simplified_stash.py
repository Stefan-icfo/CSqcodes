#from experiments.cs_experiment import CSExperiment

from experiments.cs_experiment import CSExperiment

import experiment_parameters as params

from instruments import *
from experiments.zurich_experiments.spectrum_0925B import run_thermomech_temp_meas
import time

from utils.zurich_data_fkt import *

from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import copy
from experiments.zurich_experiments.takedemodtimetrace211025 import *




class CS_meta(CSExperiment):
    def __init__(self):
        super().__init__()
        #from experiments.zurich_experiments.spectrum_0925B import run_thermomech_temp_meas
        # Add any additional initialization here if needed
    def load_parameters(self):
    # First, call the parent's load_parameters
        super().load_parameters()
        
        # The parent has already imported and reloaded params
        import importlib
        importlib.reload(params)  
        # Add CS_meta specific parameters
        self.therm_reps = params.therm_reps
        self.temp_meas_counts = params.temp_meas_counts
        self.increment_meta= params.increment_meta
        self.startpos_gate_meta=params.startpos_gate_meta
        self.startpos_auxgate_meta=params.startpos_auxgate_meta
        self.freq_bands=params.freq_bands
        self.softening_pitch=params.softening_pitch
        self.softening_reps=params.softening_reps
        self.background_reps=params.background_reps
        self.pos_list=params.pos_list
        self.autocorr_reps=params.autocorr_reps
        self.manual_background_set=params.manual_background_set
        self.manual_thermomech_frequency=params.manual_thermomech_frequency
        self.update_therm_freq=params.update_therm_freq
        self.therm_autocorr_pitch=params.therm_autocorr_pitch
        self.autocorr_Vg_pitch=params.autocorr_Vg_pitch
        self.driven_avg_num_meta=params.driven_avg_num_meta
        self.driven_range_meta=params.driven_range_meta
        self.driven_pitch_meta=params.driven_pitch_meta
        self.driven_amp_meta=params.driven_amp_meta
        self.DQD_stability_start_vg1=params.DQD_stability_start_vg1
        self.DQD_stability_start_vg2=params.DQD_stability_start_vg2
        self.pos_listg3h2g1=params.pos_listg3h2g1
        self.mech_freq_list=params.mech_freq_list
        self.pos_list_5g_freq=params.pos_list_5g_freq
        self.cs_ranges=params.cs_ranges


    





    def therm_vs_sitpos(self,f_mech,demod_only=False,Vg_cs_adjustment_during_measurement=False):#maybe add separate nr of reps for background here...
        self.load_parameters()
        reps_nodrive=self.softening_reps
        softening_pitch=self.softening_pitch
        max_detuning=0.7e-3
        Vg,G,sens=self.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=True)
        G_avg=centered_moving_average(a=G,n=10)
        peakpos=Vg[np.argmax(G_avg)]
        print(f"initial peakpos cs {peakpos*1e3:.5g} mV")

        start_vg=peakpos-max_detuning#for non-adjustment case
        stop_vg=peakpos+max_detuning

        if Vg_cs_adjustment_during_measurement:#now adjust a narrow range for the regular GVg between thermal sweeps
            self.set_params(start_vg_cs=peakpos-5e-3)
            self.set_params(stop_vg_cs=peakpos+5e-3)
            self.set_params(step_num_cs=10*100)

        
        
        zurich.set_mixdown(f_mech-1e6)
        qdac.ramp_multi_ch_slowly([6],[start_vg])
        time.sleep(100)
        background_id=run_thermomech_temp_meas(exp_name=f'background_',reps_nodrive=reps_nodrive,background_id=None)

        
        next_detuning=-copy.copy(max_detuning)
        continue_loop_condition=True
        while continue_loop_condition:
            current_V=qdac.ch06.dc_constant_V()
            print(f"set ch06  to {current_V:6g} mV")
            time.sleep(5)
            zurich.set_mixdown(f_mech)
            time.sleep(100)
            if not demod_only:
                run_thermomech_temp_meas(exp_name=f'thermalV_gcs_={current_V*1e3:6g} mV',reps_nodrive=reps_nodrive,background_id=background_id)
            else:
                 for m in range(self.autocorr_reps):
                      takedemodtimetrace()

            if Vg_cs_adjustment_during_measurement:
                Vg,G,sens=self.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=True)
                G_avg=centered_moving_average(a=G,n=10)
                peakpos=Vg[np.argmax(G_avg)]
                next_detuning+=self.softening_pitch
                qdac.ch06.dc_constant_V(peakpos+next_detuning)
                continue_loop_condition=(next_detuning<max_detuning)
            else:
                qdac.ch06.dc_constant_V(current_V+softening_pitch)
                time.sleep(1)
                continue_loop_condition=(qdac.ch06.dc_constant_V()<stop_vg)
        qdac.ramp_multi_ch_slowly([6],[start_vg])
        time.sleep(10)
        self.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=False)#doublecheck


    def measure_singledot_config(self,
                                 thermal_spectra=True,
                                 driven_traces=False,
                                 temp_meas_counts=None,
                                 therm_reps=None,                   ##########                             
                                 thermal_softening=False,
                                 find_freq_range=None,                             
                                 background_id=1,
                                 name_addition=None,
                                 softening_demod_only=False,
                                 Vg_cs_adjustment_during_measurement=True,
                                 adjustment_linesweep=False):
        #init
        self.load_parameters()
        if therm_reps==None:
            therm_reps=self.therm_reps
        if temp_meas_counts==None:
            temp_meas_counts=self.temp_meas_counts
        if name_addition is not None:
             exp_name=f"Spectrum"+name_addition
        else:
             exp_name=f"Spectrum"
             name_addition="_"
        autocorr_reps=self.autocorr_reps

        if adjustment_linesweep:
            _,_,I_sens_sit=self.linesweep_parallel_LFsens_extended(costum_prefix='adjustment_linesweep',
                                                sitside="right",
                                                check_around_current_V=True,
                                                check_V_range=[-0.03,0.03],    ##########
                                                check_pt_pitch=2e-3,          ###########
                                                set_best_sitpos=True,
                                                find_startpos=True,
                                                main_gate=qdac.ch01.dc_constant_V)
        
        #find_mode# could do else statement and put for I_Sens_sit one of the linesweep returns
        else:
            _,I_sens_sit=self.sit_at_max_Isens(side=self.sitside)#changed evening 181025
        print(f"FINDING MECHANICAL MODE")
        f_max,_=self.find_mech_mode()
        if f_max==None:
             print("MOVING TO NEXT POS")
             return
             
        #possibly overriding frequency by manual input
        if self.manual_thermomech_frequency is not None:#override f_max
             updated_freq=self.manual_thermomech_frequency
        else:
             updated_freq=f_max

        #driven trances if desired
        if driven_traces:
            driven_avg_num_meta=self.driven_avg_num_meta
            driven_range_meta=self.driven_range_meta
            driven_pitch_meta=self.driven_pitch_meta
            driven_amp_meta=self.driven_amp_meta
            zurich.output1_amp1(driven_amp_meta)
            for n in range(driven_avg_num_meta):
                  self.mech_simple_fun_db(costum_prefix="for_avg_singledot_config_"+name_addition,start_f=updated_freq-driven_range_meta/2,stop_f=updated_freq+driven_range_meta/2,step_num_f=abs(round(driven_range_meta/driven_pitch_meta)))
            zurich.output1_amp1(0)

        zurich.set_mixdown(f_max+1e3)
        time.sleep(100)
        if thermal_spectra:
                print("THERMOMECHANICAL SPECTRUM")
                for n in range(temp_meas_counts):
                    updated_freq=run_thermomech_temp_meas(reps_nodrive=therm_reps,
                                             take_time_resolved_spectrum=True,
                                             background_id=background_id,
                                             add_to_metadata=[I_sens_sit],
                                             metadata_entry_names=['I_sens_sit'],
                                             exp_name=exp_name)
                    if self.update_therm_freq:
                        zurich.set_mixdown(updated_freq+1e3)
                    print("demod_timetraces")
                for m in range(autocorr_reps):
                        takedemodtimetrace()
             
        if thermal_softening:
                print("SOFTENING, THERMAL")
                #self.therm_vs_g2(f_mech=f_max)
                self.therm_vs_sitpos(f_mech=f_max,demod_only=softening_demod_only,Vg_cs_adjustment_during_measurement=Vg_cs_adjustment_during_measurement)

        
        return updated_freq
        


    def go_through_gate_pos_softening_only(self,pos_list=None,name_addition=None,
                            gate=qdac.ch02,auxgate=qdac.ch01,increment=-0.4,startpos_gate=0.3,startpos_auxgate=0.8,
                            pos_list_5g=True,#in this case all 5 gate positions are given in the pos list, not just one auxgate and one compensation gate; hence, the above values arent used
                            Vg_cs_adjustment_during_measurement=True):
        self.load_parameters()
        if pos_list is None:
             pos_list=self.pos_list
        for i, pos in enumerate(pos_list):
            self.load_parameters()
            if name_addition is None:     
                name_addition_full=f"step_{i+1}"
            else:
               name_addition_full=f"step_{i+1}" +name_addition
            if self.freq_bands is not None:
                     freq_bands=self.freq_bands
            zurich.sigout1_amp1_enabled_param.value(0)#switch off gate just incase it's on
            if pos_list_5g:
                qdac.ramp_multi_ch_slowly([1,2,3,4,5],pos)
            else:#auxgate-adjustment method
                auxgate_pos=startpos_auxgate+increment*(pos-startpos_gate)
                print(f"ramping to next step nr {i+1} at gate={pos} and auxgate={auxgate_pos}")
                qdac.ramp_multi_ch_slowly([gate,auxgate],[pos,auxgate_pos])
            qdac.read_channels()
            print(f"i={i},softening")
              
            for freq_band in freq_bands:
                    self.load_parameters()
                    if self.freq_bands is not None:
                        freq_bands=self.freq_bands
                
                #softening_pitch=self.softening_pitch
                #softening_reps=self.softening_reps
                
            self.measure_singledot_config(thermal_spectra=False,
                                 temp_meas_counts=0,
                                 therm_reps=0,
                                 find_freq_range=freq_band,                  ##########                             
                                 thermal_softening=True,
                                 driven_traces=False,
                                 background_id=self.manual_background_set,
                                 name_addition=name_addition_full,
                                 softening_demod_only=False,
                                 Vg_cs_adjustment_during_measurement=Vg_cs_adjustment_during_measurement)#for now only demod
                    


    def go_through_gate_pos(self,pos_list=None,name_addition=None,
                            gate=qdac.ch02,auxgate=qdac.ch01,increment=-0.4,startpos_gate=0.3,startpos_auxgate=0.8,
                            ):
        self.load_parameters()
        if pos_list is None:
             pos_list=self.pos_list
        for i, pos in enumerate(pos_list):
            self.load_parameters()
            if name_addition is None:     
                name_addition_full=f"step_{i+1}"
            else:
               name_addition_full=f"step_{i+1}" +name_addition
            
            if self.freq_bands is not None:
                     freq_bands=self.freq_bands
            #softening_pitch=self.softening_pitch
            #softening_reps=self.softening_reps
            therm_reps=self.therm_reps
            #temp_meas_count=self.temp_meas_count
            background_reps=self.background_reps
            temp_meas_counts=self.temp_meas_counts
            
            zurich.sigout1_amp1_enabled_param.value(0)#switch off gate just incase it's on
            auxgate_pos=startpos_auxgate+increment*(pos-startpos_gate)
            print(f"ramping to next step nr {i+1} at gate={pos} and auxgate={auxgate_pos}")
            time.sleep(10)
            qdac.ramp_multi_ch_slowly([gate,auxgate],[pos,auxgate_pos],step_size=4e-2,ramp_speed=4e-3)
            time.sleep(10)
            qdac.read_channels()
            softening=False
            if i % 5 == 0:
                softening=False#no softening at all

                print(f"i={i},softening={softening}")
               # if i==0:
                #     softening=False
                
                if self.manual_background_set is None:
                    print("BACKGROUND SPECTRUM")
                    self.sit_at_max_Isens(side="left")
                    zurich.set_mixdown(130e6)
                    time.sleep(100)
                    background_id=run_thermomech_temp_meas(exp_name=f"backgroundspecat_{pos}",reps_nodrive=background_reps,take_time_resolved_spectrum=True,background_id=None)
                    
                else:
                     background_id=self.manual_background_set
                #background_id=
            for freq_band in freq_bands:
                self.load_parameters()
                if self.freq_bands is not None:
                     freq_bands=self.freq_bands
                
                #softening_pitch=self.softening_pitch
                #softening_reps=self.softening_reps
                
                self.measure_singledot_config(thermal_spectra=True,
                                 temp_meas_counts=temp_meas_counts,
                                 therm_reps=therm_reps,
                                 find_freq_range=freq_band,                  ##########                             
                                 thermal_softening=softening,
                                #softening_reps=softening_reps, 
                              #softening_pitch=softening_pitch,               ##########              
                                 background_id=background_id,
                                 name_addition=name_addition_full,
                                 driven_traces=False)

    


    def movedot_g2g3(self,pos_listg3h2g1=None,name_addition=None,softening=True
                            ):
        self.load_parameters()
        if pos_listg3h2g1 is None:
             pos_listg3h2g1=self.pos_listg3h2g1
        for i, pos in enumerate(pos_listg3h2g1):
            self.load_parameters()
            if name_addition is None:     
                name_addition_full=f"step_{i+1}"
            else:
               name_addition_full=f"step_{i+1}" +name_addition
            
            if self.freq_bands is not None:
                     freq_bands=self.freq_bands
            #softening_pitch=self.softening_pitch
            #softening_reps=self.softening_reps
            therm_reps=self.therm_reps
            #temp_meas_count=self.temp_meas_count
            background_reps=self.background_reps
            temp_meas_counts=self.temp_meas_counts
            
            zurich.sigout1_amp1_enabled_param.value(0)#switch off gate just incase it's on
            qdac.ramp_multi_ch_slowly([3,2,1],pos)
            print(f"ramping to next step nr {i+1}")
            time.sleep(50)
            qdac.read_channels()
            softening=softening
            if i % 5 == 0:
               

                print(f"i={i},softening={softening}")
               # if i==0:
                #     softening=False
                
                if self.manual_background_set is None:
                    print("BACKGROUND SPECTRUM")
                    self.sit_at_max_Isens(side="left")
                    zurich.set_mixdown(130e6)
                    time.sleep(100)
                    background_id=run_thermomech_temp_meas(exp_name=f"backgroundspecat_{pos}",reps_nodrive=background_reps,take_time_resolved_spectrum=True,background_id=None)
                    
                else:
                     background_id=self.manual_background_set
                #background_id=
            for freq_band in freq_bands:
                self.load_parameters()
                if self.freq_bands is not None:
                     freq_bands=self.freq_bands
                
                #softening_pitch=self.softening_pitch
                #softening_reps=self.softening_reps
                
                self.measure_singledot_config(thermal_spectra=True,
                                 temp_meas_counts=temp_meas_counts,
                                 therm_reps=therm_reps,
                                 find_freq_range=freq_band,                  ##########                             
                                 thermal_softening=softening,
                                #softening_reps=softening_reps, 
                              #softening_pitch=softening_pitch,               ##########              
                                 background_id=background_id,
                                 name_addition=name_addition_full,
                                 driven_traces=False,
                                 adjustment_linesweep=True)


    def movedot_g2g3_with_g1_sweep(self,pos_listg3h2g1=None,mech_freq_list=None
                            ):
        #this requires the frequencies to be found precisely, separately, before the run - because it might be that the g1 start position there is nothing, so better to find it first in another code
        self.load_parameters()
        if pos_listg3h2g1 is None:
             pos_listg3h2g1=self.pos_listg3h2g1
        if mech_freq_list is None:
             mech_freq_list=self.mech_freq_list

             i=-1
        for  pos,freq in zip(pos_listg3h2g1,mech_freq_list):
            self.load_parameters()
            i+=1
                      
            zurich.sigout1_amp1_enabled_param.value(0)#switch off gate just incase it's on
            #pos[2]=pos[2]-100e-3#symmetrize around g1 value, now done in therm_vs_g1()
            qdac.ramp_multi_ch_slowly([3,2,1],pos)
            print(f"ramping to next step nr {i+1}")
            qdac.read_channels()
            
            self.therm_vs_g1(freq)
            
           
               
                
                
    
    def ramp_to_gate_pos(self,pos,
                            gate=qdac.ch02,auxgate=qdac.ch01,increment=None,startpos_gate=None,startpos_auxgate=None,
                           ):
            self.load_parameters()
            if increment==None:
                 increment=self.increment_meta
            if startpos_gate==None:
                 startpos_gate=self.startpos_gate_meta
            if startpos_auxgate==None:
                 startpos_auxgate=self.startpos_auxgate_meta
            auxgate_pos=startpos_auxgate+increment*(pos-startpos_gate)
            print(f"ramping to  gate={pos} and auxgate={auxgate_pos}")
            time.sleep(10)
            qdac.ramp_multi_ch_slowly([gate,auxgate],[pos,auxgate_pos],step_size=4e-2,ramp_speed=4e-3)
            time.sleep(10)
            self.sit_at_max_Isens(side=self.sitside)




            
            
    def therm_vs_g2(self,f_mech,reps_nodrive=50,g2_pitch=5e-3,demod_only=False,compensate_g1=False):#maybe add separate nr of reps for background here...
        self.load_parameters()
        Vg,G,sens=self.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=True)
        peakpos=Vg[np.argmax(G)]
        print(f"setting cs to {peakpos*1e3:.5g} mV")
        self.set_params(start_vg_cs=peakpos-5e-3)
        self.set_params(stop_vg_cs=peakpos+5e-3)
        self.set_params(step_num_cs=10*50)
        #approx_maxpos=peakpos
        start_vg2=qdac.ch02.dc_constant_V()
        stop_vg2=start_vg2+200e-3
        mech_freq=f_mech
        
        self.sit_at_max_Isens(side="left")
        zurich.set_mixdown(mech_freq-1e6)
        time.sleep(100)
        background_id=run_thermomech_temp_meas(exp_name=f'background_',reps_nodrive=reps_nodrive,background_id=None)#and here...and in calling fkt
        time.sleep(100)
        i=0
        if not demod_only:
            while qdac.ch02.dc_constant_V()<stop_vg2:
                i+=1
                zurich.set_mixdown(mech_freq)
                current_Vg2=qdac.ch02.dc_constant_V()
                current_Vg1=qdac.ch01.dc_constant_V()
                time.sleep(1)
                qdac.ch02.dc_constant_V(current_Vg2+g2_pitch)
                if compensate_g1:
                    qdac.ch01.dc_constant_V(current_Vg1-0.4*g2_pitch)
                time.sleep(5)
                self.sit_at_max_Isens(side="left")
                zurich.set_mixdown(mech_freq)
                time.sleep(100)
                run_thermomech_temp_meas(exp_name=f'g2_thermalV_gcs_={current_Vg2*1e3:6g} mV',reps_nodrive=reps_nodrive,background_id=background_id)
                print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
                zurich.set_mixdown(mech_freq)
                print(f"setting ch02  to {current_Vg2:6g} mV")
                time.sleep(5)
        self.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=False)#doublecheck


    def therm_vs_g1(self,f_mech,reps_nodrive=50,g1_pitch=10e-3,g1_range=198e-3):#maybe add separate nr of reps for background here...
        self.load_parameters()
        current_Vg1=qdac.ch01.dc_constant_V()
        time.sleep(1)
        qdac.ch01.dc_constant_V(current_Vg1-g1_range/2)
        time.sleep(20)

        Vg,G,sens=self.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=True)
        peakpos=Vg[np.argmax(G)]
        print(f"setting cs to {peakpos*1e3:.5g} mV")
        self.set_params(start_vg_cs=peakpos-5e-3)
        self.set_params(stop_vg_cs=peakpos+5e-3)
        self.set_params(step_num_cs=10*50)
        #approx_maxpos=peakpos
        start_vg1=qdac.ch01.dc_constant_V()
        stop_vg1=start_vg1+g1_range
        mech_freq=f_mech
        zurich.set_mixdown(mech_freq-1e6)
        self.sit_at_max_Isens(side="left")
        time.sleep(100)
        background_id=run_thermomech_temp_meas(exp_name=f'background_',reps_nodrive=reps_nodrive,background_id=None)#and here...and in calling fkt
        time.sleep(100)
        i=0
        
        while qdac.ch01.dc_constant_V()<stop_vg1:
                i+=1
                zurich.set_mixdown(mech_freq)
                current_Vg1=qdac.ch01.dc_constant_V()
                time.sleep(1)
                qdac.ch01.dc_constant_V(current_Vg1+g1_pitch)
                time.sleep(5)
                self.sit_at_max_Isens(side="left")
                zurich.set_mixdown(mech_freq)
                time.sleep(100)
                run_thermomech_temp_meas(exp_name=f'g1_thermalV_gcs_={current_Vg1*1e3:6g} mV',reps_nodrive=reps_nodrive,background_id=background_id)
                print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
                zurich.set_mixdown(mech_freq)
                print(f"setting ch01  to {current_Vg1:6g} mV")
                time.sleep(5)
        self.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=False)#doublecheck


######################
    def find_freq_only_5g(self,first_freq_guess,freq_range,name_addition=None,
                            pos_list_5g_freq=None):
        self.load_parameters()
        if pos_list_5g_freq is None:
             pos_list_5g_freq=self.pos_list_5g_freq
        #if cs_range is None:
             #cs_range=self.cs_range
        self.set_params(start_f=first_freq_guess-freq_range/2)
        self.set_params(stop_f=first_freq_guess+freq_range/2)#pitch from find_M
              
        for i, pos in enumerate(pos_list_5g_freq):
            self.load_parameters()
            if name_addition is None:     
                name_addition_full=f"step_{i+1}"
            else:
               name_addition_full=f"step_{i+1}" +name_addition
            
            zurich.sigout1_amp1_enabled_param.value(0)#switch off gate just incase it's on
            qdac.ramp_multi_ch_slowly([1,2,3,4,5],pos)
            #qdac.ramp_multi_ch_slowly([],cs_range[0])
            
            qdac.read_channels()
            
            
            self.load_parameters()
                    
                
                #softening_pitch=self.softening_pitch
                #softening_reps=self.softening_reps
                
            new_freq=self.measure_singledot_config(thermal_spectra=False,
                                 temp_meas_counts=0,
                                 therm_reps=0,
                                 find_freq_range=None,                  ##########                             
                                 thermal_softening=False,
                                 driven_traces=False,
                                 background_id=self.manual_background_set,
                                 name_addition=name_addition_full,
                                 softening_demod_only=False,
                                 )#for now only demod
            self.set_params(start_f=new_freq-freq_range/2)
            self.set_params(stop_f=new_freq+freq_range/2)#pitch from find_M



##################################
    def find_freq_only_cs(self,first_freq_guess,freq_range,name_addition=None,
                            cs_ranges=None):
        self.load_parameters()
        #if pos_list_5g_freq is None:
        #     pos_list_5g_freq=self.pos_list_5g_freq
        if cs_ranges is None:
             cs_ranges=self.cs_ranges
        self.set_params(start_f=first_freq_guess-freq_range/2)
        self.set_params(stop_f=first_freq_guess+freq_range/2)#pitch from find_M
              
        for i, cs_pos in enumerate(cs_ranges):
            self.load_parameters()
            if name_addition is None:     
                name_addition_full=f"step_{i+1}"
            else:
               name_addition_full=f"step_{i+1}" +name_addition
    
            zurich.sigout1_amp1_enabled_param.value(0)#switch off gate just incase it's on
            qdac.ramp_multi_ch_slowly([6],cs_pos[0])
            #qdac.ramp_multi_ch_slowly([],cs_range[0])
            self.set_params(start_vg_cs = cs_pos[0],
                            stop_vg_cs = cs_pos[1]
                            )
            
            qdac.read_channels()
            
            
            self.load_parameters()
                    
                
                #softening_pitch=self.softening_pitch
                #softening_reps=self.softening_reps
                
            new_freq=self.measure_singledot_config(thermal_spectra=False,
                                 temp_meas_counts=0,
                                 therm_reps=0,
                                 find_freq_range=None,                  ##########                             
                                 thermal_softening=False,
                                 driven_traces=False,
                                 background_id=self.manual_background_set,
                                 name_addition=name_addition_full,
                                 softening_demod_only=False,
                                 )#for now only demod
            self.set_params(start_f=new_freq-freq_range/2)
            self.set_params(stop_f=new_freq+freq_range/2)#pitch from find_M




#############################
    def repeat_linesweep(self,run_id,gate_nr=2, auxgate_nr=1,step_nr=430*5,adjust_constant_gates=False):
         metadata=get_metadata(run_id,return_data=True)
         time.sleep(5)
         if metadata is None:
            print(f"Error: No metadata found for run_id {run_id}")
            return
         start_Vg=metadata[f'qdac_ch0{gate_nr}_dc_constant_V']
         stop_Vg=metadata[f'endVg{gate_nr}']
         aux_startVg=metadata[f'qdac_ch0{auxgate_nr}_dc_constant_V']
         aux_stopVg=metadata[f'endVg{auxgate_nr}']

         if adjust_constant_gates:
              constant_gate_nrs=[1,2,3,4,5]
              start_voltages=[]
              for constant_gate_nr in constant_gate_nrs:
                   start_voltages.append(metadata[f'qdac_ch0{constant_gate_nr}_dc_constant_V'])
              qdac.ramp_multi_ch_slowly(constant_gate_nrs, start_voltages, step_size=5e-2, ramp_speed=5e-3) 
                   


         self.linesweep_parallel_LFsens_extended(costum_prefix=f'rep_{run_id}',
                                                  main_gate=qdac.channel(gate_nr).dc_constant_V,
                                                  aux_gates=[qdac.channel(auxgate_nr).dc_constant_V,],
                                                  aux_gate_start_stop=[aux_startVg,aux_stopVg],
                                                  start_vgo=start_Vg,stop_vgo=stop_Vg,
                                                  step_vgo_num=step_nr,
                                                  check_around_current_V=False,
                                                  set_best_sitpos=False)
         












