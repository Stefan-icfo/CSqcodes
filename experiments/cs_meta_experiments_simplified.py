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
        


    


    def therm_vs_sitpos(self,f_mech,reps_nodrive=10,softening_pitch=0.5e-4):#maybe add separate nr of reps for background here...
        Vg,G,sens=self.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=True)
        peakpos=Vg[np.argmax(G)]
        print(f"setting cs to {peakpos*1e3:.5g} mV")
        approx_maxpos=peakpos
        start_vg=approx_maxpos-1.5e-3
        stop_vg=approx_maxpos+1.5e-3
        mech_freq=f_mech
        zurich.set_mixdown(mech_freq-1e6)
        qdac.ch06.dc_constant_V(start_vg)
        background_id=run_thermomech_temp_meas(exp_name=f'background_',reps_nodrive=reps_nodrive,background_id=None)#and here...and in calling fkt
        time.sleep(100)
        while qdac.ch06.dc_constant_V()<stop_vg:
            zurich.set_mixdown(mech_freq)
            current_V=qdac.ch06.dc_constant_V()
            qdac.ch06.dc_constant_V(current_V+softening_pitch)
            run_thermomech_temp_meas(exp_name=f'thermalV_gcs_={current_V*1e3:6g} mV',reps_nodrive=reps_nodrive,background_id=background_id)
            print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
            zurich.set_mixdown(mech_freq)
            print(f"setting ch06  to {current_V:6g} mV")
            time.sleep(5)


    def measure_singledot_config(self,
                                 thermal_spectra=True,
                                 temp_meas_counts=3,
                                 therm_reps=10,                   ##########                             
                                 thermal_softening=False,
                                 softening_reps=5,
                                 softening_pitch=1e-4, 
                                 find_freq_range=None,                             
                                 background_id=2,
                                 name_addition=None):
        self.load_parameters()
        if therm_reps==None:
            therm_reps=self.therm_reps
        if temp_meas_counts==None:
            temp_meas_counts=self.temp_meas_counts
        if name_addition is not None:
             exp_name=f"Spectrum 185mK"+name_addition
        else: exp_name=f"Spectrum 185mK"
        autocorr_reps=self.autocorr_reps
        
        _,I_sens_sit=self.sit_at_max_Isens(side=self.sitside)#changed evening 181025
        print("FINDING MECHANICAL MODE")
        f_max,_=self.find_mech_mode(start_drive=75e-3,end_drive=200e-6,freq_range=find_freq_range,found_range=1e6,start_step_pitch=None,div_factor=4,div_f=2,
                                    min_sig_I=None,#1.5e-12,
                                    min_initial_sig_I=None,#1.5e-12,
                                    avg_num=None#1
                                    )
        if f_max==None:
             print("MOVING TO NEXT POS")
             return
             
        zurich.set_mixdown(f_max+1e3)
        if thermal_spectra:
                print("THERMOMECHANICAL SPECTRUM")
                for n in range(temp_meas_counts):
                    updated_freq=run_thermomech_temp_meas(reps_nodrive=therm_reps,
                                             take_time_resolved_spectrum=True,
                                             background_id=background_id,
                                             add_to_metadata=[I_sens_sit],
                                             metadata_entry_names=['I_sens_sit'],
                                             exp_name=exp_name)
                    zurich.set_mixdown(updated_freq+1e3)
                    print("demod_timetraces")
                for m in range(autocorr_reps):
                        takedemodtimetrace()
        if thermal_softening:
                print("SOFTENING, THERMAL")
                self.therm_vs_sitpos(f_mech=updated_freq,reps_nodrive=softening_reps,softening_pitch=softening_pitch)
        


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
                softening=False

                print(f"i={i},softening={softening}")
               # if i==0:
                #     softening=False
                print("BACKGROUND SPECTRUM")
                self.sit_at_max_Isens(side="left")
                zurich.set_mixdown(120e6)
                time.sleep(100)
                background_id=run_thermomech_temp_meas(exp_name=f"backgroundspecat_{pos}",reps_nodrive=background_reps,take_time_resolved_spectrum=True,background_id=None)
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
                                 name_addition=name_addition_full)

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




            
            
    def go_along_line(self,start_gate,stop_gate,start_auxgate,stop_auxgate,
                            step_num,
                            gate=qdac.ch01,auxgate=qdac.ch02,
                            background_reps=80,
                            therm_reps=40,temp_meas_counts=3,
                            #softening_pitch=1e-4,
                            #softening_reps=1
                            ):
        self.load_parameters()
        gate_pos=np.linspace(start_gate,stop_gate,step_num)
        auxgate_pos=np.linspace(start_auxgate,stop_auxgate,step_num)
        i=0
        for gateV,auxgateV in zip(gate_pos,auxgate_pos):

            print(f"ramping to next step nr {i} at gate={gateV} and auxgate={auxgateV}")
            time.sleep(10)
            qdac.ramp_multi_ch_slowly([gate,auxgate],[gateV,auxgateV],step_size=4e-2,ramp_speed=4e-3)
            time.sleep(10)
            if i % 5 == 0:
                 #softening=True
                 print("BACKGROUND SPECTRUM")
                 self.sit_at_max_Isens(side="left")
                 zurich.set_mixdown(120e6)
                 background_id=run_thermomech_temp_meas(exp_name=f"backgroundspecat_{gateV}",reps_nodrive=background_reps,take_time_resolved_spectrum=True,background_id=None)
            self.measure_singledot_config(thermal_spectra=True,
                                 temp_meas_counts=temp_meas_counts,
                                 therm_reps=therm_reps,                   ##########                             
                                 #thermal_softening=softening,
                                #softening_reps=softening_reps, 
                             # softening_pitch=softening_pitch,               ##########              
                                 background_id=background_id)

#############################
    def repeat_linesweep(self,run_id,gate_nr=2, auxgate_nr=1,step_nr=100,adjust_constant_gates=False):
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
         












