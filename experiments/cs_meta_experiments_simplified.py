#from experiments.cs_experiment import CSExperiment

from experiments.cs_experiment import CSExperiment

import experiment_parameters as params

from instruments import *
from experiments.zurich_experiments.spectrum_0925B import run_thermomech_temp_meas
import time

from utils.zurich_data_fkt import *

from dataprocessing.extract_fkts import *



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


    


    def therm_vs_sitpos(self,f_mech,reps_nodrive=10,softening_pitch=0.5e-4):
        Vg,G,sens=self.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=True)
        peakpos=Vg[np.argmax(G)]
        print(f"setting cs to {peakpos*1e3:.5g} mV")
        approx_maxpos=peakpos
        start_vg=approx_maxpos-1.5e-3
        stop_vg=approx_maxpos+1.5e-3
        mech_freq=f_mech
        zurich.set_mixdown(mech_freq-1e6)
        qdac.ch06.dc_constant_V(start_vg)
        run_thermomech_temp_meas(exp_name=f'background_')
        time.sleep(100)
        while qdac.ch06.dc_constant_V()<stop_vg:
            zurich.set_mixdown(mech_freq)
            current_V=qdac.ch06.dc_constant_V()
            qdac.ch06.dc_constant_V(current_V+softening_pitch)
            run_thermomech_temp_meas(exp_name=f'thermalV_gcs_={current_V*1e3:6g} mV',reps_nodrive=reps_nodrive)
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
                                 softening_pitch=1e-4,                ##########              
                                 background_id=2):
        if therm_reps==None:
            therm_reps=self.therm_reps
        if temp_meas_counts==None:
            temp_meas_counts=self.temp_meas_counts

        
        _,I_sens_sit=self.sit_at_max_Isens(side="left")#changed evening 181025
        print("FINDING MECHANICAL MODE")
        f_max,_=self.find_mech_mode(start_drive=75e-3,end_drive=200e-6,freq_range=None,found_range=1e6,start_step_pitch=0.5e3,div_factor=4,div_f=2,min_sig_I=1.5e-12,min_initial_sig_I=2e-12,avg_num=1)
        if f_max==None:
             ("MOVING TO NEXT POS")
             return
             
        zurich.set_mixdown(f_max+1e3)
        if thermal_spectra:
                print("THERMOMECHANICAL SPECTRUM")
                for n in range(temp_meas_counts):
                    updated_freq=run_thermomech_temp_meas(reps_nodrive=therm_reps,
                                             take_time_resolved_spectrum=True,
                                             background_id=background_id,
                                             add_to_metadata=[I_sens_sit],
                                             metadata_entry_names=['I_sens_sit'])
                    zurich.set_mixdown(updated_freq+1e3)
        #if thermal_softening:
        #        print("SOFTENING, THERMAL")
        #        self.therm_vs_sitpos(f_mech=updated_freq,reps_nodrive=softening_reps,softening_pitch=softening_pitch)
        


        

    #midpoints_old_sweep_1946 = [0.405507, 0.555695, 0.685857, 0.821026, 0.966208, 1.11139, 1.25657, 1.39675, 1.54193, 1.68711, 1.83229, 1.97747, 2.11264, 2.25282, 2.398, 2.53317, 2.66834, 2.80851, 2.94869, 3.08385, 3.21402, 3.34919, 3.48436, 3.60951, 3.73467, 3.86984]
    
    #pos_list = [0.416102, 0.571977, 0.702712, 0.838475, 0.984294, 1.13011, 1.27593, 1.41672, 1.55751, 1.70333, 1.84915]
    pos_list1 = [0.412826, 0.563126, 0.693387, 0.828657, 0.973948, 1.11924, 1.26453, 1.40481, 1.5501, 1.69539, 1.84068, 1.98597, 2.12124, 2.26152, 2.4018, 2.53707, 2.67234, 2.81263, 2.95291, 3.08818, 3.22345, 3.35872, 3.49399, 3.61924, 3.74449, 3.87976]#midpoints of g2 linesweep on 1710
    #choosing some reliable configs
    pos_list1_simple = [0.828657, 0.973948, 1.11924, 1.26453, 1.40481, 1.5501, 1.69539, 1.84068, 1.98597, 2.12124, 2.26152, 2.4018, 2.53707, 2.67234, 2.81263, 2.95291, 3.08818, 3.22345, 3.35872, 3.49399, 3.61924, 3.74449, 3.87976]#midpoints of g2 linesweep on 1710
    #pos_list1 = [ 0.693387, 0.828657, 0.973948, 1.11924, 1.26453, 1.40481, 1.5501, 1.69539, 1.84068, 1.98597, 2.12124, 2.26152, 2.4018, 2.53707, 2.67234, 2.81263, 2.95291, 3.08818, 3.22345, 3.35872, 3.49399, 3.61924, 3.74449, 3.87976]#midpoints of g2 linesweep on 1710
    pos_list2 =[0.452906, 0.598196, 0.723447, 0.866232, 1.00902, 1.15681, 1.2996, 1.43988, 1.58768, 1.73046, 1.87826, 2.02104, 2.15381, 2.2991, 2.43437, 2.57214, 2.70491, 2.8502, 2.98547, 3.12325, 3.25601, 3.39379, 3.52655, 3.6493, 3.77705, 3.91483]#closer to top asymmetry=3
    pos_list2_simple =[ 0.866232, 1.00902, 1.15681, 1.2996, 1.43988, 1.58768, 1.73046, 1.87826, 2.02104, 2.15381, 2.2991, 2.43437, 2.57214, 2.70491, 2.8502, 2.98547, 3.12325, 3.25601, 3.39379, 3.52655, 3.6493, 3.77705, 3.91483]#closer to top asymmetry=3
    pos_list2 =[ 0.723447, 0.866232, 1.00902, 1.15681, 1.2996, 1.43988, 1.58768, 1.73046, 1.87826, 2.02104, 2.15381, 2.2991, 2.43437, 2.57214, 2.70491, 2.8502, 2.98547, 3.12325, 3.25601, 3.39379, 3.52655, 3.6493, 3.77705, 3.91483]#closer to top asymmetry=3 removed first 2
    pos_list2 =[ 0.723447, 0.866232, 1.00902, 1.15681, 1.2996, 1.43988, 1.58768, 1.73046, 1.87826, 2.02104, 2.15381, 2.2991, 2.43437, 2.57214, 2.70491, 2.8502, 2.98547, 3.12325, 3.25601, 3.39379, 3.52655, 3.6493, 3.77705, 3.91483]#closer to top asymmetry=3 removed first 2

    pos_list1 = [ 1.84068, 1.98597, 2.12124, 2.26152, 2.4018, 2.53707, 2.67234, 2.81263, 2.95291, 3.08818, 3.22345, 3.35872, 3.49399, 3.61924, 3.74449, 3.87976,0.412826, 0.563126, 0.693387]#midpoints of g2 linesweep on 1710

    pos_list1_remaining = [ 2.95291, 3.08818, 3.22345, 3.35872, 3.49399, 3.61924, 3.74449, 3.87976,0.412826, 0.563126, 0.693387]#midpoints of g2 linesweep on 1710

    pos_list1_softening = [ 1.69539,0.973948,2.4018,3.08818]#midpoints of g2 linesweep on 1710
   
    pos_list1 = [0.394398, 0.541463, 0.667519, 0.808582, 0.956647, 1.10071, 1.24678, 1.38884, 1.5329, 1.67697, 1.82103, 1.96309, 2.10116, 2.24322, 2.38628, 2.52234, 2.6574, 2.79847, 2.93853, 3.07159, 3.20265, 3.33871, 3.47377, 3.60082, 3.72688, 3.85894]

    def go_through_gate_pos(self,pos_list=pos_list1,
                            gate=qdac.ch02,auxgate=qdac.ch01,increment=-0.4,startpos_gate=4,startpos_auxgate=-0.67,
                            background_reps=80,
                            therm_reps=40,temp_meas_counts=3,
                            softening_pitch=1e-4,softening_reps=1):
        for i, pos in enumerate(pos_list):
            auxgate_pos=startpos_auxgate+increment*(pos-startpos_gate)
            print(f"ramping to next step nr {i} at gate={pos} and auxgate={auxgate_pos}")
            time.sleep(10)
            qdac.ramp_multi_ch_slowly([gate,auxgate],[pos,auxgate_pos],step_size=4e-2,ramp_speed=4e-3)
            time.sleep(10)
            if i % 5 == 0:
                 softening=True
                 print("BACKGROUND SPECTRUM")
                 self.sit_at_max_Isens(side="right")
                 zurich.set_mixdown(120e6)
                 background_id=run_thermomech_temp_meas(exp_name=f"spectrumvsenrg2at_{pos}",reps_nodrive=background_reps,take_time_resolved_spectrum=True,background_id=None)
            else:
                 softening=False
            if i==0:
                 softening=False
            self.measure_singledot_config(thermal_spectra=True,
                                 temp_meas_counts=temp_meas_counts,
                                 therm_reps=therm_reps,                   ##########                             
                                 thermal_softening=softening,
                                softening_reps=softening_reps, 
                              softening_pitch=softening_pitch,               ##########              
                                 background_id=background_id)
            
            
    def go_through_gate_pos_softening_only(self,pos_list=pos_list1_softening,
                            gate=qdac.ch02,auxgate=qdac.ch01,increment=-0.4,startpos_gate=4,startpos_auxgate=-0.67,
                            softening_pitch=1e-4,softening_reps=1):
        for i, pos in enumerate(pos_list):
            auxgate_pos=startpos_auxgate+increment*(pos-startpos_gate)
            print(f"ramping to next step nr {i} at gate={pos} and auxgate={auxgate_pos}")
            time.sleep(10)
            qdac.ramp_multi_ch_slowly([gate,auxgate],[pos,auxgate_pos],step_size=4e-2,ramp_speed=4e-3)
            time.sleep(10)
            _,I_sens_sit=self.sit_at_max_Isens(side="left")#changed evening 181025
            print("FINDING MECHANICAL MODE")
            f_max,_=self.find_mech_mode(start_drive=75e-3,end_drive=200e-6,freq_range=None,found_range=1e6,start_step_pitch=0.5e3,div_factor=4,div_f=2,min_sig_I=1.5e-12,min_initial_sig_I=2e-12,avg_num=1)
            if f_max==None:
                ("MOVING TO NEXT POS")
                return
            zurich.set_mixdown(f_max+1e3)
            updated_freq=run_thermomech_temp_meas(reps_nodrive=20,
                                             take_time_resolved_spectrum=True,
                                             background_id=654,
                                             add_to_metadata=[I_sens_sit],
                                             metadata_entry_names=['I_sens_sit'])
            zurich.set_mixdown(updated_freq)
            print("SOFTENING, THERMAL")
            self.therm_vs_sitpos(f_mech=updated_freq,reps_nodrive=softening_reps,softening_pitch=softening_pitch)

            
            














