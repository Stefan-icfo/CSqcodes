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


    


    def therm_vs_sitpos(self,f_mech,reps_nodrive=10):
        Vg,G,sens=exp.GVG_fun_sensitivity(return_only_Vg_G_and_Isens=True,return_data=True)
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
            qdac.ch06.dc_constant_V(current_V+0.5e-4)
            run_thermomech_temp_meas(exp_name=f'thermalV_gcs_={current_V*1e3:6g} mV',reps_nodrive=reps_nodrive)
            print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
            zurich.set_mixdown(mech_freq)
            print(f"setting ch06  to {current_V:6g} mV")
            time.sleep(5)


    def measure_singledot_config(self,
                                 thermal_spectra=True,
                                 temp_meas_counts=3,
                                 therm_reps=40,                   ##########
                                 driven_spectra=True,               
                                 driven_reps=5,                   #########                                  
                                 thermal_softening=True,
                                 softening_reps=9,                ##########              
                                 power_sweep=True):
        if therm_reps==None:
            therm_reps=self.therm_reps
        if temp_meas_counts==None:
            temp_meas_counts=self.temp_meas_counts

        
        self.sit_at_max_Isens(side="right")
        print("FINDING MECHANICAL MODE")
        f_max,_=self.find_mech_mode()
        zurich.set_mixdown(f_max+5e3)
        if thermal_spectra:
                print("THERMOMECHANICAL SPECTRUM")
                for n in range(temp_meas_counts):
                    run_thermomech_temp_meas(reps_nodrive=therm_reps,take_time_resolved_spectrum=True)
            

        

    midpoints_old_sweep_1946 = [0.405507, 0.555695, 0.685857, 0.821026, 0.966208, 1.11139, 1.25657, 1.39675, 1.54193, 1.68711, 1.83229, 1.97747, 2.11264, 2.25282, 2.398, 2.53317, 2.66834, 2.80851, 2.94869, 3.08385, 3.21402, 3.34919, 3.48436, 3.60951, 3.73467, 3.86984]

    pos_list = [0.416102, 0.571977, 0.702712, 0.838475, 0.984294, 1.13011, 1.27593, 1.41672, 1.55751, 1.70333, 1.84915]

    def go_through_gate_pos(self,pos_list=pos_list,gate=qdac.ch02,auxgate=qdac.ch01,increment=-0.4,startpos_gate=0.22,startpos_auxgate=0.828):
        for pos in pos_list:
            auxgate_pos=startpos_auxgate+increment*(pos-startpos_gate)
            print(f"ramping to next step at gate={pos} and auxgate={auxgate_pos}")
            time.sleep(10)
            qdac.ramp_multi_ch_slowly([gate,auxgate],[pos,auxgate_pos],step_size=4e-2,ramp_speed=4e-3)
            time.sleep(10)
            self.measure_singledot_config()
            














