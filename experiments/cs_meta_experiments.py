#from experiments.cs_experiment import CSExperiment

from experiments.cs_experiment import CSExperiment

import experiment_parameters as params

from instruments import *
from experiments.zurich_experiments.spectrum_0925B import run_thermomech_temp_meas
import time

from utils.zurich_data_fkt import *

from dataprocessing.extract_fkts import *

lower_sens_limit=4e-12

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

    def therm_driven_spectra(self,f_mech,reps_nodrive=10):
        Temp=0.16
                          
        #avg_num = 14035

        # Drive loop parameters
        drive = 50e-6      # Start at 10 μV                ######################
        end_drive = 75e-3   # End at 30 mV                ########################
        factor = 1.2      # +100% increment per step

        self.sit_at_max_Isens(side="left")
        #qdac.ch06.dc_constant_V(0.966)#setting to the other side

        zurich.set_mixdown(f_mech)
        time.sleep(50)

        while drive <= end_drive:
            zurich.output1_amp1(drive)
            time.sleep(5)  # Let system settle

            drive_uV = int(round(drive * 1e6))
            exp_name = f"Spectrum_{Temp:.3f}__{drive_uV}u"

            print(f"\n▶️ Running measurement at {drive_uV} μV drive...")

            run_thermomech_temp_meas(
                reps_nodrive=reps_nodrive,
                exp_name=exp_name,
                fit_lorentzian=False
            )
            
            drive *= factor

            exp_name = f"Spectrum_{Temp:.3f}_thermal"
            zurich.output1_amp1(0)
            run_thermomech_temp_meas(
                    reps_nodrive=reps_nodrive,
                    exp_name=exp_name,
                    fit_lorentzian=False
                )
            zurich.move_mixdown(-8e6)
            exp_name = f"Spectrum_{Temp:.3f}_background"
            zurich.output1_amp1(0)
            run_thermomech_temp_meas(
                    reps_nodrive=reps_nodrive,
                    exp_name=exp_name,
                    fit_lorentzian=False
                )
            
    def power_sweep(self,f_mech):
        amplitude_V =20e-3                                #####################
        amplitude_min_V = 50e-6#0.1e-3
        scale_factor = 0.9#<1!                           ###############
        #scale_substract=20e-6
        #reps=10
        #initial_reps=5
        step_index = 0
        self.sit_at_max_Isens(side="left")


        while amplitude_V >= amplitude_min_V:
            print(f"\n==============================")
            print(f" Step {step_index}: Measuring at amplitude = {amplitude_V * 1e3:.2f} mV")
            print(f"==============================\n")

            # Set amplitude using your custom Zurich driver
            zurich.output1_amp1(amplitude_V)

            # Optional: short delay
            time.sleep(0.5)

            # Run measurement
            
            
            #for n in range(reps):
            exp.mech_simple_fun_db(costum_prefix='g2_dot_g2_drivepwsweep',start_f=f_mech-100e3,stop_f=f_mech+1e6,step_num_f=(1.1e6/250)) ######################
            first_step=False
            amplitude_V *= scale_factor
            #amplitude_V -= scale_substract
            step_index += 1
        self.sit_at_max_Isens(side="left")


    def measure_singledot_config(self,
                                 thermal_spectra=True,
                                 temp_meas_counts=3,
                                 therm_reps=50,                   ##########
                                 driven_spectra=True,               
                                 driven_reps=5,                   #########                                  
                                 thermal_softening=True,
                                 softening_reps=9,                ##########              
                                 power_sweep=True):
        if therm_reps==None:
            therm_reps=self.therm_reps
        if temp_meas_counts==None:
            temp_meas_counts=self.temp_meas_counts

        maxmax_sens_vgo,maxmax_sens_Vgcs,maxmax_sens=self.linesweep_parallel_LFsens_extended(costum_prefix='adjustment_linesweep',
                                                sitside="right",
                                                check_around_current_V=True,
                                                check_V_range=[-0.02,0.02],    ##########
                                                check_pt_pitch=2e-3,          ###########
                                                set_best_sitpos=True,
                                                find_startpos=True,
                                                main_gate=qdac.ch02.dc_constant_V)

        #self.sit_at_max_Isens(side="left")
        if maxmax_sens>lower_sens_limit:
            print("FINDING MECHANICAL MODE")
            f_max,_=self.find_mech_mode()
            zurich.set_mixdown(f_max+5e3)
            if thermal_spectra:
                print("THERMOMECHANICAL SPECTRUM")
                for n in range(temp_meas_counts):
                    run_thermomech_temp_meas(reps_nodrive=therm_reps,take_time_resolved_spectrum=True)
            
            #if driven_spectra:
            #    print("SPECTRUM VS DRIVE")
            #    self.therm_driven_spectra(self,f_mech=f_max)
            #if thermal_softening:
            #    print("SOFTENING, THERMAL")
            #    self.therm_vs_sitpos(self,f_mech=f_max,reps_nodrive=softening_reps)
            #if power_sweep:
            #    print("POWER SWEEP")
            #    self.power_sweep(f_mech=f_max)

        else:
            print(f"sensitivity too low with maxmax_sens={maxmax_sens}")

    pos_list=[2.8,2.55,2.26,1.98,1.7,1.42,1.15,0.86,0.57,0.32] 
                                 #################
   # pos_list2=[2.8,2.7,2.6,2.5,2.4,2.3,2.18,2.08,1.95,1.85,1.74,1.64,1.51,1.41,1.3,1.2,1.06,0.97]

    def go_through_gate_pos(self,pos_list=pos_list,gate=qdac.ch02,auxgate=qdac.ch01,increment=-0.4,startpos_gate=2,startpos_auxgate=0.8-0.4*(1.7)):
        for pos in pos_list:
            auxgate_pos=startpos_auxgate+increment*(pos-startpos_gate)
            print(f"ramping to next step at gate={pos} and auxgate={auxgate_pos}")
            time.sleep(10)
            qdac.ramp_multi_ch_slowly([gate,auxgate],[pos,auxgate_pos],step_size=4e-2,ramp_speed=4e-3)
            time.sleep(10)
            self.measure_singledot_config()
            

            





















    def measure_singledot_slice(self,temp_meas_counts=3,
                                 therm_reps=20,                   ##########
                                 driven_spectra=True,               
                                 driven_reps=5,                   #########                                  
                                 thermal_softening=True,
                                 softening_reps=9,                ##########              
                                 power_sweep=True):
        self.sit_at_max_Isens(side='left')
        
        #self.sit_at_max_Isens(side="left")
        print("FINDING MECHANICAL MODE")
        f_max,_=self.find_mech_mode()
        zurich.set_mixdown(f_max+5e3)
        print("THERMOMECHANICAL SPECTRUM")
        for n in range(temp_meas_counts):
            run_thermomech_temp_meas(reps_nodrive=therm_reps,take_time_resolved_spectrum=True)
        #if driven_spectra:
        #    print("SPECTRUM VS DRIVE")
        #    self.therm_driven_spectra(self,f_mech=f_max)
        #if thermal_softening:
        #    print("SOFTENING, THERMAL")
        #    self.therm_vs_sitpos(self,f_mech=f_max,reps_nodrive=softening_reps)
        #if power_sweep:
        #    print("POWER SWEEP")
        #    self.power_sweep(f_mech=f_max)
    def go_through_gate_slices(self,start_V=-2.91,stop_V=-3.71,step_V=10e-3,gate=qdac.ch02):
        current_V=start_V
        qdac.ramp_multi_ch_slowly([gate],[current_V],step_size=9e-2,ramp_speed=9e-3)
        #gate.dc_constant_V(current_V)
        while current_V<stop_V:
            
            print(f"ramping to next step at {current_V}")
            gate.dc_constant_V(current_V)
            time.sleep(10)
            self.measure_singledot_slice()
            current_V+=step_V