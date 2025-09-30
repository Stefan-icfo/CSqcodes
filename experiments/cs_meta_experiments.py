#from experiments.cs_experiment import CSExperiment

from experiments.cs_experiment import CSExperiment

from instruments import *
from experiments.zurich_experiments.spectrum_0925 import run_thermomech_temp_meas
import time

class CS_meta(CSExperiment):
    def __init__(self):
        super().__init__()
        # Add any additional initialization here if needed


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
            #mech_freq=exp.max_thermomech_freq#here setting for real
            print(f"setting drive to thermal max:{mech_freq/1e6:6g} MHz")
            zurich.set_mixdown(mech_freq)
            print(f"setting ch06  to {current_V:6g} mV")
            time.sleep(5)

    def therm_driven_spectra(self,f_mech,reps_nodrive=10):
        Temp=0.13
        nr_bursts = 7
        reps_nodrive = 10                            
        demod_ch = 3
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

    def measure_singledot_config(self,temp_meas_counts=3,
                                 therm_reps=50,                   ##########
                                 driven_spectra=True,               
                                 driven_reps=5,                   #########                                  
                                 thermal_softening=True,
                                 softening_reps=9,                ##########              
                                 power_sweep=True):

        self.linesweep_parallel_LFsens_extended(costum_prefix='adjustment_linesweep',
                                                sitside="right",
                                                check_around_current_V=True,
                                                check_V_range=[-0.03,0.03],    ##########
                                                check_pt_pitch=2e-3,          ###########
                                                set_best_sitpos=True,
                                                find_startpos=True,
                                                main_gate=qdac.ch02.dc_constant_V)

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

    pos_list=[3.0,2.9,2.8,2.7,2.6,2.5,2.4,2.3,2.18,2.08,1.95,1.85,1.74,1.64,1.51,1.41,1.3,1.2,1.06,0.97] 
                                 #################
    #pos_list2=[2.8,2.7,2.6,2.5,2.4,2.3,2.18,2.08,1.95,1.85,1.74,1.64,1.51,1.41,1.3,1.2,1.06,0.97]
    def go_through_gate_pos(self,pos_list=pos_list,gate=qdac.ch02):
        for pos in pos_list:
            print(f"ramping to next step at {pos}")
            qdac.ramp_multi_ch_slowly([gate],[pos],step_size=4e-2,ramp_speed=4e-3)
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