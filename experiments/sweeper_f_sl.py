#from instruments import *
 
#measure nonlinearity 310725
import time
import copy
from instruments import *
import qcodes as qc
from tqdm import tqdm
#Vg=0#init
#from experiments.zurich_experiments.spectrum_0825 import *
exp_name="sweeper"
device_name=exp.device_name
 
start_fsl=0
stop_fsl=100e6
step_num=5000
tc=30e-3
vsdac=10e-6#random filler
 
slf_sweep = zurich.freq0.sweep(start=start_fsl, stop=stop_fsl, num=step_num)
 
 
 
experiment = new_experiment(name=exp_name, sample_name=device_name)
meas = Measurement(exp=experiment)
meas.register_parameter(slf_sweep.parameter)
#meas.register_custom_parameter('G', unit='S', setpoints=[vgdc_sweep.parameter])
meas.register_custom_parameter('V_r', unit='V', setpoints=[slf_sweep.parameter])
meas.register_custom_parameter('Phase', unit='rad', setpoints=[slf_sweep.parameter])
#meas.register_custom_parameter('R', unit='Ohm', setpoints=[vgdc_sweep.parameter])
#meas.register_custom_parameter('I', unit='A', setpoints=[vgdc_sweep.parameter])
 
with meas.run() as datasaver:
                #varnames = [str(name) for name in [tc, vsd_dB, amp_lvl, x_avg, y_avg]]
                #save_metadata_var(datasaver.dataset, varnames, [tc, vsd_dB, amp_lvl, x_avg, y_avg])
                qdac.add_dc_voltages_to_metadata(datasaver=datasaver)
                zurich.save_config_to_metadata(datasaver=datasaver)
               
                for slf_value in tqdm(slf_sweep, desc=f'sl frequency Sweep{zurich.freq1()}'):
                    zurich.freq0(slf_value)
                    time.sleep(1.1 * tc)
 
                    # Some measurement
                    #_ = measured_parameter()
                    theta_calc, v_r_calc, I, G = zurich.phase_voltage_current_conductance_compensate(vsdac)
                    
 
                    datasaver.add_result(
                        #('R', R),
                        #('G', G),
                        ('V_r', v_r_calc),
                        ('Phase', theta_calc),
                        #('I', I),
                        (slf_sweep.parameter, slf_value)
                    )
 
                    
     
            
                
       
 
   
 