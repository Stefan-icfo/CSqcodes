import matplotlib.pyplot as plt
from dataprocessing.extract_fkts_save import *
from utils.CS_utils import *
import qcodes as qc
from database import *


def voltage_to_psd(v_rms, rbw, impedance=50):
  
    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd


#qc.config["core"]["db_location"] ='.\\Data\\Raw_data\\CD12_B5_F4v27_06_11_25.db'#1e
#print("Opening DB:", qc.config["core"]["db_location"])

#run_id=195
#f_bounds=[142.8275e6,142.8315e6]
#background_level=-20e-9

qc.config["core"]["db_location"] ='.\\Data\\Raw_data\\CD12_B5_F4v30_14_11_25.db'#1e
print("Opening DB:", qc.config["core"]["db_location"])

run_id=979
f_bounds=[142.828e6,142.832e6]
background_level=620e-9

#qc.config["core"]["db_location"] ='.\\Data\\Raw_data\\CD12_B5_F4v28_10_11_25.db'#234 e
#print("Opening DB:", qc.config["core"]["db_location"])

###Values for 2nd electron; tensioned config

#run_id=201
#f_bounds=[142.922e6,142.926e6]
#background_level=28e-9

###Values for 3rd electron; tensioned config
#run_id=250
#f_bounds=[142.616e6,142.620e6]
#background_level=25e-9


###Values for 4th electron; tensioned config
#run_id=697
#f_bounds=[141.103e6,141.107e6]
#background_level=40e-9

metadata=get_metadata(meas_id=run_id-1,return_data=True)
rbw=metadata['rbw']
Isens=11e-12#metadata['I_sens_sit']

freq_axis,spectrum=extract_1d(run_id, data_1d_name = "V_fft_avg_avg", setpoint_name = 'freq_param',  plot = True,return_exp_name=False)

mask = (freq_axis >= f_bounds[0]) & (freq_axis <= f_bounds[1])

substracted_Vfft=spectrum[mask]-background_level

plt.plot(freq_axis[mask],substracted_Vfft)
plt.show()

subtracted_spectrum=voltage_to_psd(substracted_Vfft, rbw)

sum_substracted_spectrum=np.sum(subtracted_spectrum)

print(f"spectral_area={sum_substracted_spectrum}")
print(f"sensitivity={Isens}")
