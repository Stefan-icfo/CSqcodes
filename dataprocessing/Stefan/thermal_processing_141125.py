import matplotlib.pyplot as plt
from dataprocessing.extract_fkts_save import *
from utils.CS_utils import *
import qcodes as qc
from database import *


def voltage_to_psd(v_rms, rbw, impedance=50):
  
    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd


def analyze_fit_quality_g(x, y, popt):
    """
    This tells you how good your fit actually is
    """
    y_fit = gaussian_fkt(x, *popt)
    residuals = y - y_fit
    
   
    

    return residuals

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

#qc.config["core"]["db_location"] ='.\\Data\\Raw_data\\CD12_B5_F4v18_171025.db'#4e untensioned
#print("Opening DB:", qc.config["core"]["db_location"])

#run_id=540
#f_bounds=[142.0645e6,142.0685e6]
#background_level=-5e-9


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



qc.config["core"]["db_location"] ='.\\Data\\Raw_data\\CD12_B5_F4v45_05_12_25.db'#1e
print("Opening DB:", qc.config["core"]["db_location"])

run_id=118
f_bounds=[137.464e6,137.474e6]
background_level=620e-9

run_id=440
f_bounds=[137.464e6,137.474e6]
background_level=110e-9


run_id=496#60
f_bounds=[137.464e6,137.474e6]
background_level=110e-9


#run_id=499#80
f_bounds=[137.464e6,137.474e6]
background_level=110e-9

run_id=517#200
f_bounds=[137.464e6,137.474e6]
background_level=110e-9

#run_id=575#205mK
#run_id=577#205mK
run_id=579#205mK
#f_bounds=[137.464e6,137.474e6]
#background_level=50e-9


#run_id=588#225mK
#run_id=586#225mK
#run_id=584#225mK
#f_bounds=[137.464e6,137.474e6]
#background_level=50e-9


#run_id=597#250mK
#run_id=599#250mK
#run_id=584#250mK
#f_bounds=[137.464e6,137.474e6]
#background_level=50e-9

#run_id=609#275mK
#run_id=611#275mK
#run_id=613#275mK
#f_bounds=[137.464e6,137.474e6]
#background_level=50e-9


#run_id=621#300mK
#run_id=623#300mK
#run_id=#300mK
#f_bounds=[137.465e6,137.475e6]
#background_level=50e-9

#run_id=633#350mK
#run_id=635#350mK
#run_id=637#350mK
#f_bounds=[137.465e6,137.475e6]
#background_level=50e-9

#run_id=645#400mK
#run_id=647#400mK
#run_id=649#400mK
#f_bounds=[137.466e6,137.476e6]
#background_level=50e-9


#run_id=657#450mK
#run_id=659#450mK
#run_id=661#450mK
#f_bounds=[137.466e6,137.476e6]
#background_level=50e-9


#run_id=669#500mK
#run_id=671#500mK
#run_id=673#500mK
#f_bounds=[137.4675e6,137.4775e6]
#background_level=50e-9

run_id=683#600mK
f_bounds=[137.469e6,137.479e6]
background_level=50e-9



#run_id=695#700mK
#f_bounds=[137.472e6,137.482e6]
#background_level=50e-9

#run_id=715#8-900mK
#f_bounds=[137.474e6,137.484e6]
#background_level=50e-9


#run_id=730#1K
#f_bounds=[137.476e6,137.486e6]
#background_level=50e-9



#qc.config["core"]["db_location"] ='.\\Data\\Raw_data\\CD12_B5_F4v46_07_12_25.db'#1e
#print("Opening DB:", qc.config["core"]["db_location"])
#run_id=12#600mK
#f_bounds=[137.466e6,137.476e6]
#background_level=50e-9

metadata=get_metadata(meas_id=run_id-1,return_data=True)
rbw=metadata['rbw']
Isens=11e-12#metadata['I_sens_sit']

freq_axis,spectrum=extract_1d(run_id, data_1d_name = "V_fft_avg_avg", setpoint_name = 'freq_param',  plot = False,return_exp_name=False)

mask = (freq_axis >= f_bounds[0]) & (freq_axis <= f_bounds[1])

substracted_Vfft=spectrum[mask]-background_level

plt.plot(freq_axis[mask],substracted_Vfft)
plt.show()

subtracted_spectrum=voltage_to_psd(substracted_Vfft, rbw)

sum_substracted_spectrum=np.sum(subtracted_spectrum)

print(f"spectral_area={sum_substracted_spectrum}")
print(f"sensitivity={Isens}")



############fitting##################
frequency=(f_bounds[0]+f_bounds[1])/2
Gamma_guess=2e3

initial_guess=[frequency,Gamma_guess,max(centered_moving_average(spectrum,n=5)),0]

popt, pcov = scp.optimize.curve_fit(gaussian_fkt, freq_axis[mask], centered_moving_average(subtracted_spectrum,n=5), p0=initial_guess)
gaussian,g_area=gaussian_fkt_w_area(freq_axis[mask],popt[0],popt[1],popt[2],popt[3])


manual_psd_offset=6.1e-19

area_corrected=sum_substracted_spectrum-len(subtracted_spectrum)*manual_psd_offset

plt.plot(freq_axis[mask],gaussian-manual_psd_offset)
plt.plot(freq_axis[mask],centered_moving_average(subtracted_spectrum,n=5)-manual_psd_offset)
plt.plot([popt[0]-3*popt[1],popt[0]+3*popt[1]],[0,0],"r*")
plt.show()


reduced_freq_axis=freq_axis[mask]
fit_quality_mask=(reduced_freq_axis >= popt[0]-3*popt[1]) & (reduced_freq_axis <= popt[0]+3*popt[1])



#residuals_g = analyze_fit_quality_g(reduced_freq_axis[fit_quality_mask], subtracted_spectrum[fit_quality_mask], popt)
residuals_g = analyze_fit_quality_g(reduced_freq_axis, subtracted_spectrum, popt)
area_error_g=sum(abs(residuals_g))

plt.plot(reduced_freq_axis, subtracted_spectrum)
plt.plot(reduced_freq_axis,residuals_g)
#plt.plot([popt[0]-3*popt[1],popt[0]+3*popt[1]],[0,0],"r*")
plt.show()


print(f"area_error_g={area_error_g}")