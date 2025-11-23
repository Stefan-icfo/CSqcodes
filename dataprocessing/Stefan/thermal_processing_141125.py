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