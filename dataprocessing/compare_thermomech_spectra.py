import numpy as np

from utils.CS_utils import lorentzian_fkt

from dataprocessing.extract_fkts import *
from dataprocessing.processing_fkt import *
from dataprocessing.specialized_extract_fkts import *
import matplotlib.pyplot as plt

freq_20mK,rep_axis_20mK,real_time_axis_20mK,freq_spectrum_real_20mK,time_avg_psd_20mK=extract_2d_spectra_w_repeatingtimeaxis(97, plot = False)

freq_axis_30mK, reps_axis_30mK, data_2d_30mK=extract_2d(run_id=2940,data_2d_name='psd', setpoints1_name='freq_param',setpoints2_name='time_param',plot = False)

freq_axis_50mK, reps_axis_50mK, data_2d_50mK=extract_2d(run_id=2760,data_2d_name='psd', setpoints1_name='freq_param',setpoints2_name='time_param',plot = False)

freq_axis_100mK, reps_axis_100mK, data_2d_100mK=extract_2d(run_id=2802,data_2d_name='psd', setpoints1_name='freq_param',setpoints2_name='time_param',plot = False)

freq_axis_200mK, reps_axis_200mK, dataV_2d_200mK=extract_2d(run_id=2538,data_2d_name='Voltage_fft_avg', setpoints1_name='freq_param',setpoints2_name='time_param',plot = False)
data_2d_200mK=voltage_to_psd(dataV_2d_200mK, 0.209, impedance=50)

freq_axis_400mK, reps_axis_400mK, data_2d_400mK=extract_2d(run_id=2864,data_2d_name='psd', setpoints1_name='freq_param',setpoints2_name='time_param',plot = False)

freq_axis_700mK, reps_axis_700mK, data_2d_700mK=extract_2d(run_id=2912,data_2d_name='psd', setpoints1_name='freq_param',setpoints2_name='time_param',plot = False)

#List_of_lolos30mK,data_interp_list30mK,avg_shifted_data30mK,avg_freq_fit30mK,avg_hgamma_fit30mK,avg_amp_fit30mK,avg_backg_fit30mK,area30mK=fit_lorentzians_and_average(data_2d_30mK,freq_axis_30mK,reps_axis_30mK)
#actual_data_interpolation_30mK,freq_fit_30mK, hgamma_fit_30mK, amp_fit_30mK,backg_fit_30mK,area_30mK=average_and_then_fitlorentizian(data_2d_30mK,freq_axis_30mK,reps_axis_30mK)

List_of_lolos50mK, data_interp_list50mK, avg_shifted_data50mK, avg_freq_fit50mK, avg_hgamma_fit50mK, avg_amp_fit50mK, avg_backg_fit50mK, area50mK = fit_lorentzians_and_average(data_2d_50mK, freq_axis_50mK, reps_axis_50mK)

#List_of_lolos100mK, data_interp_list100mK, avg_shifted_data100mK, avg_freq_fit100mK, avg_hgamma_fit100mK, avg_amp_fit100mK, avg_backg_fit100mK, area100mK = fit_lorentzians_and_average(data_2d_100mK, freq_axis_100mK, reps_axis_100mK)
#actual_data_interpolation_100mK,freq_fit_100mK, hgamma_fit_100mK, amp_fit_100mK,backg_fit_100mK,area_100mK=average_and_then_fitlorentizian(data_2d_100mK,freq_axis_100mK,reps_axis_100mK)


List_of_lolos200mK, data_interp_list200mK, avg_shifted_data200mK, avg_freq_fit200mK, avg_hgamma_fit200mK, avg_amp_fit200mK, avg_backg_fit200mK, area200mK = fit_lorentzians_and_average(data_2d_200mK, freq_axis_200mK, reps_axis_200mK)

List_of_lolos400mK, data_interp_list400mK, avg_shifted_data400mK, avg_freq_fit400mK, avg_hgamma_fit400mK, avg_amp_fit400mK, avg_backg_fit400mK, area400mK = fit_lorentzians_and_average(data_2d_400mK, freq_axis_400mK, reps_axis_400mK)

#List_of_lolos700mK, data_interp_list700mK, avg_shifted_data700mK, avg_freq_fit700mK, avg_hgamma_fit700mK, avg_amp_fit700mK, avg_backg_fit700mK, area700mK = fit_lorentzians_and_average(data_2d_700mK, freq_axis_700mK, reps_axis_700mK)
actual_data_interpolation_700mK,freq_fit_700mK, hgamma_fit_700mK, amp_fit_700mK,backg_fit_700mK,area_700mK=average_and_then_fitlorentizian(data_2d_700mK,freq_axis_700mK,reps_axis_700mK)


crosscap_values=[2.25,3,3.75,2.9,4,4.5]#from 120MHz
crosscap_values=[28,52,35,120,80,35,380]#from 1uV excitation without considering q
#crosscap_values=[1,1,1,1,1,1,1]
# Plot data for 30mK
plt.plot(freq_axis_30mK+3300, (np.mean(data_2d_30mK,axis=0)-0.2e-14)/crosscap_values[0], label='30mK')
# Plot data for 50mK
plt.plot(freq_axis_50mK, avg_shifted_data50mK/crosscap_values[1], label='50mK')
# Plot data for 100mK
plt.plot(freq_axis_100mK+2700, (np.mean(data_2d_100mK,axis=0)-0.3e-14)/crosscap_values[2], label='100mK')
# Plot data for 200mK
plt.plot(freq_axis_200mK, avg_shifted_data200mK/crosscap_values[3], label='200mK')
# Plot data for 400mK
plt.plot(freq_axis_400mK, avg_shifted_data400mK/crosscap_values[4], label='400mK')
# Plot data for 700mK
plt.plot(freq_axis_700mK-4300, (np.mean(data_2d_700mK,axis=0)-0.9e-14)/crosscap_values[5], label='700mK')
# Plot data for 20mK
plt.plot(freq_20mK, (time_avg_psd_20mK-0.2e-14)/crosscap_values[6], label='20mK')
# Add labels and title
plt.xlim(-5000, 5000)
plt.xlabel('Frequency (Hz)')  # Adjust the unit if necessary
plt.ylabel('Avg Shifted psd')
plt.title('Avg Shifted normalized psd vs Frequency for Different Temperatures')

# Add a legend
plt.legend()

# Show the plot
plt.show()
