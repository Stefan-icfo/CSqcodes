# further data processing fkt
import scipy as scp
import numpy as np
from utils.CS_utils import lorentzian_fkt,lorentzian_fkt_w_area
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def voltage_to_psd(v_rms, rbw, impedance=50):
  
    # Calculate PSD using the formula
    psd = (v_rms ** 2) / (impedance * rbw)
    return psd


def fit_lorentzians_and_average(data_2d,freq_axis,reps_axis,initial_guess_gamma=200,initial_guess_offset=2e-15,lower_peak_bound=5e-15):
    List_of_lolos=[]
    data_interp_list=[]
    #sum_freq_fit=0
    #sum_hgamma_fit=0
    #sum_amp_fit=0
    #sum_backg_fit=0
    #bounds = ([-np.inf, 0, lower_peak_bound], [np.inf, upper_gamma_bound, np.inf])
    shifted_freq_axes = np.zeros([len(reps_axis), len(freq_axis)])  


    for i in range(len(reps_axis)):
        peak_guess_index = np.argmax(data_2d[i])
        initial_guess=[freq_axis[peak_guess_index], initial_guess_gamma, max(data_2d[i])/2,initial_guess_offset]
        #plt.plot(freq_axis,lorentzian_fkt(freq_axis,initial_guess[0],initial_guess[1],initial_guess[2],initial_guess[3]))
        #plt.plot(freq_axis,data_2d[i])
        #plt.show()
        #popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, freq_axis, data_2d[i],bounds=([-np.inf, 100, 1e-15, 1e-15], [np.inf, 500, 50e-15, 3e-15]))
        popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, freq_axis, data_2d[i])
        freq_fit, hgamma_fit, amp_fit,backg_fit=popt
        List_of_lolos.append(popt)
        #sum_freq_fit+=freq_fit
        #sum_hgamma_fit+=hgamma_fit
        #sum_amp_fit+=amp_fit
        #sum_backg_fit+=backg_fit
        shifted_freq_axes[i]=freq_axis-freq_fit
        interpolation_fun=interp1d(shifted_freq_axes[i], data_2d[i], kind='linear', fill_value="extrapolate")
        actual_data_interpolation = interpolation_fun(freq_axis)-backg_fit
        data_interp_list.append(actual_data_interpolation)

        
    Lolosarray=np.array(List_of_lolos)
    #Shifted_freq_axes_array=np.array(shifted_freq_axes)
            

    avg_freq_fit=np.mean(Lolosarray[:,0])
    avg_hgamma_fit=np.mean(Lolosarray[:,1])
    avg_amp_fit=np.mean(Lolosarray[:,2])
    avg_backg_fit=np.mean(Lolosarray[:,3])
    avg_shifted_data=sum(data_interp_list)/i
    temp_lorentzian, area = lorentzian_fkt_w_area(freq_axis,0,avg_hgamma_fit,avg_amp_fit,avg_backg_fit)


    return List_of_lolos,data_interp_list,avg_shifted_data,avg_freq_fit,avg_hgamma_fit,avg_amp_fit,avg_backg_fit,area

def average_and_then_fitlorentizian(data_2d,freq_axis,reps_axis,initial_guess_gamma=200,initial_guess_offset=2e-15,lower_peak_bound=5e-15):
    #List_of_lolos=[]
    #data_interp_list=[]
    #sum_freq_fit=0
    #sum_hgamma_fit=0
    #sum_amp_fit=0
    #sum_backg_fit=0
    #bounds = ([-np.inf, 0, lower_peak_bound], [np.inf, upper_gamma_bound, np.inf])
    shifted_freq_axis = np.zeros(len(freq_axis))  

    avg_data=np.mean(data_2d,axis=0)
        
    peak_guess_index = np.argmax(avg_data)
    initial_guess=[freq_axis[peak_guess_index], initial_guess_gamma, max(avg_data)/2,initial_guess_offset]
        #plt.plot(freq_axis,lorentzian_fkt(freq_axis,initial_guess[0],initial_guess[1],initial_guess[2],initial_guess[3]))
        #plt.plot(freq_axis,data_2d[i])
        #plt.show()
        #popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, freq_axis, data_2d[i],bounds=([-np.inf, 100, 1e-15, 1e-15], [np.inf, 500, 50e-15, 3e-15]))
    popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, freq_axis, avg_data)
    freq_fit, hgamma_fit, amp_fit,backg_fit=popt
        
        #sum_freq_fit+=freq_fit
        #sum_hgamma_fit+=hgamma_fit
        #sum_amp_fit+=amp_fit
        #sum_backg_fit+=backg_fit
    shifted_freq_axis=freq_axis-freq_fit
    interpolation_fun=interp1d(shifted_freq_axis, avg_data, kind='linear', fill_value="extrapolate")
    actual_data_interpolation = interpolation_fun(freq_axis)-backg_fit
       
        
    
    
    temp_lorentzian, area = lorentzian_fkt_w_area(freq_axis,0,hgamma_fit,amp_fit,backg_fit)


    return actual_data_interpolation,freq_fit, hgamma_fit, amp_fit,backg_fit,area
