# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 16:19:10 2024

@author: rtormo
"""
r'''Some of the information in this script is taken from: 
    W:\Electromechanics\Projects\DQD-mech-qubit\Devices\Chiral_August2024\Data\FFT_AM_signal
    The PhD thesis of William Nielsen (NBI), 2016
    '''
from tqdm import tqdm
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, hbar, e, h
from scipy.constants import Boltzmann as kB
from statsmodels.graphics.tsaplots import plot_acf
from scipy.signal import savgol_filter
from statsmodels.tsa.stattools import acf
from scipy.optimize import curve_fit
from dataprocessing.extract_fkts import *

# fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
# fig.set_size_inches(10,2)

# for fname in os.listdir('.'):
#     # if fnmatch.fnmatch(fname, 'dev2102_demods_7_sample_r_avg_00000.csv'):
#     if fnmatch.fnmatch(fname, 'dev2102_demods_4_sample_y_avg_00000.csv'):        
#         if fnmatch.fnmatch(fname, '*header*'):
#             print(f'Found file {fname}: for header')
#             pass
#         else:
#             print(f'Found file {fname}: for data')
            # data = np.genfromtxt(fname,skip_header=1,delimiter=';',dtype=float)             
            # demod = fname.split('_')[2]
            # variable_name = fname.split('_')[4] 
            
            # chunk = data[:,0]  
            # chunk_idx = np.array([0])
            # for val in np.where(chunk[:-1] != chunk[1:])[0]:
            #     chunk_idx = np.append(chunk_idx,val+1) 
            
            # chunk_idx = np.append(chunk_idx,len(chunk))
            
            # sampling_rate = 13.73e3
            # # sampling_rate = 1.717e3
            
            # var = data[:,2]
            # time = np.linspace(0,
            #                    len(var)*1/sampling_rate,
            #                    num=len(var))

            # for idx, idx_num in enumerate(chunk_idx):
                
            #     if idx <  len(chunk_idx)-1 :

            #         time_data = time[idx_num:chunk_idx[idx+1]]-time[idx_num]
            #         var_data = var[idx_num:chunk_idx[idx+1]]

#x_data = np.load('.\\test_with_stefan_data\\autocorr_t_5s.npy')
#r_data = np.load('.\\test_with_stefan_data\\autocorr_vr_5s.npy')
#phase_data = np.load('.\\test_with_stefan_data\\autocorr_Phase_5s.npy')
#t_data = np.load('.\\test_with_stefan_data\\autocorr_x_5s.npy')
run_id=134

decay_ms_filter_guess=0.01
decay_ms_mech_guess=6

t_data,x_data = extract_1d(run_id, data_1d_name = "x", setpoint_name = 'time_param',  plot = True)
_,y_data = extract_1d(run_id, data_1d_name = "v_r", setpoint_name = 'time_param',  plot = True)
_,phase_data = extract_1d(run_id, data_1d_name = "Phase", setpoint_name = 'time_param',  plot = True)

# var_data = x_data
# var_data = r_data
var_data = phase_data
time_data = t_data - t_data[0]

frac = 1
time_data = time_data[0:len(time_data)//frac]
var_data = var_data[0:len(var_data)//frac]


# =========================================================
# Raw time trace
# =========================================================

fig, ax = plt.subplots()
ax.plot(time_data,var_data)
ax.set_xlabel('time (s)')
fig.tight_layout()

#%% =========================================================
# DFT (Discrete Fourier Transform)   
# =========================================================
N = len(var_data)
duration = max(time_data)
sample_rate  = N/duration

from scipy.fft import fft, fftfreq
xf = fftfreq(N, 1/sample_rate)[0:len(time_data)//2]
yf = np.abs(fft(var_data))[0:len(time_data)//2]
yf_filtered = savgol_filter(yf, window_length=50, polyorder=2)

plt.figure()
plt.plot(xf,yf, label='Fourier transformed of data')
plt.plot(xf,yf_filtered, label='Filtered Fourier transformed of data')
plt.yscale('log')
plt.ylabel(r'FFT')
plt.xlabel(r'frequency (Hz)')
plt.legend()

#%% =========================================================
# Power dpesctral density                         
# =========================================================

from scipy import signal
f, Pxx_den = signal.periodogram(var_data, sample_rate, window='hamming')
plt.figure()
plt.title(f"Number of points={N}"
              f"\nTime trace={duration}s"
              f"\nSampling rate={round(sample_rate,2)}Hz")
plt.plot(f,Pxx_den)
plt.ylabel(r'$variable^2$ / Hz')
plt.xlabel(r'frequency (Hz)')
plt.yscale('log')
plt.xscale('log')

#%% =========================================================
# # Auto-correlation measurement from chatgpt v1
# # =========================================================
def autocorr_per_lag(x):
    """
    Computes the autocorrelation with per-lag normalization.
    """
    n = len(x)
    result = np.zeros(n)
    for k in tqdm(range(n)):
        # Overlapping portions for lag k
        x1 = x[:n - k]   # Original part
        x2 = x[k:]       # Lagged part
        # Calculate correlation for this lag
        result[k] = np.sum(x1 * x2) / ((n - k) * np.std(x1) * np.std(x2))
    return result

# Compute the autocorrelation
auto_corr = autocorr_per_lag(var_data)

# Plot
plt.figure()
plt.plot(time_data, auto_corr)
plt.xlabel('time (s)')
plt.ylabel('Normalized Autocorrelation')
plt.title('Autocorrelation Function (Per-Lag Normalization)')
# plt.xscale('log')
plt.tight_layout()
plt.show()

#%% =========================================================
# Auto-correlation measurement from chatgpt v2
# =========================================================
plt.figure()
def autocorr(x):
    result = np.correlate(x, x, mode='full')  # Full correlation
    return result[result.size // 2:]  # Return the second half

auto_corr = autocorr(var_data)

plt.plot(time_data, auto_corr)
plt.xlabel('time (s)')
# plt.xscale('log')
plt.ylabel('Autocorrelation')
plt.title('Autocorrelation Function of Variable Array (np.correlate)')
plt.tight_layout()

#%%
# =========================================================
# Auto-correlation measurement from Wei Yang
# =========================================================
num = len(time_data)
auto = np.zeros(num)

# version 2 Wei
# for n in tqdm(range(num)):
#     for m in range(num - n):
#         P = var_data[m] * var_data[m+n]
#         auto[n] += P
#     auto[n] /= (num - n)  # Normalize the product

# version 1 Wei
for n in tqdm(range(num)):
    # Fix: Use the correct slicing to ensure arrays are the same size
    P = var_data[0:num-n] * var_data[n:num]  # This ensures both arrays have the same length
    auto[n] = np.sum(P) / (num-n)  # Normalize the product
amplitude = np.abs(auto)

#%%
# Time in milliseconds
t_ms = time_data * 1e3
index = np.where((t_ms > -0.1) & (t_ms < 12))[0]

xData = t_ms[index]
yData = auto_corr[index]

# Initial guesses for curve fitting
begin_fit = index[0]
end_fit = index[0] + len(index)
a0 = (auto_corr[begin_fit] - auto_corr[end_fit]) 
b0 = 1 / decay_ms_mech_guess  # Initial guess for decay rate
c0 = auto_corr[end_fit] 
d0 = a0
e0 = 1 / decay_ms_filter_guess

# Define the fitting function
def fitting_func(x, a, b, c, d, e):
    return a * np.exp(-b * x) + c + d * np.exp(-e * x)

# Perform curve fitting
popt, _ = curve_fit(fitting_func, xData, yData, p0=[a0, b0, c0, d0, e0])
a, b, c, d, e = popt

decay_ms1 = 1 / b
decay_ms2 = 1 / e

Amp_fit = fitting_func(t_ms[index], a, b, c, d, e)

# Plot results
plt.figure()
plt.plot(t_ms[index], yData, label='Data')
plt.plot(t_ms[index], Amp_fit, label=r'Fit $\tau_1 =$'+ f'{decay_ms1:0.2f} ms\n'+r'$\tau_2 =$'+ f'{decay_ms2:0.2f} ms')
plt.xlabel('Time (ms)', fontsize=14)
plt.ylabel('Auto Correlation', fontsize=14)
plt.legend()
plt.grid(True)
plt.show()
plt.tight_layout()

#%% =========================================================
# DFT (Discrete Fourier Transform)   
# =========================================================
N = len(amplitude)
duration = max(time_data)
sample_rate  = N/duration

from scipy.fft import fft, fftfreq
xf = fftfreq(N, 1/sample_rate)[0:len(time_data)]
yf = np.abs(fft(amplitude))[0:len(time_data)]
yf_filtered = savgol_filter(yf, window_length=50, polyorder=2)

plt.figure()
plt.plot(xf[0:len(xf)//2],yf[0:len(xf)//2], label='Fourier transformed of autocorr')
# plt.plot(xf[len(xf)//2::],yf[len(xf)//2::],)
# plt.plot(xf,yf_filtered, label='Filtered Fourier transformed of autocorr')
plt.yscale('log')
plt.ylabel(r'FFT')
plt.xlabel(r'frequency (Hz)')
plt.legend()
