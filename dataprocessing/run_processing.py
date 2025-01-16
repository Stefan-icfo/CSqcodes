import numpy as np
import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from dataprocessing.processing_fkt import *

time,demod_trace= extract_1d(2421, data_1d_name = "x",setpoint_name = 'time_param')

#data = np.random.randn(100)  # Replace with your time series data
#time = np.linspace(0, 10, 100)  # Replace with your time array (e.g., timestamps)

#lags, autocorr = autocorrelation(demod_trace, time)
acf, time_lags = autocorrelation_cld(demod_trace, time, max_lag_seconds=1)
"""
# Plot the autocorrelation
plt.figure(figsize=(8, 4))
plt.plot(lags, autocorr, marker='o')
plt.title('Autocorrelation 10 bursts x')
plt.xlabel('Lag (time units)')
plt.ylabel('Autocorrelation')
plt.xscale('log')  # Set x-axis to logarithmic scale
plt.grid()
plt.show()


# Compute the Fourier Transform of the autocorrelation
autocorr_fft = np.fft.fft(autocorr)
frequencies = np.fft.fftfreq(len(autocorr), d=lags[1] - lags[0])

# Plot the Fourier Transform
plt.figure(figsize=(8, 4))
plt.plot(frequencies, np.abs(autocorr_fft))
plt.title('Fourier Transform of Autocorrelation')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.grid()
plt.show()
"""