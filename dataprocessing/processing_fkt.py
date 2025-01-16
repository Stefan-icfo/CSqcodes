import numpy as np
import matplotlib.pyplot as plt

def autocorrelation(time_series, time_array):
    """
    Computes the autocorrelation of a 1D numpy array.

    Parameters:
        time_series (numpy.ndarray): The input time series data.
        time_array (numpy.ndarray): The time values corresponding to the time series data.

    Returns:
        lags (numpy.ndarray): Array of lag values in terms of time.
        autocorr (numpy.ndarray): Autocorrelation values for each lag.
    """
    if len(time_series) != len(time_array):
        raise ValueError("time_series and time_array must have the same length.")

    n = len(time_series)
    mean = np.mean(time_series)
    variance = np.var(time_series)
    autocorr = np.correlate(time_series - mean, time_series - mean, mode='full') / (variance * n)
    
    # Take the second half of the autocorrelation (positive lags)
    autocorr = autocorr[n-1:]

    # Compute lag times
    time_diffs = time_array[1:] - time_array[:-1]
    avg_time_step = np.mean(time_diffs)
    lags = np.arange(len(autocorr)) * avg_time_step

    return lags, autocorr

from datetime import datetime


def autocorrelation_cld(data, timestamps, max_lag_seconds):
    """
    Calculate and plot autocorrelation for a 1D time series with timestamps
    
    Parameters:
    data: 1D numpy array containing the time series values
    timestamps: 1D numpy array containing the corresponding timestamps in seconds
    max_lag_seconds: maximum time lag to calculate (in seconds)
    """
    if len(data) != len(timestamps):
        raise ValueError("Data and timestamps must have the same length")
    
    # Sort data by timestamps (in case they're not ordered)
    sort_idx = np.argsort(timestamps)
    timestamps = timestamps[sort_idx]
    data = data[sort_idx]
    
    # Calculate mean and normalize data
    mean = np.mean(data)
    normalized_data = data - mean
    var = np.var(data)
    
    # Create time bins
    min_time = timestamps[0]
    max_time = min(timestamps[-1], min_time + max_lag_seconds)
    time_lags = np.arange(0, max_time - min_time + 1)  # 1-second bins
    
    # Calculate autocorrelation for each time lag
    acf = np.zeros(len(time_lags))
    n = len(data)
    
    for i, lag in enumerate(time_lags):
        # Find all pairs of points separated by approximately this time lag
        valid_pairs = 0
        correlation = 0
        
        for j in range(n):
            # Find points that are approximately 'lag' seconds after point j
            future_points = np.where(
                (timestamps - timestamps[j] >= lag - 0.5) & 
                (timestamps - timestamps[j] <= lag + 0.5)
            )[0]
            
            if len(future_points) > 0:
                correlation += np.sum(normalized_data[j] * normalized_data[future_points])
                valid_pairs += len(future_points)
        
        if valid_pairs > 0:
            acf[i] = correlation / (valid_pairs * var)
    
    # Plot the results
    plt.figure(figsize=(12, 6))
    plt.plot(time_lags, acf)
    plt.axhline(y=0, color='r', linestyle='--')
    plt.axhline(y=1.96/np.sqrt(n), color='r', linestyle='--')
    plt.axhline(y=-1.96/np.sqrt(n), color='r', linestyle='--')
    plt.xlabel('Time Lag (seconds)')
    plt.ylabel('Autocorrelation')
    plt.title('Autocorrelation Function')
    plt.grid(True)
    plt.show()
    
    return acf, time_lags