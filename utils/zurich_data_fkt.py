import os
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
from zhinst.toolkit import Session

# Function to average every n points
import numpy as np

def average_every_n_points(data, n=100):
    """
    Averages every n points in the data array.

    Parameters:
        data (numpy array): The data to be averaged.
        n (int): The number of points to average over.

    Returns:
        averaged_data (numpy array): The averaged data values.
    """
    # Trim data to a multiple of n
    length = len(data) // n * n
    data = data[:length]
    
    # Reshape and compute mean for each bin
    reshaped_data = np.reshape(data, (-1, n))
    averaged_data = np.mean(reshaped_data, axis=1)
    
    return averaged_data


def take_spectrum(demod_ch=3,BURST_DURATION = 4.772,SAMPLING_RATE = 13730,nr_burst=10, avg_num=100):
# Initialize session and device
    try:
        session = Session("localhost")
        device = session.connect_device("DEV20039")  # Replace with your actual device ID
        device.demods[0].enable(True)
       # print("Connected to the device successfully.")
    except Exception as e:
        print(f"Error connecting to device: {e}")
        sys.exit()  # Exit on connection error

    # Subscribe to the sample nodes (demodulator X and Y signals)
    sample_node = device.demods[demod_ch].sample.xiy.fft.abs.avg
    filter_node = device.demods[demod_ch].sample.xiy.fft.abs.filter

    # DAQ parameters
    TOTAL_DURATION = BURST_DURATION * nr_burst # Total duration in seconds for data collection
    #SAMPLING_RATE = 13730  # Samples per second
    #BURST_DURATION = 4.772  # Duration of each burst in seconds

    # Calculate number of columns and bursts for the DAQ module
    num_cols = int(np.ceil(SAMPLING_RATE * BURST_DURATION))
    num_bursts = int(np.ceil(TOTAL_DURATION / BURST_DURATION))

    # Configure the DAQ module
    daq_module = session.modules.daq
    daq_module.device(device)
    daq_module.type(0)  # Continuous acquisition
    daq_module.grid.mode(2)
    daq_module.count(num_bursts)  # Set number of bursts to collect
    daq_module.duration(BURST_DURATION)  # Set duration of each burst
    daq_module.grid.cols(num_cols)  # Set number of columns

    # Subscribe to sample node

    daq_module.subscribe(sample_node)
    daq_module.subscribe(filter_node)
    # Retrieve clock base for timing
    clockbase = device.clockbase()

    #   Start data acquisition
    daq_module.execute()
    time.sleep(2)  # Allow some time for the system to warm up

    # Prepare to collect results
    #results = {node: [] for node in sample_nodes}
    averaged_data_per_burst = []
    full_data=[]
    filter_data=[]

# Collect data over the specified duration
    for burst_idx in range(num_bursts):
        try:
            # Read data from the DAQ
            daq_data = daq_module.read(raw=False, clk_rate=clockbase)
            #print(f"Reading burst {burst_idx + 1}...")

            #for node in sample_nodes:
            if sample_node in daq_data.keys():
                for sig_burst in daq_data[sample_node]:
                    #results[sample_node].append(sig_burst)  # Store the signal burst
                        
                    value = sig_burst.value[0, :]  # Get the value of the signal
                    if value.size > 0:  # Check if value is not empty
                        # Apply averaging every 100 points
                        full_data.append(value)
                        averaged_value, averaged_frequency = average_every_n_points(data=value, n=avg_num)
                        averaged_data_per_burst.append(averaged_value)  # Store the averaged burst
                        #print(f"Burst {burst_idx + 1}: Appended averaged value with shape {averaged_value.shape}")
                   # else:
                      #  print(f"Burst {burst_idx + 1}: No data in this value.")
            else:
                print(f"Burst {burst_idx + 1}: No data available for node {sample_node}")

            if filter_node in daq_data.keys():
                #print("filter node available")
                for filter_burst in daq_data[filter_node]:
                    #results[sample_node].append(sig_burst)  # Store the signal burst
                        
                    filter_value = filter_burst.value[0, :]  # Get the value of the signal
                    if filter_value.size > 0:  # Check if value is not empty
                        #print("reading filter")
                        filter_data.append(value)
                        #averaged_value = average_every_n_points(value, n=100)
                        #averaged_data_per_burst.append(averaged_value)  # Store the averaged burst
#                        print(f"adding filter data")
                   # else:
                   #     print(f"Burst {burst_idx + 1}: No data in this filter.")

            time.sleep(BURST_DURATION)  # Wait for the next burst

        except Exception as e:
            print(f"Error reading data for burst {burst_idx + 1}: {e}")

    if averaged_data_per_burst:
        #print("äverage data")
        averaged_data = np.mean(averaged_data_per_burst, axis=0)

    num_samples = np.shape(full_data[-1])[0]  # Number of FFT bins (samples per burst)
   # print(num_samples)
    
        # Generate both positive and negative frequencies
    freqs = np.fft.fftfreq(num_samples, 1 / SAMPLING_RATE)
    freqs_sorted = np.sort(freqs)
    compressed_freqs = average_every_n_points(freqs_sorted, n=avg_num)



    return full_data, averaged_data_per_burst, averaged_data, freqs_sorted, compressed_freqs, filter_data


def demod_xy_timetrace(sample_nodes,daq_module,device,demod_ch=3,BURST_DURATION = 1,SAMPLING_RATE = 13730,nr_burst=10):
# Initialize session and device
#    try:
#        session = Session("localhost")
#        device = session.connect_device("DEV20039")  # Replace with your actual device ID
#        device.demods[0].enable(True)
#        print("Connected to the device successfully.")
#    except Exception as e:
#        print(f"Error connecting to device: {e}")
#        sys.exit()  # Exit on connection error

    # Subscribe to the sample nodes (demodulator X and Y signals)
    #sample_nodes = [
    #device.demods[demod_ch].sample.x,
    #device.demods[demod_ch].sample.y
    #]

    #sample_node = device.demods[demod_ch].sample.xiy.fft.abs.avg
    #filter_node = device.demods[demod_ch].sample.xiy.fft.abs.filter

    # DAQ parameters
    TOTAL_DURATION = BURST_DURATION * nr_burst # Total duration in seconds for data collection
    #SAMPLING_RATE = 13730  # Samples per second
    #BURST_DURATION = 4.772  # Duration of each burst in seconds

    # Calculate number of columns and bursts for the DAQ module
    num_cols = int(np.ceil(SAMPLING_RATE * BURST_DURATION))
    num_bursts = int(np.ceil(TOTAL_DURATION / BURST_DURATION))

    # Configure the DAQ module
#    daq_module = session.modules.daq
#    daq_module.device(device)
#    daq_module.type(0)  # Continuous acquisition
#    daq_module.grid.mode(2)
    daq_module.count(num_bursts)  # Set number of bursts to collect
    daq_module.duration(BURST_DURATION)  # Set duration of each burst
    daq_module.grid.cols(num_cols)  # Set number of columns

    # Subscribe to sample node
    #for node in sample_nodes:
    #    daq_module.subscribe(node)
    #daq_module.subscribe(sample_nodes[1])
    # Retrieve clock base for timing
    clockbase = device.clockbase()

    #   Start data acquisition
    #daq_module.execute()
    time.sleep(2)  # Allow some time for the system to warm up

    # Prepare to collect results
    #results = {node: [] for node in sample_nodes}
    #averaged_data_per_burst = []
    full_x_data=[]
    full_y_data=[]
    full_t_data=[]
    time_offset=0
# Collect data over the specified duration
    for burst_idx in range(num_bursts):
        try:
            # Read data from the DAQ
            daq_data = daq_module.read(raw=False, clk_rate=clockbase)
            #print(f"Reading burst {burst_idx + 1}...")
            
            for node in sample_nodes:
                if node in daq_data.keys():
                    print("node in keys")
                    for sig_burst in daq_data[node]:
                    #results[sample_node].append(sig_burst)  # Store the signal burst
                        
                        value = sig_burst.value[0, :]  # Get the value of the signal
                        #t0_burst = sig_burst.header['createdtimestamp'][0] / clockbase
                        t = (sig_burst.time)+time_offset   
                        #print(f"t0={t0_burst}")
                        full_t_data.append(t)
                        if value.size > 0:  # Check if value is not empty
                        # Apply averaging every 100 points
                            if node==device.demods[demod_ch].sample.x:
                                full_x_data.append(value)
                                print("appending x")
                            if node==device.demods[demod_ch].sample.y:
                                full_y_data.append(value)
                                print("appending y")
                        else:
                            print(f"Burst {burst_idx + 1}: No data in this value.")
                else:
                    print(f"Burst {burst_idx + 1}: No data available for node {node}")

           

        except Exception as e:
                print(f"Error reading data for burst {burst_idx + 1}: {e}")
        time.sleep(BURST_DURATION)
        time_offset+=BURST_DURATION
    
    


    return full_t_data,full_x_data,full_y_data
