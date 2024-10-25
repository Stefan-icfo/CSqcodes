import numpy as np


from instruments import zurich
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


measured_parameter = zurich.demods.demods0.sample

timeout=100#s
tc=0.006
window_size=14000

#values=np.zeros(window_size)
values = np.zeros(window_size, dtype=complex)

#fill initial window
print("filling initial window")
for i in range(window_size):
    measured_value=measured_parameter()
    x = measured_value['x'][0]#SF: COMMENTED OUT 
    y = measured_value['y'][0]#SF: COMMENTED OUT
    xy_complex = complex(x,y)
    #xiy=x+y*j
    v = np.absolute(xy_complex)
    if i%1000==1:
        print(i)
    values[i]=xy_complex
    #time.sleep(tc)#stefan I love you. Roger****amazing
print("done filling initial window")
fig, ax = plt.subplots()
start_time=time.time()

while time.time()-start_time<timeout:
    signal_fft = np.fft.fft(values)
    magnitude_spectrum = np.abs(signal_fft) / window_size
    frequencies = np.fft.fftfreq(window_size, tc)

    ax.cla()  # Clear the previous plot
    #ax.plot(frequencies, magnitude_spectrum)  # Plot both positive and negative frequencies
    ax.scatter(frequencies, magnitude_spectrum , s=5)
    # Step 4: Set y-axis to log scale and limit y-axis bounds
    ax.set_yscale('log')
    ax.set_ylim(1e-9, 1e-3)  # Set bounds to 10^-9 to 10^-4 V/vHz
    #ax.set_xlim(0, 300)
    # Optional: You can set x-axis limits if needed
    ax.set_xlim(min(frequencies), max(frequencies))

    # Step 5: Force the plot to update and pause briefly
    plt.pause(0.01)  # Small pause for real-time update

    values=np.roll(values,-1)

    measured_value=measured_parameter()
    x = measured_value['x'][0]
    y = measured_value['y'][0]
    xy_complex = complex(x,y)
    v = np.absolute(xy_complex)
    values[-1]=xy_complex


#plt.plot(frequencies,magnitude_spectrum)
plt.show()

