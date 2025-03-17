import numpy as np
import pyvisa
import time
from scipy.signal import butter, filtfilt
from instruments import qdac

import random
channel = 8
# Define a SCPI command (note the space at the end if needed)
command = "trace:data CH8, "

traclen = 1000
# Generate white noise values between -1 and 1
values = [random.uniform(-0.1,0.1) for _ in range(traclen)]

# Write the white noise trace to channel 8 of QDAC 2.
#qdac.write_binary_values("trace:data \"Ch8 White Noise\",", values)
qdac.write_floats(cmd=command, values=values)

print(f"White noise (100 kHz bandwidth, 0.1 V RMS) applied to QDAC channel {channel} with a 0.01 s trace.")
