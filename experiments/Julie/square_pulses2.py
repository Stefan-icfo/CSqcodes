import numpy as np
import numpy as np
import pyvisa
import time

from instruments import qdac


import random

reps = 10
pts = 5000
voltage = 612e-6

# Define square wave values
values = np.tile(np.concatenate([np.full(pts, voltage), np.full(pts, -voltage)]), reps)

# Define trace
trace_name = 'squarepulses'
qdac.write(f"trace:define '{trace_name}',{len(values)}")

# Upload waveform data
qdac.write_floats(f"trace:data '{trace_name}',", values)

# Assign trace to AWG on channel 8
qdac.write('sour4:func:mode awg')
qdac.write(f"sour4:awg:define '{trace_name}'")

# Trigger and start
qdac.write('sour4:awg:trig:sour imm')
qdac.write('sour4:awg:stat on')
qdac.write('sour4:awg:init')
qdac.write('outp4 on')

#qdac.write('sour4:awg:abor')
#to close the channel and delete the trace 
#qdac.write("trace:remove:all") # Deletes all traces in trace memory
