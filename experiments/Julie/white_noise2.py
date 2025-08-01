import numpy as np
import pyvisa
import time

from instruments import qdac


import random




qdac.write('trace:define "WhiteNoise",40000') # Defines a trace of 40 points (1 μs)
qdac.ask("trac:cat?") # check that trace has been created


qdac.write('sour8:awg:define "WhiteNoise"') # assign the trace to channel AWG

traclen= 40000
values = [random.uniform(-0.5,0.5) for _ in range(traclen)]


qdac.write_floats \
("trace:data \"WhiteNoise\",",values)
qdac.write('sour8:awg:trig:sour imm') 

qdac.write('sour8:awg:init') 

#qdac.write('sour8:awg:abor')
#to close the channel and delete the trace 
#qdac.write("trace:remove:all") # Deletes all traces in trace memory
#qdac.ch08.dc_constant_V(0) #force set it to zero 

# then reset and to be safe set all slewrates.. 
# for some reason even force setting it to zero doesnt work to close the channel
"""
qdac.write("*RST")

qdac.set_all_slewrates()
"""