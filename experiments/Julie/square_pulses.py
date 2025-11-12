import numpy as np
import pyvisa
import time

from instruments import qdac


import random








reps=10
pts=2000
voltage=0.1e-3
qdac.write(f'trace:define "squarepulses",{reps*pts*2}') # Defines a trace of 40 points (1 μs)
qdac.ask("trac:cat?") # check that trace has been created


qdac.write('sour8:awg:define "squarepulses"') # assign the trace to channel AWG


values=[]
for n in range(reps):
    values.extend(pts*[voltage])
    values.extend(pts*[-voltage])



qdac.write_floats 
("trace:data \"squarepulses\",",values)
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