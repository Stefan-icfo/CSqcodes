import numpy as np
import os

from instruments import qdac, manual, station, k2400
from qcodes.dataset import Measurement, new_experiment
from utils.sample_name import sample_name
import drivers.k2400 as k2
import time
from tqdm import tqdm
import scipy as scp
import matplotlib.pyplot as pt

print("time taken to read current")
timezero=time.time()
k2400.curr()
print((time.time()-timezero)*1000)

#print("time taken to set voltage to 100uV and read current")
#timezero=time.time()
#k2400.volt(1e-4)
#k2400.curr()
#print((time.time()-timezero)*1000)
