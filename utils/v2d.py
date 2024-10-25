import numpy as np
def v2d(vrms):
    dBm=10*np.log10(20)+20*np.log10(vrms)
    return dBm