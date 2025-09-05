import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit
from scipy.special import polygamma

# constants
kB = 1.380649e-23
e  = 1.602176634e-19
ALPHA = 0.25

def cb_convoluted_trigamma(Vg, V0, Gamma, T, A, Goff, alpha=ALPHA):
    Delta = alpha*e*(Vg-V0)
    z = 0.5 + (Gamma/2.0 + 1j*Delta) / (2*np.pi*kB*T)
    return Goff + A * np.real(polygamma(1, z))
