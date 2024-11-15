import numpy as np
import time
from utils.CS_utils import lorentzian_fkt,thermal_CB_peak,breit_wigner_fkt,breit_wigner_detuning,zurich_phase_voltage_current_conductance_compensate,breit_wigner_derivative_analytical
from tqdm import tqdm
import scipy as scp
import matplotlib.pyplot as plt


e_C=1.60217e-19
hbar = 1.054571817e-34  # Reduced Planck constant (hbar) in Joule seconds (J·s)
kB = 1.380649e-23       # Boltzmann constant (kB) in Joules per Kelvin (J/K)
kB_rad=kB/hbar
kB_eV=kB/e_C


def GVG_simple(gate_sweep,measured_parameter,step_sleep_time,vsdac,x_avg,y_avg,reverse=False):
    Glist=[]
    Vlist=[]
    Rlist=[]
    Phaselist=[]
    Ilist=[]
        
        
    for gate_value in tqdm(gate_sweep, leave=False, desc='inner gate Sweep', colour = 'blue'): #fast axis loop (source) #TD: REVERSE DIRECTION
        gate_sweep.set(gate_value)
        time.sleep(step_sleep_time) # Wait 3 times the time contanst of the lock-in plus gate ramp speed
        measured_value = measured_parameter()
        theta_calc, v_r_calc, I,  G = zurich_phase_voltage_current_conductance_compensate(measured_value, vsdac,x_avg, y_avg)        
        Glist=Glist+[G]
        Vlist=Vlist+[v_r_calc]
        Ilist=Ilist+[I] 
        Phaselist=Phaselist+[theta_calc]
    if reverse: #if the sweep is reversed then the measurement lists have to be reversed too, since fast_axis_unreversible_list has the values for the unreversed sweep. double-check on a measurement if it really works as intended!
            Glist.reverse()
            Vlist.reverse()
            Ilist.reverse()
            Phaselist.reverse()
            #GIVlist.reverse()
            #VRlist.reverse()
            #PHASElist.reverse()
    return Glist,Vlist,Ilist,Phaselist

def fit_and_find_sitpos_singlepeak(gate_sweep,Glist,initial_guess=None, sitfraction="l_max_slope",return_full_fit_data=False):
     if initial_guess==None:
        initial_guess = [(gate_sweep[-1]+gate_sweep[0])/2, 5e-4, 5e-6,20e-9]#initial guess for peakV, Gamma,height,offset for first GVg
        #initial_guess_thermal = [(gate_sweep[-1]+gate_sweep[0])/2, 5e-4, 5e-6,20e-9]
     
      
     popt, pcov = scp.optimize.curve_fit(breit_wigner_fkt, gate_sweep, Glist, p0=initial_guess)
     
     if sitfraction=="l_max_slope":
          sitpos=popt[0]-popt[1]/np.sqrt(3)
     elif sitfraction=="r_max_slope":
          sitpos=popt[0]+popt[1]/np.sqrt(3)
     elif sitfraction=="max":
          sitpos=popt[0]          
     elif isinstance(sitfraction, (int, float)):
          sitpos=popt[0]-breit_wigner_detuning(popt[2]*sitfraction,popt[2],popt[1],popt[2])
     else:
        raise ValueError("sitpos must be a string or a number")
     
     if return_full_fit_data:
        return popt, pcov,sitpos
     else:
        return sitpos

def fit_and_find_sitpos_singlepeak_thermal(gate_sweep,Glist,initial_guess=None, sitfraction="l_max_slope",return_full_fit_data=False,alpha=0.15,Delta_E_V=0.091):
     
     
     if initial_guess==None:
        Glist_np = np.array(Glist)
        max_index_G = np.argmax(Glist_np)
        initial_guess = [gate_sweep[max_index_G],15e-6,100e-3]#initial guess for peakV, G_infty,T
     
     popt, pcov = scp.optimize.curve_fit(thermal_CB_peak, gate_sweep, Glist, p0=initial_guess)
     

     
     if sitfraction=="l_max_slope":
          sitpos=popt[0] - (2 * kB_eV * popt[2] / alpha) * np.arcsinh(1 / np.sqrt(2))
          slope = +popt[1] * (alpha * alpha * Delta_E_V) / (6 * kB_eV**2 * popt[2]**2 * np.sqrt(3))
     elif sitfraction=="r_max_slope":
          sitpos=popt[0] + (2 * kB_eV * popt[2] / alpha) * np.arcsinh(1 / np.sqrt(2))
          slope = -popt[1] * (alpha * alpha * Delta_E_V) / (6 * kB_eV**2 * popt[2]**2 * np.sqrt(3))
     elif sitfraction=="max":
          sitpos=popt[0]
          slope=0          
     #elif isinstance(sitfraction, (int, float)):
     #     sitpos=popt[0]-breit_wigner_detuning(popt[2]*sitfraction,popt[2],popt[1],popt[2])
     else:
        raise ValueError("sitpos must be a string")
     
     amp_at_max_slope = popt[1] * (alpha * Delta_E_V) / (6 * kB_eV * popt[2])


     if return_full_fit_data:
        return popt, pcov, amp_at_max_slope, slope, sitpos
     else:
        return amp_at_max_slope, slope, sitpos