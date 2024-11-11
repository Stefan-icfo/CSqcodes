import numpy as np
import time
from utils.CS_utils import breit_wigner_fkt,breit_wigner_detuning,zurich_phase_voltage_current_conductance_compensate,breit_wigner_derivative_analytical
from tqdm import tqdm
import scipy as scp


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

def fit_and_find_sitpos_singlepeak(gate_sweep,Glist,initial_guess=None, sitfraction="l_max_slope",return_full_fit_data=False,return_slope_and_sitamp=False):
     if initial_guess==None:
        Glist_np = np.array(Glist)
        max_index_G = np.argmax(Glist_np)
        initial_guess = [gate_sweep[max_index_G], 5e-4, 5e-6,20e-9]#initial guess for peakV, Gamma,height,offset for first GVg
     
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
     
     if return_full_fit_data or return_slope_and_sitamp:
         derivative=breit_wigner_derivative_analytical(sitpos,popt[0],popt[1],popt[2])
         amplitude_at_sitpos=breit_wigner_fkt(sitpos,popt[0],popt[1],popt[2],popt[3])
         
     if return_full_fit_data:
        return popt, pcov,derivative,amplitude_at_sitpos,sitpos
     elif return_slope_and_sitamp:
        return derivative,amplitude_at_sitpos,sitpos
     else:
        return sitpos
          