import numpy as np
import time
from utils.CS_utils import *
from tqdm import tqdm
import scipy as scp
import matplotlib.pyplot as plt
#from experiments.GVg_qdac_zurich_general import GVG_fun

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

def fit_and_find_sitpos_singlepeak_tunnel(gate_sweep,Glist,initial_guess=None, sitfraction="l_max_slope",return_full_fit_data=False):
     if initial_guess==None:
        Glist_np = np.array(Glist)
        max_index_G = np.argmax(Glist_np)
        initial_guess = [gate_sweep[max_index_G], 5e-4, 5e-6,20e-9]#initial guess for peakV, Gamma,height,offset for first GVg
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
     
     slope=breit_wigner_derivative_analytical(sitpos, popt[0], popt[1], popt[2])
     if return_full_fit_data:
        return popt, pcov,slope,sitpos
     else:
        return slope,sitpos

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
     elif isinstance(sitfraction, (int, float)):
          sitpos=popt[0]-thermal_CB_detuning(sitfraction*max(Glist_np), popt[1], popt[2])
          slope=dG_thermal_dx(sitpos, popt[0], popt[1], popt[2])
     else:
        raise ValueError("sitpos must be a string or a number")
     


     if return_full_fit_data:
        return popt, pcov,  slope, sitpos
     else:
        return slope, sitpos
     

def find_sitpos_from_avg_data(Vg,G_vals,sitfraction=0.5,data_avg_num=7,sit_side="left",return_avg_data=False):
     avg_G=centered_moving_average(G_vals,data_avg_num)
     max_avg=max(avg_G)
     deriv_avg=avg_G[:-1] - avg_G[1:]

     if isinstance(sitfraction, (int, float)):
          print("sitfraction is a number")
          if sit_side=="left":
               pos_idx = np.argmax(avg_G > max_avg*sitfraction)
          if sit_side=="left":
               pos_idx = len(avg_G) - 1 - np.argmax((avg_G[::-1] > max_avg * sitfraction))
          else:
              raise ValueError("sit_side is neither left nor right")
          sitpos=Vg[pos_idx]
          x=[Vg[pos_idx-data_avg_num:pos_idx+data_avg_num]]
          y=[G_vals[pos_idx-data_avg_num:pos_idx+data_avg_num]]
          if pos_idx == 0:
               raise ValueError("left_idx=0. probably no peak found")
          result=scp.stats.linregress(x,y)#try y*1e7, result/1e7
          slope=result.slope#
     elif sitfraction=="r_max_slope":
            rmax_id=np.argmax(deriv_avg)
            sitpos=(Vg[rmax_id]+Vg[rmax_id+1])/2
            slope=deriv_avg[rmax_id]/(Vg[rmax_id-1]-Vg[rmax_id])
            pos_idx=rmax_id
     elif sitfraction=="l_max_slope":
            lmax_id=np.argmin(deriv_avg)
            sitpos=(Vg[lmax_id]+Vg[lmax_id+1])/2
            slope=deriv_avg[lmax_id]/(Vg[lmax_id-1]-Vg[lmax_id])
            pos_idx=lmax_id
     elif sitfraction=="max":
            max_id=np.argmax(avg_G)
            sitpos=Vg[max_id]
            slope=0
            pos_idx=max_id
     else:
            raise ValueError("sitpos must be a string or a number")
     
     if return_avg_data:
          return avg_G,slope,sitpos,pos_idx
     else:
          return slope,sitpos,pos_idx