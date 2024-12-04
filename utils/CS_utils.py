import numpy as np
import math
import time
from qcodes.dataset import Measurement, new_experiment
import qcodes as qc
import inspect

import pandas as pd
import matplotlib.pyplot as plt

import scipy as scp


#------------------define constants----------------------
e_C=1.60217e-19
hbar = 1.054571817e-34  # Reduced Planck constant (hbar) in Joule seconds (J·s)
kB = 1.380649e-23       # Boltzmann constant (kB) in Joules per Kelvin (J/K)
kB_rad=kB/hbar
kB_eV=kB/e_C

#------------------metadata functions----------------------

def get_metadata(meas_id,return_data=False):
    qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1_part2.db'
    experiments=qc.experiments()
    dataset=qc.load_by_id(meas_id)
    print(dataset.metadata)
    if return_data:
          return dataset.metadata

def get_gate_Vs_from_metadata(meas_id,pre_str='qdac_ch0',post_str='_dc_constant_V',gate_nrs=[1,2,3,4,5]):
    qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1_part2.db'
    experiments=qc.experiments()
    dataset=qc.load_by_id(meas_id)
    metadata=dataset.metadata
    gates_dict = {f"gate{nr}": None for nr in gate_nrs}
    for gate_nr in gate_nrs:
        get_str=pre_str+f'{gate_nr}'+post_str
        gates_dict[f"gate{gate_nr}"] = metadata[get_str]
    
    return gates_dict

        
def get_var_name(var):
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

def save_metadata_var(dataset,varnamelist,varlist):
       for varname,var in zip(varnamelist,varlist):
            dataset.add_metadata(varname[0],var)
        #print(temp_name[0])
        #print(temp_name[1])



#------------------CB peak fitting functions----------------------

def thermal_CB_peak(x, peak_V, G_infty, T,Delta_E_V=0.091,alpha=0.15,offset=0):
    Delta_E=alpha*Delta_E_V
    delta=alpha*(x-peak_V)
    term = (Delta_E / (4 * kB_eV * T)) * np.cosh(delta / (2 * kB_eV * T)) ** -2
    G = G_infty * term
    return G

def thermal_CB_detuning(G, G_infty, T, Delta_E_V=0.091, alpha=0.15):

    Delta_E = alpha * Delta_E_V

    # Calculate the argument for arccosh
    argument = np.sqrt((Delta_E * G_infty) / (4 * kB_eV * T * G))
    
    # Ensure the argument is within valid range for arccosh
    if argument < 1:
        raise ValueError("Invalid argument for arccosh; G may be too large.")

    # Calculate delta using arccosh
    delta = 2 * kB_eV * T * np.arccosh(argument)

    # Solve for x
    x= delta / alpha
    return x

def dG_thermal_dx(x, peak_V, G_infty, T, Delta_E_V=0.091, alpha=0.15):
    # Compute Delta_E and delta
    Delta_E = alpha * Delta_E_V
    delta = alpha * (x - peak_V)
    
    # Compute the term inside the hyperbolic functions
    cosh_term = np.cosh(delta / (2 * kB_eV * T))
    sinh_term = np.sinh(delta / (2 * kB_eV * T))
    
    # Compute the derivative of G with respect to x
    dG_dx_value = -G_infty * (alpha * Delta_E) / (4 * kB_eV**2 * T**2) * cosh_term**-3 * sinh_term
    
    return dG_dx_value

def breit_wigner_fkt(x, peak_V, gamma ,peak_G,offset=0):
                return np.sqrt((peak_G*(gamma**2 / (gamma**2 + ((x-peak_V)**2))))**2+offset**2)



def breit_wigner_detuning(G, peak_G, gamma): #plus/minus
                return gamma*np.sqrt(peak_G/G-1)

def breit_wigner_derivative_analytical(x, peak_V, gamma, peak_G):
    A = peak_G
    V = peak_V
    # Analytical derivative of the Breit-Wigner function
    derivative = - (2 * A * (x - V) * gamma**2) / ((gamma**2 + (x - V)**2) ** 2)
    return derivative



def get_slope_at_given_sitpos_tunnel(gate_sweep,Glist,sitpos,initial_guess=None,return_full_fit_data=False,plot=True):
     if initial_guess==None:
        Glist_np = np.array(Glist)
        max_index_G = np.argmax(Glist_np)
        initial_guess = [gate_sweep[max_index_G], 5e-4, 5e-6,20e-9]#initial guess for peakV, Gamma,height,offset for first GVg
     
     popt, pcov = scp.optimize.curve_fit(breit_wigner_fkt, gate_sweep, Glist, p0=initial_guess)
     
     
     derivative=breit_wigner_derivative_analytical(sitpos,popt[0],popt[1],popt[2])
     amplitude_at_sitpos=breit_wigner_fkt(sitpos,popt[0],popt[1],popt[2],popt[3])

     if plot:
         plt.plot(gate_sweep,Glist)
         plt.plot(gate_sweep,breit_wigner_fkt(gate_sweep,popt[0],popt[1],popt[2],popt[3]))
         plt.plot([sitpos-100e-6,sitpos,sitpos+100e-6],[amplitude_at_sitpos-100e-6*derivative,amplitude_at_sitpos,amplitude_at_sitpos+100e-6*derivative])
         plt.show()
                  
         
     if return_full_fit_data:
        return popt, pcov,derivative,amplitude_at_sitpos
     
     else:
        return derivative,amplitude_at_sitpos 

def get_slope_at_given_sitpos_thermal(gate_sweep,Glist,sitpos,initial_guess=None,return_full_fit_data=False,plot=True):#in constr
     if initial_guess==None:
        Glist_np = np.array(Glist)
        max_index_G = np.argmax(Glist_np)
        initial_guess = [gate_sweep[max_index_G],15e-6,100e-3]#initial guess for peakV, G_infty,T
     
     popt, pcov = scp.optimize.curve_fit(thermal_CB_peak, gate_sweep, Glist, p0=initial_guess)
     
     
     derivative=dG_thermal_dx(sitpos,popt[0],popt[1],popt[2])
     #print(derivative)
     amplitude_at_sitpos=thermal_CB_peak(sitpos,popt[0],popt[1],popt[2])

     if plot:
         plt.plot(gate_sweep,Glist)
         plt.plot(gate_sweep,thermal_CB_peak(gate_sweep,popt[0],popt[1],popt[2]))
         plt.plot([sitpos-100e-6,sitpos,sitpos+100e-6],[amplitude_at_sitpos-100e-6*derivative,amplitude_at_sitpos,amplitude_at_sitpos+100e-6*derivative])
         plt.show()
                  
         
     if return_full_fit_data:
        return popt, pcov,derivative,amplitude_at_sitpos
     
     else:
        return derivative,amplitude_at_sitpos 


#---------mechanical fit functions------------------


def lorentzian_fkt(x, peak_V, gamma, peak_G, offset=0):
  
    # Lorentzian function
    lorentzian = np.sqrt((peak_G * (gamma / (np.pi * (gamma**2 + (x - peak_V)**2))))**2+offset**2)
    
    
    return lorentzian

def lorentzian_fkt_w_area(x, peak_V, gamma, peak_G):
  
    # Lorentzian function
    lorentzian = peak_G * (gamma / (np.pi * (gamma**2 + (x - peak_V)**2)))
    
    # Calculate the area under the Lorentzian peak
    # Area = amplitude * gamma * pi
    area = np.pi * peak_G * gamma
    
    return lorentzian, area

#---------------idt fit functions---------------------------------

def idt_perpendicular_angle(x1,y1,x2,y2):
        x=x2-x1
        y=y2-y1
        return math.atan(-y/x) 

def make_detuning_axis(x1,y1,x2,y2,delta=500e-6):
    beta=idt_perpendicular_angle(x1,y1,x2,y2)
    start_x=(x1+x2)/2+delta*math.cos(beta) 
    start_y=(y1+y2)/2+delta*math.sin(beta) 
    stop_x=(x1+x2)/2-delta*math.cos(beta) 
    stop_y=(y1+y2)/2-delta*math.sin(beta)
    return start_x,start_y,stop_x,stop_y

def make_detuning_axis_noncenter(x1,y1,x2,y2,delta=500e-6,xi=0,epsilon_0=0):
    beta=idt_perpendicular_angle(x1,y1,x2,y2)
    start_x=((1+xi)*x1+(1-xi)*x2)/2+(delta-epsilon_0)*math.cos(beta) 
    start_y=((1+xi)*y1+(1-xi)*y2)/2+(delta+epsilon_0)*math.sin(beta) 
    stop_x=((1+xi)*x1+(1-xi)*x2)/2-(delta-epsilon_0)*math.cos(beta) 
    stop_y=((1+xi)*y1+(1-xi)*y2)/2-(delta+epsilon_0)*math.sin(beta)
    return start_x,start_y,stop_x,stop_y


def make_detuning_axis_noncenterM(x1,y1,x2,y2,delta=500e-6,xi=0,epsilon_0=0):
    beta=idt_perpendicular_angle(x1,y1,x2,y2) 
    start_x=((1+xi)*x1+(1-xi)*x2)/2+(delta-epsilon_0)*math.cos(beta) 
    start_y=((1+xi)*y1+(1-xi)*y2)/2+(delta-epsilon_0)*math.sin(beta) 
    stop_x=((1+xi)*x1+(1-xi)*x2)/2-(delta+epsilon_0)*math.cos(beta) 
    stop_y=((1+xi)*y1+(1-xi)*y2)/2-(delta+epsilon_0)*math.sin(beta)
    return start_x,start_y,stop_x,stop_y
    
def make_detuning_axis_noncenterM2(x1,y1,x2,y2,delta=500e-6,xi=0,epsilon_0=0):
    beta=idt_perpendicular_angle(x1,y1,x2,y2)
    xm=(x1+x2)/2 + xi*((x1+x2)/2)
    ym=(y1+y2)/2 + xi*((y1+y2)/2)
    start_x=xm+(delta-epsilon_0)*math.cos(beta) 
    start_y=ym+(delta-epsilon_0)*math.sin(beta) 
    stop_x=xm-(delta+epsilon_0)*math.cos(beta) 
    stop_y=ym-(delta+epsilon_0)*math.sin(beta)
    
    return start_x,start_y,stop_x,stop_y,xm,ym
    
def idt_shape_energy(epsilon, t,Te):
    return (1/2) * (1 - ((epsilon) / np.sqrt((epsilon)**2 + (4*t**2)))) * np.tanh(np.sqrt((epsilon)**2 + (4*t**2)) / (2 * kB_rad * Te))

def idt_shape_voltage(detuning,leverarm, t,Te):#not done/correct yet!
    epsilon=detuning*leverarm
    M=idt_shape_energy(epsilon, t,Te)

#----------general functions--------------------------       


def in_range_2d(point, x_range, y_range):
    x, y = point
    return x_range[0] <= x <= x_range[1] and y_range[0] <= y <= y_range[1]

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)  # Cumulative sum of the array
    ret[n:] = ret[n:] - ret[:-n]     # Adjust cumulative sums to window sums
    ret = ret[n - 1:] / n            # Compute moving averages
    return np.concatenate((a[:n - 1], ret))  # Pad with original values for the first n-1 elements

def centered_moving_average(a, n=3):
    # Calculate offsets for even and odd window sizes
    offset = n // 2
    # Start by calculating the cumulative sum of the array
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    centered_avg = ret[n - 1:] / n
    # Use original values where the window doesn't fit
    if n % 2 == 0:
        centered_avg = np.concatenate((a[:offset], centered_avg, a[-offset + 1:]))
    else:
        centered_avg = np.concatenate((a[:offset], centered_avg, a[-offset:]))
    return centered_avg
 

