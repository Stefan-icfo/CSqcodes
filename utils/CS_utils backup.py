import numpy as np
import math
import time
from qcodes.dataset import Measurement, new_experiment
import qcodes as qc
import inspect

import pandas as pd


#define constants
hbar = 1.054571817e-34  # Reduced Planck constant (hbar) in Joule seconds (J·s)
kB = 1.380649e-23       # Boltzmann constant (kB) in Joules per Kelvin (J/K)
kB_rad=kB/hbar

def get_metadata(meas_id):
    qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1.db'
    experiments=qc.experiments()
    dataset=qc.load_by_id(meas_id)
    print(dataset.metadata)
    


def in_range_2d(point, x_range, y_range):
    x, y = point
    return x_range[0] <= x <= x_range[1] and y_range[0] <= y <= y_range[1]

def get_var_name(var):
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

def save_metadata_var(dataset,varnamelist,varlist):
       for varname,var in zip(varnamelist,varlist):
            dataset.add_metadata(varname[0],var)
        #print(temp_name[0])
        #print(temp_name[1])
        

def breit_wigner_fkt(x, peak_V, gamma ,peak_G,offset=0):
                return peak_G*(gamma**2 / (gamma**2 + ((x-peak_V)**2)))+offset

def breit_wigner_detuning(G, peak_G, gamma): #plus/minus
                return gamma*np.sqrt(peak_G/G-1)

def breit_wigner_derivative_analytical(x, peak_V, gamma, peak_G):
    A = peak_G
    V = peak_V
    # Analytical derivative of the Breit-Wigner function
    derivative = - (2 * A * (x - V) * gamma**2) / ((gamma**2 + (x - V)**2) ** 2)
    return derivative

def lorentzian_fkt(x, peak_V, gamma, peak_G):
  
    # Lorentzian function
    lorentzian = peak_G * (gamma / (np.pi * (gamma**2 + (x - peak_V)**2)))
    
    
    return lorentzian

def lorentzian_fkt_w_area(x, peak_V, gamma, peak_G):
  
    # Lorentzian function
    lorentzian = peak_G * (gamma / (np.pi * (gamma**2 + (x - peak_V)**2)))
    
    # Calculate the area under the Lorentzian peak
    # Area = amplitude * gamma * pi
    area = np.pi * peak_G * gamma
    
    return lorentzian, area

def zurich_phase_voltage_conductance(measured_value, vsdac, gain_RT = 200 ,gain_HEMT = 5.64,Z_tot = 7521):
                x = measured_value['x'][0] #SF: COMMENTED OUT 
                y = measured_value['y'][0]#SF: COMMENTED OUT
                #xy_complex = measured_value
                xy_complex = complex(x,y)
                v_r = np.absolute(xy_complex)
                theta = np.angle(xy_complex)
                    
                #G calculation
                I = v_r/(gain_RT*gain_HEMT*Z_tot)
                G = 1/((vsdac/I)-Z_tot)
                return theta, v_r, G

def zurich_phase_voltage_current_conductance(measured_value, vsdac, gain_RT = 200 ,gain_HEMT = 5.64,Z_tot = 7521):
                x = measured_value['x'][0] #SF: COMMENTED OUT 
                y = measured_value['y'][0]#SF: COMMENTED OUT
                #xy_complex = measured_value
                xy_complex = complex(x,y)
                v_r = np.absolute(xy_complex)
                theta = np.angle(xy_complex)
                    
                #G calculation
                I = v_r/(gain_RT*gain_HEMT*Z_tot)
                G = 1/((vsdac/I)-Z_tot)
                return theta, v_r, I, G

def zurich_phase_voltage_current_conductance_compensate(measured_value, vsdac, x_avg, y_avg, gain_RT = 200 ,gain_HEMT = 5.64,Z_tot = 7521):
                x = measured_value['x'][0]-x_avg #SF: COMMENTED OUT 
                y = measured_value['y'][0]-y_avg#SF: COMMENTED OUT
                #xy_complex = measured_value
                xy_complex = complex(x,y)
                v_r = np.absolute(xy_complex)
                theta = np.angle(xy_complex)
                    
                #G calculation
                I = v_r/(gain_RT*gain_HEMT*Z_tot)
                G = 1/((vsdac/I)-Z_tot)
                return theta, v_r, I, G

def zurich_x_y_avg(measured_parameter, tc=100e-3,avg_nr=100):
            x_sum=0
            y_sum=0
            for n in range(avg_nr):
                time.sleep(tc)
                measured_value=measured_parameter()
                x_sum+= measured_value['x'][0]  
                y_sum+= measured_value['y'][0]#SF: COMMENTED OUT
            x_avg=x_sum/avg_nr
            y_avg=y_sum/avg_nr     
            return x_avg,y_avg

def zurich_working(measured_parameter, tc=100e-3,avg_nr=10,cutoff_x=4e-6,cutoff_y=4e-6):#zurich.demods.demods0.sample
            x_sum=0
            y_sum=0
            for n in range(avg_nr):
                time.sleep(tc)
                measured_value=measured_parameter()
                x_sum+= measured_value['x'][0]  
                y_sum+= measured_value['y'][0]#SF: COMMENTED OUT
            x_avg=x_sum/avg_nr
            y_avg=y_sum/avg_nr 
            if abs(x_avg)<cutoff_x and abs(y_avg)<cutoff_y:
                return False
            else:
                return True


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
