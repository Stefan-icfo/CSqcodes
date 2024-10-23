import numpy as np
import time
from utils.CS_utils import zurich_phase_voltage_current_conductance_compensate
from tqdm import tqdm

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