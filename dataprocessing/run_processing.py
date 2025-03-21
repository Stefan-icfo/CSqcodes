import numpy as np
import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from dataprocessing.processing_fkt import *
from utils.CS_utils import *




gain_RT=200
gain_HEMT=5.64
Z_tot=7521
vsdac=10.9e-6#rms at device

def V_to_G(V):

    I = V / (gain_RT * gain_HEMT * Z_tot)
    G = 1 / ((vsdac / I) - Z_tot)

    return G

dB_part2_loc=dB_location="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1_part2.db'

run_id1=994 #in dB part2
run_id2=821 #in dB part3

Vg1,power,V1_2d= extract_2d(run_id1, 
               data_2d_name="V_aux",
               setpoints1_name='QDAC_ch06_dc_constant_V',  # cs gate
               setpoints2_name='drive_mag',  # drive amp
               plot=False,
               dB_location=dB_part2_loc)

V1_1d=V1_2d[:,1]
maxindex1=np.argmax(V1_1d)
Vg_of_max1=Vg1[maxindex1]
print(f"peak of GVg1 at Vg= {Vg_of_max1} V")


Vg2,V2=extract_1d(run_id2, data_1d_name = "V_r", setpoint_name = 'QDAC_ch06_dc_constant_V',  plot = False)

maxindex2=np.argmax(V2)
Vg_of_max2=Vg2[maxindex2]
print(f"peak of GVg2 at Vg= {Vg_of_max2} V")

#correct voltages to be able to plot together
Vg1-=Vg_of_max1
Vg2-=Vg_of_max2

step_size1=Vg1[1]-Vg1[0]
step_size2=Vg2[1]-Vg2[0]

print(f"step size 1 is {step_size1*1e6} uV and step size 2 is {step_size2*1e6} uV")

plt.figure(1)
plt.plot(Vg1,V1_1d)
plt.plot(Vg2,V2)
plt.show()

#compute conductances
G1=V_to_G(V1_1d)
G2=V_to_G(V2)

plt.figure(2)
plt.title("comparison of Conductances")
plt.plot(Vg1*1e3,G1*1e6,label='peak used to measure mechanics on ICT')
plt.plot(Vg2*1e3,G2*1e6,label='good peak')
plt.legend()
plt.xlabel("Gate Voltage (mV)")  
plt.ylabel("Conductance (µS)")   

plt.show()

# Compute numerical derivatives
dG1_dVg = np.gradient(G1, Vg1)  # First derivative of V1 with respect to Vg1
dG2_dVg = np.gradient(G2, Vg2)  # First derivative of V2 with respect to Vg2

dG1_dVg=centered_moving_average(dG1_dVg,21)
dG2_dVg=centered_moving_average(dG2_dVg,21)

plt.figure(3)
plt.title("comparison of transconductances")
plt.plot(Vg1*1e3,dG1_dVg,label='peak used to measure mechanics on ICT')
plt.plot(Vg2*1e3,dG2_dVg,label='good peak')
plt.legend()
plt.xlabel("Gate Voltage (mV)")  
plt.ylabel("Transconductance (S/V)")
plt.show()

print(f"max transconductance1 is {max(dG1_dVg):3g} and max transconductance 2 is {max(dG2_dVg):3g} uV, fraction is {(max(dG2_dVg)/max(dG1_dVg)):3g}")