import matplotlib.pyplot as plt
import numpy as np
temp=[205,225,250,275,300,350,400,450,500,600,700,800,1000]
areas=np.array([430,498,529,690,675,820,766,839,776,643,340,347,180])
sensitivities=np.array([3.88,3.73,3.74,3.64,3.56,3.52,3.42,3.30,3.14,2.92,2.5,2.2,2.0])#from 700 on manual extraction
errors=np.array([60,70,70,90,86,108,95,109,105,100,74,76,53])

temp_ac=[300,400,500,600,800,1000]
autocorr=[0.133,0.107,0.106,0.055,0.0754,0.0554]

plt.plot(temp,areas,"g*")
plt.title("thermal areas")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("temp mK")
plt.ylabel("area fW/Hz")
plt.show()


plt.plot(temp,areas/sensitivities**2,"g*")
plt.title("sensitivity-scaled thermal areas")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("temp mK")
plt.ylabel("scaled area a.u.")
plt.show()


plt.errorbar(temp,areas, yerr=errors,fmt='g*', label='just areas w errorbars')
plt.title(" thermal areas w error")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("temp mK")
plt.ylabel("area fW/Hz")
plt.show()

scaled_errors = errors / sensitivities**2
plt.errorbar(temp,areas/sensitivities**2, yerr=scaled_errors,fmt='g*', label='just areas w errorbars')
plt.title(" thermal scaled areas w error")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("temp mK")
plt.ylabel("area fW/Hz")
plt.show()



plt.plot(temp_ac,autocorr,"r*")
plt.title("autocorr")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("temp mK")
plt.ylabel("ac time [ms]")
plt.show()