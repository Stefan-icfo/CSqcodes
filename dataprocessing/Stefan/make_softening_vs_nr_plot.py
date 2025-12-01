import matplotlib.pyplot as plt

plt.plot([16,20,24],[7.4,8.7,10.7],"g*")
plt.title("left softening")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("holenr")
plt.ylabel("softening [kHz]")
plt.show()


plt.plot([16,20,24,28],[6.5,7.8,13,13.8],"g*")
plt.title("avg softening")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("holenr")
plt.ylabel("softening [kHz]")
plt.show()