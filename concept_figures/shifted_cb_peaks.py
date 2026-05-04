import numpy as np
import matplotlib.pyplot as plt

# --- Parameters ---
Gamma = 0.0008         # Full width at half maximum, FWHM (V)
V0_1 = 0.830-Gamma/2         # Center of first peak (V)
G_max = 12.0           # Peak conductance (µS)


# Second peak shifted by ~1 linewidth
V0_2 = V0_1 + Gamma

# --- Gate voltage axis ---
V = np.linspace(0.825, 0.835, 5000)

def breit_wigner(V, V0, G_max, Gamma):
    """Breit-Wigner (Lorentzian) resonance lineshape."""
    return G_max * (Gamma / 2)**2 / ((V - V0)**2 + (Gamma / 2)**2)

G1 = breit_wigner(V, V0_1, G_max, Gamma)
G2 = breit_wigner(V, V0_2, G_max, Gamma)

# --- Plot ---
fig, ax = plt.subplots(figsize=(6, 4.5))

ax.plot(V, G1, color='black',  lw=1.5,label=r'$M$ ')#, label=r'Peak 1: $V_0 = %.4f$ V' % V0_1)
ax.plot(V, G2, color='black',  lw=1.5, linestyle='--',
        label=r'$M+1$ ')

ax.set_xlabel(r'$V_{GS}$ (V)', fontsize=13)
ax.set_ylabel(r'$G$ ($\mu$S)', fontsize=13)
ax.set_xlim(0.825, 0.835)
ax.set_ylim(bottom=-0.3)
ax.legend(fontsize=15)

plt.tight_layout()
plt.savefig('breit_wigner.png', dpi=150)
plt.show()