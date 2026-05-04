import numpy as np
import matplotlib.pyplot as plt

# --- Parameters ---
Gamma = 0.0008
V0_1 = 0.830 - Gamma / 8
G_max = 12.0
V0_2 = V0_1 + Gamma / 4

# --- Gate voltage axis ---
V = np.linspace(0.828, 0.832, 5000)

def breit_wigner(V, V0, G_max, Gamma):
    return G_max * (Gamma / 2)**2 / ((V - V0)**2 + (Gamma / 2)**2)

G1 = breit_wigner(V, V0_1, G_max, Gamma)
G2 = breit_wigner(V, V0_2, G_max, Gamma)

# --- Sinusoidal line ---
# Start from the left half-maximum point of the right peak (G2)
V_start = V0_2 - Gamma / 2 -Gamma/8       # left HWHM point of peak 2
V_sine = V[V >= V_start]           # only rightward from that point

sine_center = G_max / 2           # vertical centre = half maximum (~6 µS)
sine_amp    = 1.2                  # amplitude (µS)
sine_freq   = 1 / (Gamma * 0.6)   # several oscillations across the range
G_sine = sine_center + sine_amp * np.sin(2 * np.pi * sine_freq * (V_sine - V_start))

# --- Plot ---
fig, ax = plt.subplots(figsize=(6, 4.5))

ax.plot(V, G1, color='black', lw=1.5, linestyle='--', label=r'$\Phi_m = \pi/4$')
ax.plot(V, G2, color='black', lw=1.5, linestyle='--', label=r'$\Phi_m = 3\pi/4$')
ax.plot(V_sine, G_sine, color='red', lw=1.5, label=r'I_sig ($\omega_m$)')

ax.set_xlabel(r'$V_{GS}$ (V)', fontsize=13)
ax.set_ylabel(r'$G$ ($\mu$S)', fontsize=13)
ax.set_xlim(0.828, 0.832)
ax.set_ylim(bottom=-0.3)
ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig('breit_wigner_with_sine.png', dpi=150)
plt.show()