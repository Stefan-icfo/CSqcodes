import numpy as np
import matplotlib.pyplot as plt


def chi_m(Omega, Omega0, Gamma, m_eff):
    return 1.0 / (m_eff * (Omega0**2 - Omega**2) - 1j * m_eff * Gamma * Omega)


# --- Parameters ----------------------------------------------------------
thetadeg=30 #phase of homodyne
A_mag   = 2.0 # amplitude of sxf term
A_phase = np.deg2rad(0)          # phase in radians of sxf term
A       = A_mag * np.exp(1j * A_phase)
#Amp_Sxf_term=2


Omega0   = 2 * np.pi * 1.1e6
Gamma    = 2 * np.pi * 5.0e3       # wider (5 kHz) so asymmetry is visible on plot
m_eff    = 5.0e-13
S_FF     = 1.0e-32
S_xx_imp = 2.0e-32

# S_xF coupling: S_xF = A * sin(theta) * cos(theta) * S_xx_imp
#A        = Amp_Sxf_term                     # dimensionless amplitude, |A| <= 1
theta    = np.deg2rad(thetadeg)        # detection quadrature angle (radians)
S_xF     = A * np.sin(theta) * np.cos(theta) * S_xx_imp

# Linear frequency grid, a few linewidths around resonance
f   = np.linspace(1.08e6, 1.12e6, 4000)
Om  = 2 * np.pi * f
df_kHz = (f - Omega0 / (2 * np.pi)) / 1e3   # detuning in kHz for x-axis

# --- Components ----------------------------------------------------------
chi       = chi_m(Om, Omega0, Gamma, m_eff)
floor     = np.full_like(Om, S_xx_imp)
lorentz   = np.abs(chi)**2 * S_FF
cross     = 2 * np.real(chi) * S_xF                 # dispersive
total     = floor + lorentz + cross

# --- Plot ----------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# (a) the three components
ax1.plot(df_kHz, lorentz,                lw=1.8, label=r'$|\chi_m|^2\,S_{FF}$  (Lorentzian)')
ax1.plot(df_kHz, cross,                  lw=1.8, label=r'$2\,\mathrm{Re}[\chi_m]\,S_{xF}$  (dispersive)')
ax1.axhline(S_xx_imp, color='gray', ls=':', lw=1.2, label=r'$S_{xx}^{\mathrm{imp}}$  (floor)')
ax1.axhline(0, color='k', lw=0.5)
ax1.axvline(0, color='r', ls='--', lw=0.8, alpha=0.6)
ax1.set_xlabel(r'$(\Omega-\Omega_0)/2\pi$  (kHz)')
ax1.set_ylabel(r'PSD  [m$^2$/(rad/s)]')
ax1.set_title('Components')
ax1.legend(loc='upper right', fontsize=9)
ax1.grid(alpha=0.3)

# (b) total, plus no-correlation reference
total_no_corr = floor + lorentz
ax2.plot(df_kHz, total_no_corr, '--', lw=1.6, label='without $S_{xF}$ (symmetric)')
ax2.plot(df_kHz, total,                lw=2.0,
         label=fr'with $S_{{xF}}$: $A$={A:.2f}, $\theta$={np.rad2deg(theta):.0f}$^\circ$')
ax2.axvline(0, color='r', ls='--', lw=0.8, alpha=0.6)
ax2.set_xlabel(r'$(\Omega-\Omega_0)/2\pi$  (kHz)')
ax2.set_ylabel(r'$S_{xx}^{\mathrm{tot}}$  [m$^2$/(rad/s)]')
ax2.set_title('Total lineshape')
ax2.legend(loc='upper right', fontsize=9)
ax2.grid(alpha=0.3)

fig.tight_layout()
plt.show()