import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit
from scipy.integrate import quad

# =============================
# User settings
# =============================
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v7.db"
RUN_ID  = 687
G_KEY   = "G"
ALPHA   = 0.25   # lever arm

# Physical constants
kB = 1.380649e-23  # J/K
e  = 1.602176634e-19  # C

# =============================
# Helper: derivative of Fermi
# (exactly your form)
# =============================
def fermi_derivative(E, T):
    beta = 1/(kB*T)
    sech = 1/np.cosh(0.5*beta*E)
    return 0.25*beta*sech**2   # -df/dE

# =============================
# Convoluted model (numeric convolution) — BLUE
# (exactly your implementation)
# =============================
def cb_convoluted_numeric(Vg, V0, Gamma, T, A, Goff, alpha=ALPHA):
    Delta = alpha * e * (Vg - V0)  # energy detuning array
    Gvals = []
    for d in Delta:
        # integrand: Lorentzian(E-delta) * (-df/dE)
        integrand = lambda E: ((Gamma/2/np.pi) / (E**2 + (Gamma/2)**2)) * fermi_derivative(E-d, T)
        val, _ = quad(integrand, -50*kB*T, 50*kB*T)  # integrate over window ~100 kB T
        Gvals.append(val)
    return Goff + A * np.array(Gvals)

# =============================
# Extra models (fit separately)
# =============================
def cb_thermal(Vg, V0, T, A, Goff, alpha=ALPHA):
    """Thermal broadening (sech^2)."""
    w = (kB*T)/(alpha*e)
    return Goff + A/np.cosh((Vg - V0)/(2*w))**2

def cb_lorentz(Vg, V0, gammaV, A, Goff):
    """Lorentzian in gate voltage (gammaV=HWHM in volts)."""
    return Goff + A*gammaV**2 / ((Vg - V0)**2 + gammaV**2)

# =============================
# Helpers
# =============================
def get_xy_from_dataset(ds, g_key=G_KEY):
    pd = ds.get_parameter_data()
    if g_key in pd:
        block = pd[g_key]
    else:
        block = None
        for name, blk in pd.items():
            if g_key in blk:
                block = blk
                break
        if block is None:
            raise KeyError(f"'{g_key}' not present in dataset.")
    vg_key_local = [k for k in block.keys() if k != g_key][0]
    y = np.asarray(block[g_key]).flatten()
    x = np.asarray(block[vg_key_local]).flatten()
    order = np.argsort(x)
    return x[order], y[order]

def r2_score(y, yfit):
    ss_res = np.sum((y - yfit)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    return 1 - ss_res/ss_tot if ss_tot > 0 else np.nan

# =============================
# Load data
# =============================
qc.config["core"]["db_location"] = DB_PATH
ds = qc.load_by_id(RUN_ID)
Vg, G = get_xy_from_dataset(ds)

# =============================
# 1) Fit — Convolution (numeric)  [BLUE]
#     (kept exactly like yours)
# =============================
V0_guess = Vg[np.argmax(G)]
Gamma_guess = 1e-22      # J (≈ µeV scale)
T_guess = 0.2            # K
A_guess = 1e-23          # your original scale
Goff_guess = np.min(G)

p0_conv = [V0_guess, Gamma_guess, T_guess, A_guess, Goff_guess]
popt_cv, pcov_cv = curve_fit(cb_convoluted_numeric, Vg, G, p0=p0_conv, maxfev=20000)
V0_cv, Gamma_cv, T_cv, A_cv, Goff_cv = popt_cv
G_cv_fit = cb_convoluted_numeric(Vg, *popt_cv)

# =============================
# 2) Fit — Thermal (sech^2)     [RED]
# =============================
p0_th = [V0_guess, 0.2, np.ptp(G), np.min(G)]
popt_th, pcov_th = curve_fit(cb_thermal, Vg, G, p0=p0_th, maxfev=20000)
V0_th, T_th, A_th, Goff_th = popt_th
G_th_fit = cb_thermal(Vg, *popt_th)

# =============================
# 3) Fit — Lorentzian           [GREEN]
# =============================
p0_lo = [V0_guess, (Vg.max()-Vg.min())/20.0, np.ptp(G), np.min(G)]
popt_lo, pcov_lo = curve_fit(cb_lorentz, Vg, G, p0=p0_lo, maxfev=20000)
V0_lo, gammaV_lo, A_lo, Goff_lo = popt_lo
G_lo_fit = cb_lorentz(Vg, *popt_lo)

# =============================
# Print results
# =============================
print("\n=== Convoluted fit (numeric, BLUE) ===")
print(f"V0    = {V0_cv:.6g} V")
print(f"Γ     = {Gamma_cv/e*1e6:.3g} µeV")
print(f"T_e   = {T_cv*1e3:.3g} mK")
print(f"A     = {A_cv:.3g} (S·J)")
print(f"Goff  = {Goff_cv:.3g} S")
print(f"R^2   = {r2_score(G, G_cv_fit):.5f}")

print("\n=== Thermal fit (sech^2, RED) ===")
print(f"V0    = {V0_th:.6g} V")
print(f"T     = {T_th*1e3:.3g} mK")
print(f"A     = {A_th:.3g} S")
print(f"Goff  = {Goff_th:.3g} S")
print(f"R^2   = {r2_score(G, G_th_fit):.5f}")

print("\n=== Lorentzian fit (GREEN) ===")
print(f"V0      = {V0_lo:.6g} V")
print(f"gammaV  = {gammaV_lo:.6g} V (HWHM)")
print(f"FWHM_V  = {2*gammaV_lo:.6g} V")
print(f"A       = {A_lo:.3g} S")
print(f"Goff    = {Goff_lo:.3g} S")
print(f"R^2     = {r2_score(G, G_lo_fit):.5f}")

# =============================
# Plot all together
# =============================
V_dense = np.linspace(Vg.min(), Vg.max(), 600)
plt.figure(figsize=(9, 6))
plt.plot(Vg, G, 'o', ms=3, color='k', label='Data')
plt.plot(V_dense, cb_convoluted_numeric(V_dense, *popt_cv), 'b-', lw=2, label='Convolution (Lorentz ⊗ Fermi)')
plt.plot(V_dense, cb_thermal(V_dense, *popt_th), 'r-', lw=2, label='Thermal (sech$^2$)')
plt.plot(V_dense, cb_lorentz(V_dense, *popt_lo), 'g-', lw=2, label='Lorentzian')
plt.xlabel("Gate voltage $V_g$ (V)")
plt.ylabel("Conductance $G$ (S)")
plt.title(f"Run {RUN_ID} — Coulomb peak fits (α = {ALPHA})")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

