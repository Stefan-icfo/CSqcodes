import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib as mpl

mpl.rcParams.update({
    "font.family": "Calibri",
    "font.size": 18,          # default per testo
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 18,
    # Usa Calibri anche dentro il mathtext (per V_g, ecc.)
    "mathtext.fontset": "custom",
    "mathtext.rm": "Calibri",
    "mathtext.it": "Calibri:italic",
    "mathtext.bf": "Calibri:bold",
})

# =============================
# User settings
# =============================
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v7.db"
#DB_PATH = r"Z:\Electromechanics\Projects\chargesensor\triton measurements\Raw_data\cd11_d7_c1\CD11_D7_C1.db"

RUN_ID  = 683
G_KEY   = "G"
ALPHA   = 0.25

# Physical constants
kB = 1.380649e-23  # J/K
e  = 1.602176634e-19  # C

# =============================
# Helper: derivative of Fermi
# =============================
def fermi_derivative(E, T):
    beta = 1/(kB*T)
    sech = 1/np.cosh(0.5*beta*E)
    return 0.25*beta*sech**2   # -df/dE

# =============================
# Convoluted model (numeric convolution)
# =============================
def cb_convoluted_numeric(Vg, V0, Gamma, T, A, Goff, alpha=ALPHA):
    Delta = alpha * e * (Vg - V0)  # energy detuning array
    Gvals = []
    for d in Delta:
        # integrand: Lorentzian(E-delta) * (-df/dE)
        integrand = lambda E: ( (Gamma/2/np.pi) / (E**2 + (Gamma/2)**2) ) * fermi_derivative(E-d, T)
        val, _ = quad(integrand, -50*kB*T, 50*kB*T)  # integrate over window ~100 kB T
        Gvals.append(val)
    return Goff + A * np.array(Gvals)

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

# =============================
# Load data
# =============================
qc.config["core"]["db_location"] = DB_PATH
ds = qc.load_by_id(RUN_ID)
Vg, G = get_xy_from_dataset(ds)

# =============================
# Fit
# =============================
V0_guess = Vg[np.argmax(G)]
Gamma_guess = 1e-22      # J (ordine di µeV)
T_guess = 0.2            # K
A_guess = 4e-21
Goff_guess = np.min(G)

p0 = [V0_guess, Gamma_guess, T_guess, A_guess, Goff_guess]

popt, pcov = curve_fit(cb_convoluted_numeric, Vg, G, p0=p0, maxfev=20000)
perr = np.sqrt(np.diag(pcov))

V0_fit, Gamma_fit, T_fit, A_fit, Goff_fit = popt
print("\n=== Convoluted fit (numeric) ===")
print(f"V0    = {V0_fit:.6g} V")
print(f"Γ     = {Gamma_fit/e*1e6:.3g} µeV")
print(f"T_e   = {T_fit*1e3:.3g} mK")
print(f"A     = {A_fit:.3g}")
print(f"Goff  = {Goff_fit:.3g}")

# =============================
# Plot
# =============================
# =============================
# Plot
# =============================
V_dense = np.linspace(Vg.min(), Vg.max(), 200)
G_fit = cb_convoluted_numeric(V_dense, *popt)

# Converte i dati in microSiemens per il plot
G_plot = G * 1e6
G_fit_plot = G_fit * 1e6

plt.figure(figsize=(9, 6))
plt.plot(Vg, G_plot, 'o', ms=3, color='k', label='data')
plt.plot(V_dense, G_fit_plot, '-', lw=2, color='red', label='fit')

plt.xlabel(" Voltage GCS (V)")
plt.ylabel("$G$ (µS)")
plt.title(f"Run {RUN_ID} — Convoluted Coulomb peak (numeric conv.)")
plt.legend()
plt.tight_layout()
plt.show()



