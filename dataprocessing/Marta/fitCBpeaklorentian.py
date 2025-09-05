import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit

# =============================
# User settings
# =============================
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v7.db"
RUN_ID  = 683
G_KEY   = "G"   # conductance for y-axis
GMAX_FIXED = 28.36e-6  # fixed maximum conductance in Siemens

# =============================
# Model
# =============================
def lorentz_fixedGmax(Vg, V0, gammaV, Goff, Gmax=GMAX_FIXED):
    """
    Lorentzian Coulomb peak with fixed maximum conductance Gmax.
    gammaV = HWHM in volts.
    """
    A = Gmax - Goff
    return Goff + A * (gammaV**2) / ((Vg - V0)**2 + gammaV**2)

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

# Initial guesses
V0_guess = Vg[np.argmax(G)]
gammaV_guess = (Vg.max()-Vg.min())/20
Goff_guess = np.min(G)

# =============================
# Fit Lorentzian with fixed Gmax
# =============================
popt, pcov = curve_fit(
    lorentz_fixedGmax, Vg, G,
    p0=[V0_guess, gammaV_guess, Goff_guess]
)
perr = np.sqrt(np.diag(pcov))

V0_fit, gammaV_fit, Goff_fit = popt
dV0, dgammaV, dGoff = perr
FWHM = 2*gammaV_fit
dFWHM = 2*dgammaV

print("\n=== Lorentzian fit (Gmax fixed at 28.36 µS) ===")
print(f"V0       = {V0_fit:.6g} ± {dV0:.2g} V")
print(f"gammaV   = {gammaV_fit:.6g} ± {dgammaV:.2g} V (HWHM)")
print(f"FWHM     = {FWHM:.6g} ± {dFWHM:.2g} V")
print(f"Goff     = {Goff_fit:.3g} ± {dGoff:.2g} S")
print(f"Gmax     = {GMAX_FIXED:.3g} S (fixed)")

# =============================
# Plot
# =============================
V_dense = np.linspace(Vg.min(), Vg.max(), 2000)
G_fit = lorentz_fixedGmax(V_dense, *popt)

plt.figure(figsize=(9, 6))
plt.plot(Vg, G, 'o', ms=3, color='k', label='Data')
plt.plot(V_dense, G_fit, '-', lw=2.0, color='green', label='Lorentzian fit (Gmax fixed)')
plt.xlabel("Gate voltage $V_g$ (V)")
plt.ylabel("Conductance $G$ (S)")
plt.title(f"Run {RUN_ID} — Lorentzian fit with $G_{{max}}={GMAX_FIXED*1e6:.2f}$ µS fixed")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
