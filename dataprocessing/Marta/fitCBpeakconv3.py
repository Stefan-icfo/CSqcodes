import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from scipy.optimize import curve_fit
from scipy.integrate import quad

# =============================
# User settings
# =============================
DB_PATH = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v7.db"
RUN_ID  = 689
G_KEY   = "G"

# Lever arm (dimensionless)
ALPHA   = 0.30          # <- come richiesto
ENERGY_UNITS = "ueV"    # "ueV" | "eV" | "J"

# Physical constants
kB = 1.380649e-23  # J/K
e  = 1.602176634e-19  # C

# =============================
# Helpers
# =============================
def fermi_derivative(E, T):
    beta = 1.0/(kB*T)
    sech = 1.0/np.cosh(0.5*beta*E)
    return 0.25*beta*sech**2   # -df/dE  [1/J]

def cb_convoluted_numeric(Vg, V0, Gamma, T, A, Goff, alpha=ALPHA):
    """
    Convoluzione numerica:  ∫ dE  L(E) * (-df/dE)(E - Δ),  dove Δ = α e (Vg - V0)
    L(E) è Lorentziana normalizzata: (Γ/2π) / (E^2 + (Γ/2)^2)
    """
    Delta = alpha * e * (Vg - V0)  # [J]
    Gvals = np.empty_like(Delta)

    # Finestra di integrazione adattiva (copre dominio di f' e Lorentziana)
    # più robusta di ±50 kBT fisso
    def one_point(d):
        Ewin = 50.0 * max(kB*T, Gamma)     # mezzo range
        integrand = lambda E: ((Gamma/2.0/np.pi) / (E**2 + (Gamma/2.0)**2)) * fermi_derivative(E - d, T)
        val, _ = quad(integrand, -Ewin, +Ewin, limit=200)
        return val

    for i, d in enumerate(Delta):
        Gvals[i] = one_point(d)

    return Goff + A * Gvals

def get_xy_from_dataset(ds, g_key=G_KEY):
    pd = ds.get_parameter_data()
    if g_key in pd:
        block = pd[g_key]
    else:
        block = None
        for _, blk in pd.items():
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

def vg_to_energy(vg, v0, alpha=ALPHA, units=ENERGY_UNITS):
    dV = vg - v0
    if units == "J":
        return alpha * e * dV
    elif units == "eV":
        return alpha * dV
    elif units == "ueV":
        return 1e6 * alpha * dV
    else:
        raise ValueError("ENERGY_UNITS must be 'ueV', 'eV', or 'J'.")

def energy_unit_label(units):
    return {"J": "J", "eV": "eV", "ueV": r"$\mu$eV"}[units]

def estimate_initials(Vg, G, alpha=ALPHA):
    """Stime robuste per V0, FWHM, T_guess, Γ_guess, A_guess, Goff_guess."""
    # Centro al massimo dei dati
    V0 = Vg[np.argmax(G)]
    Gpk = G.max()
    Goff = np.percentile(G, 5)  # baseline robusta (coda bassa)
    half = Goff + 0.5*(Gpk - Goff)

    # FWHM sul gate (trova i due crossing sul lato sx/dx)
    # fallback se non trovati
    try:
        idx = np.arange(len(Vg))
        left  = idx[(Vg < V0) & (G <= half)]
        right = idx[(Vg > V0) & (G <= half)]
        iL = left[-1]
        iR = right[0]
        fwhm_V = Vg[iR] - Vg[iL]
        fwhm_J = alpha * e * fwhm_V
    except Exception:
        # fallback: usa 10% dell'intervallo come FWHM
        fwhm_V = 0.1*(Vg.max() - Vg.min())
        fwhm_J = alpha * e * fwhm_V

    # Ipotesi: regime termico dominante → FWHM ≈ 3.53 kBT
    T_guess = max(fwhm_J / (3.53 * kB), 0.02)  # non scendere sotto 20 mK
    # Γ_guess come frazione della FWHM in energia
    Gamma_guess = max(0.2 * fwhm_J, 1e-27)

    # Scala di ampiezza: ordine di (Gpk - Goff) * 4 kB T
    A_guess = (Gpk - Goff) * max(4.0 * kB * T_guess, 4.0 * Gamma_guess/np.pi)

    return V0, Gamma_guess, T_guess, A_guess, Goff

# =============================
# Load data
# =============================
qc.config["core"]["db_location"] = DB_PATH
ds = qc.load_by_id(RUN_ID)
Vg, G = get_xy_from_dataset(ds)

# Stime iniziali
V0_data, Gamma_g, T_g, A_g, Goff_g = estimate_initials(Vg, G, alpha=ALPHA)

# --- Selezione dati attorno al picco (±4×scala energetica) ---
eps_J_full = ALPHA * e * (Vg - V0_data)
Escale = max(kB*T_g, Gamma_g)
mask = np.abs(eps_J_full) <= 4.0 * Escale
Vg_fit = Vg[mask]
G_fit_data = G[mask]

# Parametri iniziali e vincoli
p0 = [V0_data, Gamma_g, T_g, A_g, Goff_g]
bounds = (
    [V0_data - 5e-3,   1e-27,  0.01,   0.0,      -np.inf],  # lower
    [V0_data + 5e-3,   5e-20,  5.0,    np.inf,    np.inf],  # upper
)

# Fit (pesatura più forte vicino al centro)
# sigma = 1/w -> residui pesati con 1/sigma^2 = w^2
epsJ_fit = ALPHA * e * (Vg_fit - V0_data)
w = 1.0 / (1.0 + (np.abs(epsJ_fit) / Escale)**2)   # campana centrata
sigma = 1.0 / (w + 1e-6)

popt, pcov = curve_fit(cb_convoluted_numeric, Vg_fit, G_fit_data,
                       p0=p0, bounds=bounds, sigma=sigma,
                       absolute_sigma=True, maxfev=20000)
perr = np.sqrt(np.diag(pcov))
V0_fit, Gamma_fit, T_fit, A_fit, Goff_fit = popt

print("\n=== Convoluted fit (numeric) ===")
print(f"V0(data max) = {V0_data:.6g} V")
print(f"V0(fitted)   = {V0_fit:.6g} V")
print(f"Γ            = {Gamma_fit/e*1e6:.3g} µeV")
print(f"T_e          = {T_fit*1e3:.3g} mK")
print(f"A            = {A_fit:.3g}")
print(f"Goff         = {Goff_fit:.3g}")

# Curve lisce per il plot e asse energia centrato sul picco dei dati
V_dense = np.linspace(Vg.min(), Vg.max(), 600)
G_model = cb_convoluted_numeric(V_dense, *popt)

eps_plot = vg_to_energy(V_dense, V0_data)      # x in ε−εF (zero al max dati)
eps_data = vg_to_energy(Vg, V0_data)

# =============================
# Plot
# =============================
plt.figure(figsize=(9, 6))
plt.plot(eps_data, G,  'o', ms=3, color='k', label='Data')
plt.plot(eps_plot, G_model, '-', lw=2, color='red', label='Convoluted fit')

plt.xlabel(r"$\epsilon - \epsilon_F$ (" + energy_unit_label(ENERGY_UNITS) + ")")
plt.ylabel("Conductance $G$ (S)")
plt.title(f"Run {RUN_ID} — Convoluted Coulomb peak (x: energy, zero at peak)")
plt.legend()
plt.tight_layout()
plt.show()

# Conversione comoda
print("\n--- Axis conversion ---")
print(f"Lever arm α = {ALPHA}")
print("ΔE [J]  = α·e·ΔVg")
print("ΔE [eV] = α·ΔVg")
print(f"1 mV ↔ {ALPHA*1e-3:.3g} eV = {ALPHA*1e3:.3g} µeV")

