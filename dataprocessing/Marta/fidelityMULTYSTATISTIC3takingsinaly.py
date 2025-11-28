import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at
from scipy.optimize import curve_fit

# ==================== USER SETTINGS ====================
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v37_26_11_25.db"
run_first, run_last = 165,174
bins      = 160
use_calibri = True
units_in_volts = True                 # True se il dataset è in V; False se già in µV
out_dir   = r"C:\Users\Public\Fidelity_AGG_Y_signed_only"
save_png  = True
# Guess iniziale manuale (None = auto)
p_guess = None
# =======================================================

os.makedirs(out_dir, exist_ok=True)
if use_calibri:
    mpl.rcParams['font.family'] = 'Calibri'

# ------------------------ helpers ------------------------
def _ravel(a):
    """Flatten arrays and lists-of-arrays to 1D."""
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()

def extract_Y_signed_uV(ds, units_in_volts=True, debug=False):
    """
    Estrae strettamente il canale 'Y' (Immaginario) del demod.
    - NESSUN fallback su R / magnitude.
    - Se esistono più candidati Y, prende il primo.
    Ritorna un array 1D (µV). Lancia eccezione se non trova nulla.
    """
    scale = 1e6 if units_in_volts else 1.0
    pdata = ds.get_parameter_data()

    for dep_key, block in pdata.items():
        keys = list(block.keys())
        if debug:
            print(f"[dep='{dep_key}'] keys: {keys}")

        # pattern candidati per Y (abbastanza permissivi)
        for k in keys:
            kl = k.lower().strip()
            if (
                kl == "y" or
                kl.endswith(" sample y") or
                "sample y" in kl or
                kl.endswith(":y") or
                kl.endswith("_y") or
                " imag" in kl or
                "(y" in kl or
                kl.startswith("y ")
            ):
                Y = _ravel(block[k]) * scale
                Y = Y[np.isfinite(Y)]
                if Y.size == 0:
                    continue
                return Y

    raise RuntimeError("Signed Y channel not found. Stampa le keys e aggiusta i pattern in 'extract_Y_signed_uV'.")

def make_hist(y, bins=120):
    """Istogramma non normalizzato -> (centers, counts, bin_width)."""
    counts, edges = np.histogram(y, bins=bins)
    centers = 0.5*(edges[:-1] + edges[1:])
    bw = np.diff(edges).mean()
    return centers, counts.astype(float), bw

def two_means_threshold(v):
    """ISODATA threshold 1D per separare due lobi (per guess iniziali)."""
    v = np.asarray(v); v = v[np.isfinite(v)]
    if v.size == 0: return np.nan
    t = np.median(v)
    for _ in range(100):
        lo = v[v <= t]; hi = v[v > t]
        if len(lo)==0 or len(hi)==0: break
        t_new = 0.5*(lo.mean() + hi.mean())
        if abs(t_new - t) < 1e-9:
            t = t_new; break
        t = t_new
    return float(t)

# modello Gaussiano doppio
def f(x, *p):
    # p = (x0, x1, w1, h1, w2, h2)
    x0, x1, w1, h1, w2, h2 = p
    return (h1*np.exp(-(x-x0)**2/(2*w1**2)) +
            h2*np.exp(-(x-x1)**2/(2*w2**2)))

def g1(x, x0, w, h):
    return h*np.exp(-(x-x0)**2/(2*w**2))

# ---------- soglia ottimale per Gaussiane con ampiezze/larghezze diverse ----------
def optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2):
    """
    Trova la soglia dove le due Gaussiane sono uguali:
        h1 * exp(-(t-x0)^2/(2 w1^2)) = h2 * exp(-(t-x1)^2/(2 w2^2))
    (Bayes-optimal con prior uguali).
    Risolve una quadratica in t; sceglie la soluzione fra x0 e x1 se esiste,
    altrimenti la più vicina al midpoint.
    """
    if h1 <= 0 or h2 <= 0:
        return 0.5 * (x0 + x1)

    A = 1.0 / (w2**2) - 1.0 / (w1**2)
    B = -2.0 * x1 / (w2**2) + 2.0 * x0 / (w1**2)
    C = (x1**2) / (w2**2) - (x0**2) / (w1**2) - 2.0 * np.log(h2 / h1)

    if abs(A) < 1e-18:
        if abs(B) < 1e-18:
            return 0.5 * (x0 + x1)
        t = -C / B
        return float(t)

    roots = np.roots([A, B, C])
    roots = roots[np.isreal(roots)].real
    if roots.size == 0:
        return 0.5 * (x0 + x1)

    mid = 0.5 * (x0 + x1)
    between = [r for r in roots if min(x0, x1) <= r <= max(x0, x1)]
    if between:
        return float(min(between, key=lambda r: abs(r - mid)))
    return float(min(roots, key=lambda r: abs(r - mid)))

# ---------- formatter percentuale che non sfora mai 100 ----------
def format_percent_floor(p, decimals=4):
    """
    Formatta p (0–1) come percentuale con `decimals` decimali,
    usando floor invece che round (per non stampare >100.0000%).
    """
    p = max(0.0, min(1.0, float(p)))
    factor = 10**decimals
    v = np.floor(p * 100.0 * factor) / factor
    return f"{v:.{decimals}f}"

# ==================== AGGREGATE (SIGNED Y ONLY) ====================
initialise_or_create_database_at(db_path)

all_Y = []
for rid in range(run_first, run_last+1):
    try:
        ds = load_by_id(rid)
        # DEBUG: stampa le keys per i primi due run
        debug = (rid in (run_first, run_first+1))
        Y = extract_Y_signed_uV(ds, units_in_volts=units_in_volts, debug=debug)
        if Y.size:
            all_Y.append(Y)
            print(f"[run {rid}] collected {Y.size} samples (SIGNED Y)")
        else:
            print(f"[run {rid}] empty trace (skipped)")
    except Exception as e:
        print(f"[run {rid}] skipped: {e}")

if not all_Y:
    raise RuntimeError("No Y-signed data found across the requested runs.")

Y_all = np.concatenate(all_Y)
Y_all = Y_all[np.isfinite(Y_all)]
print(f"\nAggregated samples: {Y_all.size}  |  range = [{np.min(Y_all):.2f}, {np.max(Y_all):.2f}] µV")

if np.min(Y_all) >= 0:
    print("WARNING: no negative values found. You are likely not reading a signed Y channel.")
    print("Check the printed keys and adjust the matching rules inside 'extract_Y_signed_uV'.")

# ==================== HISTOGRAM + TWO-GAUSSIAN FIT ====================
ticks, counts, bw = make_hist(Y_all, bins=bins)

# Guess iniziale (auto) se p_guess è None
if p_guess is None:
    thr = two_means_threshold(Y_all)
    lo = Y_all[Y_all <= thr]; hi = Y_all[Y_all > thr]
    if lo.size < 10 or hi.size < 10:
        p30, p70 = np.percentile(Y_all, [30, 70])
        lo = Y_all[Y_all <= p30]; hi = Y_all[Y_all >= p70]
    x0 = np.mean(lo) if lo.size else np.min(Y_all)
    x1 = np.mean(hi) if hi.size else np.max(Y_all)
    if x0 > x1: x0, x1 = x1, x0
    w1 = np.std(lo, ddof=1) if lo.size > 1 else np.std(Y_all, ddof=1)
    w2 = np.std(hi, ddof=1) if hi.size > 1 else np.std(Y_all, ddof=1)

    def peak_height(xc):
        if len(ticks) == 0: return 1.0
        i = np.argmin(np.abs(ticks - xc))
        return max(counts[i], 1.0)
    h1 = peak_height(x0); h2 = peak_height(x1)
    p0 = np.array([x0, x1, max(w1, 1e-6), h1, max(w2, 1e-6), h2], dtype=float)
else:
    p0 = np.array(p_guess, dtype=float)

# Fit (bounds per tenere width/height positivi)
lb = [-np.inf, -np.inf, 1e-12, 0.0, 1e-12, 0.0]
ub = [ np.inf,  np.inf,  np.inf, np.inf,  np.inf, np.inf]
popt, pcov = curve_fit(f, ticks, counts, p0=p0, bounds=(lb, ub), maxfev=40000)
x0, x1, w1, h1, w2, h2 = popt
if x0 > x1:
    x0, x1 = x1, x0
    w1, w2 = w2, w1
    h1, h2 = h2, h1

# Soglia ottimale
threshold = optimal_threshold_equal_height(x0, x1, w1, h1, w2, h2)

# Fractions numeriche
arr = np.linspace(ticks.min(), ticks.max(), 200_001)
F0 = np.sum(g1(arr[arr < threshold], x0, w1, h1)) / np.sum(g1(arr, x0, w1, h1))
F1 = np.sum(g1(arr[arr > threshold], x1, w2, h2)) / np.sum(g1(arr, x1, w2, h2))
Fidelity = 0.5*(F0 + F1)
SNR = 2.0 * (x1 - x0)**2 / (w1**2 + w2**2)

print("\n=== Aggregated two-Gaussian fit (SIGNED Y ONLY) ===")
print(f"y0 = {x0:.2f} µV,  w1 = {w1:.2f} µV,  h1 = {h1:.1f}")
print(f"y1 = {x1:.2f} µV,  w2 = {w2:.2f} µV,  h2 = {h2:.1f}")
print(f"Optimal threshold = {threshold:.2f} µV")

print(f"F0(left)  = {format_percent_floor(F0, decimals=4)}%")
print(f"F1(right) = {format_percent_floor(F1, decimals=4)}%")
print(f"Fidelity  = {format_percent_floor(Fidelity, decimals=4)}%   |   SNR = {SNR:.2f}")

# --------------------------- PLOT ---------------------------
xplot = np.linspace(ticks.min(), ticks.max(), 2000)
fit_total = f(xplot, x0, x1, w1, h1, w2, h2)
g_low     = g1(xplot, x0, w1, h1)
g_high    = g1(xplot, x1, w2, h2)

plt.figure(figsize=(7.4, 4.8))
plt.bar(ticks, counts, width=bw, alpha=0.65,
        label=f"Aggregated runs {run_first}-{run_last}")
plt.plot(xplot, fit_total, "-", lw=2.0, label="Fit (G1+G2)")
plt.plot(xplot, g_low,  "--", lw=1.2, label="G1")
plt.plot(xplot, g_high, "--", lw=1.2, label="G2")
plt.axvline(threshold, ls="--", lw=1.4,
            label=f"opt. thr. = {threshold:.1f} µV")
plt.xlabel("amplitude Y (µV)", fontsize=13)
plt.ylabel("counts", fontsize=13)
plt.title("Aggregated histogram (SIGNED Y) and two-Gaussian fit", fontsize=13)
plt.legend(fontsize=8)
plt.tight_layout()

if save_png:
    fig_path = os.path.join(out_dir, f"agg_hist_2G_Ysigned_runs_{run_first}_{run_last}.png")
    plt.savefig(fig_path, dpi=220)
    print("Saved figure ->", fig_path)

plt.show()
