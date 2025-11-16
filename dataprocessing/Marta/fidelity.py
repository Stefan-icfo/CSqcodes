import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ========= USER =========
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v28_10_11_25.db"
run_id  = 430
bins    = 120
use_calibri = True
# ========================

if use_calibri:
    mpl.rcParams['font.family'] = 'Calibri'

initialise_or_create_database_at(db_path)
ds = load_by_id(run_id)
pdata = ds.get_parameter_data()

def _ravel(a):
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()

def extract_time_and_amplitude(pdata):
    """
    Trova un time trace e un canale di ampiezza:
    - preferisce 'R' (magnitude),
    - altrimenti sqrt(X^2+Y^2),
    - altrimenti usa X.
    Ritorna t (s), R_uV (µV).
    """
    # stampa qualche info per capire la struttura
    print("Top-level keys:", list(pdata.keys()))
    for dep_key, block in pdata.items():
        keys = list(block.keys())
        print(f"[Block '{dep_key}'] keys:", keys)

        # trova un asse tempo
        set_keys = [k for k in keys if k != dep_key]
        if not set_keys:
            continue
        t_key = next((k for k in set_keys if "time" in k.lower() or k.lower()=="t"), set_keys[0])
        t = _ravel(block[t_key])

        # ampiezza diretta?
        amp_keys = [k for k in keys if any(s in k.lower() for s in [" sample r", "sample r", " amplitude", "magnitude", "_r", ":r", "_ampl", "amplitude"])]
        if amp_keys:
            a_key = amp_keys[0]
            R = _ravel(block[a_key]) * 1e6  # V -> µV
            return t, R

        # X/Y -> magnitude
        x_keys = [k for k in keys if any(s in k.lower() for s in [" sample x", "sample x", "_x", ":x"]) and k != t_key]
        y_keys = [k for k in keys if any(s in k.lower() for s in [" sample y", "sample y", "_y", ":y"]) and k != t_key]
        if x_keys and y_keys:
            X = _ravel(block[x_keys[0]])
            Y = _ravel(block[y_keys[0]])
            R = np.sqrt(X**2 + Y**2) * 1e6  # V -> µV
            return t, R

        # fallback: usa X
        if x_keys:
            X = _ravel(block[x_keys[0]]) * 1e6
            return t, np.abs(X)

    raise RuntimeError("Nessun canale R/X/Y trovato per istogramma.")

t, R_uV = extract_time_and_amplitude(pdata)

# pulizia
m = np.isfinite(t) & np.isfinite(R_uV)
t, R_uV = t[m], R_uV[m]

print(f"N punti: {len(R_uV)}  |  min={np.nanmin(R_uV):.1f} µV, max={np.nanmax(R_uV):.1f} µV")

# soglia automatica (due-mean)
def twomeans_threshold(v):
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    t0 = np.median(v)
    for _ in range(100):
        lo = v[v <= t0]
        hi = v[v >  t0]
        if len(lo)==0 or len(hi)==0:
            break
        t1 = 0.5*(lo.mean()+hi.mean())
        if abs(t1 - t0) < 1e-9:
            t0 = t1; break
        t0 = t1
    return float(t0)

thr = twomeans_threshold(R_uV)
low  = R_uV[R_uV <= thr]
high = R_uV[R_uV >  thr]

print(f"Threshold = {thr:.2f} µV | LOW={len(low)} HIGH={len(high)}")

# range istogramma con margine
vmin, vmax = float(np.nanmin(R_uV)), float(np.nanmax(R_uV))
pad = 0.02*(vmax - vmin) if vmax>vmin else 1.0
h_range = (vmin - pad, vmax + pad)

# ---- PLOT ----
fig = plt.figure(figsize=(10, 4.8))

# (A) trace (una finestra per non saturare)
ax1 = fig.add_subplot(1, 2, 1)
Nmax = 6000
t_plot = t[-Nmax:] if len(t) > Nmax else t
r_plot = R_uV[-Nmax:] if len(R_uV) > Nmax else R_uV
ax1.plot(t_plot, r_plot, ".", ms=2, alpha=0.7)
ax1.axhline(thr, ls="--", lw=1.4, color="tab:red", label=f"thr = {thr:.1f} µV")
ax1.set_xlabel("time (s)", fontsize=13)
ax1.set_ylabel("amplitude (µV)", fontsize=13)
ax1.set_title("Demod amplitude vs time", fontsize=13)
ax1.legend(fontsize=10)
ax1.grid(alpha=0.3)

# (B) istogrammi
ax2 = fig.add_subplot(1, 2, 2)
# istogramma globale (grigio chiaro) — così vedi qualcosa anche se LOW/HIGH sono vuoti
ax2.hist(R_uV, bins=bins, range=h_range, alpha=0.30, color="0.6", label=f"All (n={len(R_uV)})", edgecolor="none")
# LOW / HIGH
if len(low):
    ax2.hist(low,  bins=bins, range=h_range, alpha=0.60, label=f"LOW (n={len(low)})", edgecolor="k")
if len(high):
    ax2.hist(high, bins=bins, range=h_range, alpha=0.60, label=f"HIGH (n={len(high)})", edgecolor="k")
ax2.axvline(thr, ls="--", lw=1.4, color="tab:red")
ax2.set_xlabel("amplitude (µV)", fontsize=13)
ax2.set_ylabel("counts", fontsize=13)
ax2.set_title("Amplitude histogram", fontsize=13)
ax2.legend(fontsize=10)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.show()
