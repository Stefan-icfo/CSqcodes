import os, math, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ========= USER =========
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v28_10_11_25.db"
run_first, run_last = 422, 439
out_dir  = r"C:\Users\Public\Fidelity_422_439"
bins     = 120
use_calibri = True
# ========================

os.makedirs(out_dir, exist_ok=True)
if use_calibri:
    mpl.rcParams['font.family'] = 'Calibri'

initialise_or_create_database_at(db_path)

# ---------- helpers ----------
def _ravel(a):
    if isinstance(a, (list, tuple)) and len(a) and hasattr(a[0], "__len__"):
        return np.concatenate([np.asarray(x).ravel() for x in a])
    return np.asarray(a).ravel()

def extract_time_and_amplitude(ds):
    """
    Prova a trovare un time trace e un canale di ampiezza:
    - preferisce 'R' (magnitude),
    - altrimenti sqrt(X^2 + Y^2),
    - altrimenti |X|.
    Ritorna t (s), R_uV (µV). Lancia eccezione se non trova niente di utile.
    """
    pdata = ds.get_parameter_data()
    # scandisci i blocchi
    for dep_key, block in pdata.items():
        keys = list(block.keys())
        set_keys = [k for k in keys if k != dep_key]
        if not set_keys:
            continue
        # tempo
        t_key = next((k for k in set_keys if "time" in k.lower() or k.lower()=="t"), set_keys[0])
        t = _ravel(block[t_key])
        # ampiezza diretta?
        amp_keys = [k for k in keys if any(s in k.lower()
                                           for s in [" sample r","sample r"," amplitude","magnitude","_r",":r","_ampl","amplitude"])]
        if amp_keys:
            a_key = amp_keys[0]
            R = _ravel(block[a_key]) * 1e6  # V -> µV
            return t, R
        # X/Y -> magnitude
        x_keys = [k for k in keys if any(s in k.lower() for s in [" sample x","sample x","_x",":x"]) and k != t_key]
        y_keys = [k for k in keys if any(s in k.lower() for s in [" sample y","sample y","_y",":y"]) and k != t_key]
        if x_keys and y_keys:
            X = _ravel(block[x_keys[0]]); Y = _ravel(block[y_keys[0]])
            R = np.sqrt(X**2 + Y**2) * 1e6  # V -> µV
            return t, R
        # fallback: X
        if x_keys:
            X = _ravel(block[x_keys[0]]) * 1e6
            return t, np.abs(X)
    raise RuntimeError("Nessun canale R/X/Y trovato.")

def twomeans_threshold(v):
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    t0 = np.median(v)
    for _ in range(100):
        lo = v[v <= t0]; hi = v[v > t0]
        if len(lo)==0 or len(hi)==0: break
        t1 = 0.5*(lo.mean() + hi.mean())
        if abs(t1 - t0) < 1e-9: t0 = t1; break
        t0 = t1
    return float(t0)

def gaussian_fidelity(low, high):
    """Stima mu, sigma da LOW/HIGH; soglia bayesiana; overlap errors; fidelity; d'."""
    mu_lo, sd_lo = np.mean(low),  np.std(low,  ddof=1)
    mu_hi, sd_hi = np.mean(high), np.std(high, ddof=1)
    # soglia ottima (priori uguali) per gaussiane con varianze diverse
    if sd_lo <= 0 or sd_hi <= 0:
        return mu_lo, sd_lo, mu_hi, sd_hi, np.nan, np.nan, np.nan, np.nan
    a = 1/(sd_hi**2) - 1/(sd_lo**2)
    b = -2*(mu_hi/(sd_hi**2) - mu_lo/(sd_lo**2))
    c = (mu_hi**2)/(sd_hi**2) - (mu_lo**2)/(sd_lo**2) - 2*np.log(sd_hi/sd_lo)
    if abs(a) > 1e-15:
        roots = np.roots([a, b, c])
        # scegli la radice tra i due picchi
        cand = np.real(roots)
        thr_opt = cand[np.argmin(np.abs(cand - 0.5*(mu_lo+mu_hi)))]
    else:
        thr_opt = (mu_lo*sd_hi**2 - mu_hi*sd_lo**2) / (sd_hi**2 - sd_lo**2)
    # errori (cdf normale)
    from math import erf, sqrt
    def Phi(z): return 0.5*(1+erf(z/np.sqrt(2)))
    eps_lo = 1 - Phi((thr_opt - mu_lo)/sd_lo)   # low->high
    eps_hi =      Phi((thr_opt - mu_hi)/sd_hi)   # high->low
    F = 1 - 0.5*(eps_lo + eps_hi)
    dprime = abs(mu_hi - mu_lo)/math.sqrt(0.5*(sd_lo**2 + sd_hi**2))
    return mu_lo, sd_lo, mu_hi, sd_hi, thr_opt, eps_lo, eps_hi, F, dprime

def plot_hist(R, thr, low, high, title, png_path, bins=120):
    vmin, vmax = float(np.nanmin(R)), float(np.nanmax(R))
    pad = 0.02*(vmax - vmin) if vmax>vmin else 1.0
    h_range = (vmin - pad, vmax + pad)

    fig, ax = plt.subplots(figsize=(6.6,4.2))

