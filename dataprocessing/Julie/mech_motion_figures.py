# Auto-generated self-contained script.
# Reproduces the Kerr nonlinear lineshape fit from mech_motion_figures.ipynb.
#
# Requirements: qcodes, scipy, numpy, matplotlib, pandas.
# Expects the qcodes database at the path set below (update if it moves).
#
# Run with:
#     python mech_motion_figures.py
# or inside IPython:
#     %run mech_motion_figures.py

# -------- Data loading (from notebook cell 0) -------------------------------
# grc
#HERE
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd        
import qcodes as qc


run_ids = [1263, 1265, 1267, 1299, 1303, 1305, 1309, 1311, 1313, 1315,
           1319, 1321, 1323, 1325, 1329, 1331, 1333, 1335, 1339, 1341, 1343, 1345]

temperatures = [35.6766, 35.6766, 35.6766,
                50.0152, 50.0152, 50.0152,
                75.0531, 75.0531, 75.0531, 75.0531,
                100.3800, 100.3800, 100.3800, 100.3800,
                124.8840, 124.8840, 124.8840, 124.8840,
                149.9940, 149.9940, 149.9940, 149.9940]

qc.config["core"]["db_location"] = r"D:\databases CD11_D7_C1\CD11_D7_C1_part3.db"


def extract_1d(run_id,
               data_1d_name="Resistance",
               setpoint_name="single_gate",
               plot=False):
    ds       = qc.load_by_id(run_id)
    par_data = ds.get_parameter_data(data_1d_name)
    freq_raw = par_data[data_1d_name][setpoint_name]
    df_dict  = ds.to_pandas_dataframe_dict()
    psd_raw  = df_dict[data_1d_name]

    freq = np.asarray(freq_raw).flatten()
    psd  = np.asarray(psd_raw).flatten()

    if plot:
        plt.plot(freq, psd)
        plt.title(f"Measurement {run_id}")
        plt.ylabel("Avg PSD (W/Hz)")
        plt.xlabel("Frequency (MHz)")
        plt.show()

    return freq, psd


temp_data       = {}
frequency_array = None

for run_id, temp in zip(run_ids, temperatures):
    freq, psd = extract_1d(run_id,
                           data_1d_name="avg_avg_psd_nodrive",
                           setpoint_name="freq_param",
                           plot=False)

    if frequency_array is None:
        frequency_array = freq
    temp_data.setdefault(temp, []).append(psd)


# -------- Kerr fit (from notebook cell 4) -----------------------------------
from scipy.optimize import curve_fit
import scipy.constants as const
import numpy as np
import matplotlib.pyplot as plt

kB   = const.k
hbar = const.hbar

# =============================================================================
# Kerr lineshape fit — ROBUST two-stage version
#
# Why the previous fit was bad: curve_fit's Trust-Region solver stalls when
# parameters span 8+ orders of magnitude (f0 ~ 1.66e8 Hz, beta ~ 1e3 Hz,
# floor ~ 1e-15 W/Hz).  It also got trapped in local minima where a huge
# |beta| + mismatched A / Gamma reproduced the wings but missed the peak
# (hence the flat descending orange curves at 75 and 100 mK).
#
# Fix:
#   1. Internally work in kHz around a pivot frequency and in fW/Hz, so
#      every fit parameter is O(1)-O(10) and the Jacobian is well-scaled.
#   2. Reparameterize beta -> eta = beta / Gamma (dimensionless); bounded.
#   3. Stage 1: pure Lorentzian (beta=0) prefit to lock f0, Gamma, A, floor.
#   4. Stage 2: Kerr refinement, seeded from Stage 1, multi-start over a
#      handful of small eta values, with Gaussian peak-weighted residuals
#      so the fit matches the peak shape rather than the baseline.
# =============================================================================

# --- Kerr integrand weights (eps * exp(-eps), normalised) --------------------
_N_EPS    = 400
_EPS_GRID = np.linspace(0.0, 12.0, _N_EPS)
_W_EPS    = _EPS_GRID * np.exp(-_EPS_GRID)

def _kerr_scaled(df_kHz, df0, G_kHz, eta, A_fWkHz, floor_fW):
    """Kerr lineshape in SCALED units (kHz, fW/Hz). Returns fW/Hz."""
    beta_kHz = eta * G_kHz
    f_eps    = df0 + beta_kHz * _EPS_GRID
    d        = df_kHz[:, None] - f_eps[None, :]
    lor      = (G_kHz / np.pi) / (d*d + G_kHz*G_kHz)
    spec     = np.trapz(_W_EPS[None, :] * lor, _EPS_GRID, axis=1)
    return floor_fW + A_fWkHz * spec

def _lorentz_scaled(df_kHz, df0, G_kHz, A_fWkHz, floor_fW):
    d   = df_kHz - df0
    lor = (G_kHz / np.pi) / (d*d + G_kHz*G_kHz)
    return floor_fW + A_fWkHz * lor

# --- Physics-unit wrappers (used for plotting) -------------------------------
def kerr_model(f_hz, f0_hz, G_hz, beta_hz, A_hz, floor_fW):
    f_eps = f0_hz + beta_hz * _EPS_GRID
    d     = f_hz[:, None] - f_eps[None, :]
    lor   = (G_hz / np.pi) / (d*d + G_hz*G_hz)
    spec  = np.trapz(_W_EPS[None, :] * lor, _EPS_GRID, axis=1)
    return (floor_fW + A_hz * spec) * 1e-15

def lorentz_model(f_hz, f0_hz, G_hz, A_hz, floor_fW):
    d   = f_hz - f0_hz
    lor = (G_hz / np.pi) / (d*d + G_hz*G_hz)
    return (floor_fW + A_hz * lor) * 1e-15

# --- Initial guesses from smoothed spectrum ---------------------------------
def _smooth(y, k=7):
    if k < 3 or len(y) < 2*k:
        return y
    return np.convolve(y, np.ones(k)/k, mode="same")

def _guesses(df_kHz, psd_fw):
    ys   = _smooth(psd_fw, 7)
    peak = int(np.argmax(ys))
    df0  = float(df_kHz[peak])
    n_edge  = max(5, len(df_kHz) // 6)
    floor   = float(np.mean(np.concatenate([psd_fw[:n_edge], psd_fw[-n_edge:]])))
    peak_ht = max(ys[peak] - floor, 0.05 * ys[peak])
    half    = floor + peak_ht/2
    above   = ys > half
    if above.sum() >= 2:
        idx   = np.where(above)[0]
        G_kHz = 0.5 * (df_kHz[idx[-1]] - df_kHz[idx[0]])
    else:
        G_kHz = (df_kHz[-1] - df_kHz[0]) * 0.005
    A_fWkHz = peak_ht * np.pi * G_kHz
    return df0, G_kHz, A_fWkHz, floor

# --- Two-stage fit ----------------------------------------------------------
_ETA_TRIALS = (-1.5, -0.8, -0.3, 0.0, 0.3, 0.8, 1.5)

def fit_kerr(freq_hz, psd_WHz, *, n_fwhm=12, peak_weight=3.0, eta_max=2.5):
    """Two-stage Kerr lineshape fit.

    Returns dict with 'popt' = (f0_hz, Gamma_hz, beta_hz, A_hz, floor_fW),
    'perr' (1-sigma errors), 'rss' (weighted), plus diagnostics.
    """
    pivot  = float(np.median(freq_hz))
    df_kHz = (freq_hz - pivot) * 1e-3
    psd_fw = psd_WHz * 1e15

    # Stage 1 - Lorentzian prefit
    df0_g, G_g, A_g, fl_g = _guesses(df_kHz, psd_fw)
    span  = df_kHz[-1] - df_kHz[0]
    mask1 = np.abs(df_kHz - df0_g) < n_fwhm * G_g
    if mask1.sum() < 20:
        mask1 = np.ones_like(df_kHz, dtype=bool)
    lo = [df0_g - 0.5*span, G_g*0.05, 0.0,    0.0          ]
    hi = [df0_g + 0.5*span, G_g*20,   np.inf, fl_g*5 + 1   ]
    lor_popt, _ = curve_fit(_lorentz_scaled, df_kHz[mask1], psd_fw[mask1],
                            p0=[df0_g, G_g, A_g, fl_g], bounds=(lo, hi), maxfev=20000)

    # Stage 2 - Kerr refinement with peak-weighted residuals
    df0_0, G_0, A_0, fl_0 = lor_popt
    mask2 = np.abs(df_kHz - df0_0) < n_fwhm * G_0
    if mask2.sum() < 20:
        mask2 = np.ones_like(df_kHz, dtype=bool)
    df_w, psd_w = df_kHz[mask2], psd_fw[mask2]
    weight      = 1.0 + (peak_weight - 1.0) * np.exp(-((df_w - df0_0) / (2.0*G_0))**2)
    sigma       = 1.0 / weight

    lo2 = [df0_0 - 0.3*span, G_0*0.25, -eta_max, A_0*0.1,  0              ]
    hi2 = [df0_0 + 0.3*span, G_0*4.0,   eta_max, A_0*10.0, fl_0*3 + 1     ]

    best = dict(popt=None, pcov=None, rss=np.inf, eta0=None)
    for eta0 in _ETA_TRIALS:
        try:
            popt, pcov = curve_fit(
                _kerr_scaled, df_w, psd_w,
                p0=[df0_0, G_0, eta0, A_0, fl_0],
                bounds=(lo2, hi2), sigma=sigma, absolute_sigma=False, maxfev=30000,
            )
            resid = psd_w - _kerr_scaled(df_w, *popt)
            rss   = float(np.sum((resid * weight) ** 2))
            if rss < best["rss"]:
                best.update(popt=popt, pcov=pcov, rss=rss, eta0=eta0)
        except (RuntimeError, ValueError):
            continue

    if best["popt"] is None:   # fall back to Lorentzian if every Kerr start failed
        df0_k, G_k, A_k, fl_k = lor_popt
        popt_s   = np.array([df0_k, G_k, 0.0, A_k, fl_k])
        perr_s   = np.full(5, np.nan)
        rss_used = np.nan
    else:
        popt_s   = best["popt"]
        rss_used = best["rss"]
        perr_s   = np.sqrt(np.diag(best["pcov"]))

    df0_k, G_k, eta, A_fWkHz, fl_fW = popt_s
    f0_hz, G_hz   = pivot + df0_k*1e3, G_k*1e3
    beta_hz, A_hz = eta * G_hz, A_fWkHz*1e3
    popt = np.array([f0_hz, G_hz, beta_hz, A_hz, fl_fW])

    df0_e, G_e_k, eta_e, A_e_k, fl_e = perr_s
    beta_e = np.sqrt((G_hz * eta_e)**2 + (eta * G_e_k*1e3)**2)
    perr   = np.array([df0_e*1e3, G_e_k*1e3, beta_e, A_e_k*1e3, fl_e])

    return dict(popt=popt, perr=perr, rss=rss_used,
                mask=mask2, eta0_best=best.get("eta0"))


# =============================================================================
# Apply fit per temperature and plot
# =============================================================================
n_temps = len(temp_data)
ncols   = 2
nrows   = int(np.ceil(n_temps / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(13, 4.5 * nrows))
axes = np.array(axes).flatten()

kerr_fit_results = {}
print(f"{'T (mK)':>8}  {'f0 (MHz)':>12}  {'Gamma (kHz)':>12}  "
      f"{'beta (Hz)':>11}  {'beta/G':>7}  {'K/2pi (Hz)':>11}  {'eta0':>6}  {'RSS':>11}")
print("-" * 94)

for ax_idx, temp in enumerate(sorted(temp_data.keys())):
    traces  = np.vstack(temp_data[temp])
    avg_psd = traces.mean(axis=0)
    freq    = frequency_array
    T_K     = temp * 1e-3

    res = fit_kerr(freq, avg_psd)
    f0_fit, G_fit, beta_fit, A_fit, Sf_fit = res["popt"]
    omega0_fit = 2 * np.pi * f0_fit
    K_hz       = beta_fit * hbar * omega0_fit / (kB * T_K)

    kerr_fit_results[temp] = dict(popt=res["popt"], perr=res["perr"],
                                  T_K=T_K, rss=res["rss"])

    f_fine    = np.linspace(freq[0], freq[-1], 3000)
    fit_curve = kerr_model(f_fine, *res["popt"]) * 1e15

    ax = axes[ax_idx]
    ax.plot(freq * 1e-6, avg_psd * 1e15, '.', ms=2, color='steelblue',
            alpha=0.6, label='Data')
    ax.plot(f_fine * 1e-6, fit_curve, '-',  lw=2.0, color='darkorange',
            label=(f'Kerr fit\n'
                   f'$K/\\Gamma={K_hz/G_fit:.3f}$'))
    ax.axvline(f0_fit * 1e-6, ls=':', lw=0.8, color='darkorange', alpha=0.6)
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('PSD (fW/Hz)')
    ax.set_title(f'T = {temp:.1f} mK')
    ax.legend(fontsize=14, loc='upper right')
    ax.ticklabel_format(axis='x', style='plain', useOffset=False)
    ax.grid(alpha=0.3)

    eta0_str = f"{res['eta0_best']:+.1f}" if res['eta0_best'] is not None else "---"
    print(f"{temp:>8.1f}  {f0_fit*1e-6:>12.5f}  {G_fit*1e-3:>12.3f}  "
          f"{beta_fit:>11.1f}  {beta_fit/G_fit:>7.2f}  "
          f"{K_hz:>11.2f}  {eta0_str:>6}  {res['rss']:>11.3e}")

for ax in axes[n_temps:]:
    ax.set_visible(False)

fig.suptitle('Kerr nonlinear lineshape fit - CNT resonator (robust two-stage)',
             fontsize=13)
fig.tight_layout()
plt.show()

# =============================================================================
# K/2pi vs temperature - should be constant if purely Kerr
# =============================================================================
if kerr_fit_results:
    T_arr = np.array(sorted(kerr_fit_results.keys()))
    K_arr = np.array([
        kerr_fit_results[t]["popt"][2]
        * hbar * 2 * np.pi * kerr_fit_results[t]["popt"][0]
        / (kB * kerr_fit_results[t]["T_K"])
        for t in T_arr
    ])
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.plot(T_arr, K_arr, 'o-', color='darkorange', lw=1.5)
    ax2.axhline(np.mean(K_arr), ls='--', color='k', lw=1,
                label=f'Mean $K/2\\pi$ = {np.mean(K_arr):.1f} Hz')
    ax2.set_xlabel('Bath temperature (mK)')
    ax2.set_ylabel('$K/2\\pi$ (Hz)')
    ax2.set_title('Kerr coefficient vs temperature\n(flat = intrinsic, sloped = T-dependent effect)')
    ax2.legend()
    ax2.grid(alpha=0.3)
    fig2.tight_layout()
    plt.show()
