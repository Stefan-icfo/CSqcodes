import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import re

# ---------------------
# CONFIGURATION
# ---------------------
qc.config["core"]["db_location"] = (
    r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v11.db"
)

# Choose a background run and a list of data runs
background_id = 2556#2885
run_ids = list(range(2316, 2556, 2))
#run_ids = list(range(2558, 2798, 2))

# Dataset keys
signal_key = "avg_avg_psd_nodrive"   # PSD channel (W/Hz)
freq_key   = "freq_param"             # frequency axis (Hz)

# ------------- helpers -------------
def load_psd(run_id):
    """
    Load one 1D PSD trace from a QCoDeS run:
    returns f_MHz (sorted ascending), P (same order), ds
    """
    ds = qc.load_by_id(run_id)
    pd = ds.get_parameter_data()

    if signal_key not in pd:
        raise KeyError(f"'{signal_key}' not found in run {run_id}. "
                       f"Available keys: {list(pd.keys())}")
    blk = pd[signal_key]
    if freq_key not in blk:
        # show user what setpoints exist
        raise KeyError(f"'{freq_key}' not found in run {run_id}. "
                       f"Available setpoints: {list(blk.keys())}")

    y = np.asarray(blk[signal_key]).ravel().astype(float)
    x = np.asarray(blk[freq_key]).ravel().astype(float) / 1e6  # Hz -> MHz

    # keep only finite, sort by frequency
    m = np.isfinite(x) & np.isfinite(y)
    if not np.any(m):
        raise ValueError(f"All values are non-finite in run {run_id}.")
    x, y = x[m], y[m]
    idx = np.argsort(x)
    return x[idx], y[idx], ds

def extract_gcs_voltage(name):
    """
    Parse 'gcs=...' from experiment name. Returns float or None.
    Example matches: gcs=0.842, gcs=0.8420, etc.
    """
    match = re.search(r"gcs=([+-]?\d+(?:\.\d+)?)", str(name))
    return float(match.group(1)) if match else None

def to_mV_if_looks_like_volts(v_array):
    """
    If gate values look like volts (<~5 in magnitude), convert to mV for plotting.
    Returns (converted_array, unit_label_str)
    """
    arr = np.asarray(v_array, dtype=float)
    if arr.size == 0:
        return arr, "arb."
    if np.nanmax(np.abs(arr)) < 5.0:
        return arr * 1e3, "mV"
    return arr, "V"

def safe_extent_from_axes(x, y, name_x="x", name_y="y"):
    x = np.asarray(x); y = np.asarray(y)
    if x.size == 0:
        raise ValueError(f"{name_x} axis is empty; check your filters/keys.")
    if y.size == 0:
        raise ValueError(f"{name_y} axis is empty; check your filters/keys.")
    return [np.nanmin(x), np.nanmax(x), np.nanmin(y), np.nanmax(y)]

# ---------------------
# LOAD BACKGROUND
# ---------------------
x_bg, P_bg, _ = load_psd(background_id)
# Use voltage-domain subtraction robustly: V = sqrt(P clipped)
V_bg = np.sqrt(np.clip(P_bg, 0.0, None))

# ---------------------
# PROCESS MEASUREMENTS
# ---------------------
xs = []          # list of frequency axes (MHz) for each run
spectra = []     # background-subtracted PSD (W/Hz), computed via V-domain subtraction
vgates = []      # parsed gate values (as stored, typically in V)
kept_runs = []
skipped = []

for run_id in run_ids:
    try:
        x, P, ds = load_psd(run_id)

        # Interpolate background voltage onto this run's x
        V = np.sqrt(np.clip(P, 0.0, None))
        V_bg_interp = np.interp(x, x_bg, V_bg, left=V_bg[0], right=V_bg[-1])

        V_corr = V - V_bg_interp
        P_corr = V_corr**2  # back to "power-like" units

        vg = extract_gcs_voltage(getattr(ds, "exp_name", ""))
        if vg is None:
            skipped.append((run_id, "V_gcs not found in exp_name"))
            continue

        xs.append(x)
        spectra.append(P_corr)
        vgates.append(vg)
        kept_runs.append(run_id)

    except Exception as e:
        skipped.append((run_id, str(e)))

if len(spectra) == 0:
    raise RuntimeError(
        "No usable runs. Reasons:\n" +
        "\n".join([f"  - {rid}: {msg}" for rid, msg in skipped]) +
        "\nCheck 'signal_key' and 'freq_key' and that exp_name contains 'gcs='."
    )

# ---------------------
# BUILD COMMON FREQ GRID (intersection) AND INTERPOLATE
# ---------------------
# Intersection of all frequency ranges
fmins = [np.nanmin(x) for x in xs]
fmaxs = [np.nanmax(x) for x in xs]
fmin = float(np.max(fmins))
fmax = float(np.min(fmaxs))
if not np.isfinite(fmin) or not np.isfinite(fmax) or fmax <= fmin:
    raise RuntimeError("Empty frequency intersection across runs.")

# Choose a common grid length = min length among runs, and linearly space in [fmin,fmax]
n_common = int(np.min([len(x) for x in xs]))
freqs = np.linspace(fmin, fmax, n_common)

# Interpolate each corrected spectrum to the common grid
spec2d = np.vstack([
    np.interp(freqs, x_i, s_i)
    for x_i, s_i in zip(xs, spectra)
])

# Convert gate axis (likely in V) to mV if appropriate
vgates = np.asarray(vgates, dtype=float)
vgates_plot, gate_unit = to_mV_if_looks_like_volts(vgates)

# Sort by gate for a nice image
sort_idx = np.argsort(vgates_plot)
spec2d = spec2d[sort_idx]
vgates_plot = vgates_plot[sort_idx]

# ---------------------
# PLOT
# ---------------------
extent = safe_extent_from_axes(freqs, vgates_plot, name_x="freqs(MHz)", name_y=f"gate({gate_unit})")

plt.figure(figsize=(10, 6))
plt.imshow(
    spec2d,
    aspect='auto',
    extent=extent,
    origin='lower',
    cmap='viridis',
    interpolation='nearest'
)
cbar = plt.colorbar()
cbar.set_label('PSD (bg-subtracted) [same units as input]')

plt.xlabel("Frequency (MHz)")
plt.ylabel(f"V_gcs ({gate_unit})")
first, last = kept_runs[0], kept_runs[-1]
plt.title(f"Peak map (background subtracted) — runs {first}..{last} (kept {len(kept_runs)}, skipped {len(skipped)})")
plt.tight_layout()
plt.show()

# Optional: print why runs were skipped
if skipped:
    print("\nSkipped runs:")
    for rid, msg in skipped:
        print(f"  - {rid}: {msg}")








