import numpy as np
import matplotlib.pyplot as plt
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ==================== USER SETTINGS ====================
db_path = r"C:\Users\LAB-nanooptomechanic\Documents\MartaStefan\CSqcodes\Data\Raw_data\CD12_B5_F4v37_26_11_25.db"
run_start, run_end = 520, 590

BURST_DURATION = 1.0
t_center = 4.0

# finestra ringdown dopo shift (per la media)
tmin_shift, tmax_shift = -1.0e-3, 5.0e-3

DECIMATE = 2
SMOOTH_PTS = 11

# griglia comune per media (0.001s totale = 6 ms se -1..5ms)
Ngrid = 3000
tgrid = np.linspace(tmin_shift, tmax_shift, Ngrid)

# ==================== HELPERS ====================
def as_1d(a):
    return np.asarray(a).ravel()

def split_by_time_reset(t, min_len=10):
    t = np.asarray(t, dtype=float).ravel()
    dt = np.diff(t)
    reset_idx = np.where(dt < 0)[0] + 1
    boundaries = np.concatenate(([0], reset_idx, [len(t)]))
    out = []
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        if (b - a) >= min_len:
            out.append((a, b))
    return out

def extract_bursts(ds, min_len=20, decimate=1):
    d = ds.get_parameter_data()
    tsbs = np.asarray(d["time_since_burst_start"]["time_since_burst_start"])
    vr   = np.asarray(d["v_r"]["v_r"])
    xp   = np.asarray(d["x_proxy"]["x_proxy"])

    bursts = []
    if vr.ndim == 2:
        nbursts = vr.shape[0]
        sl = slice(None, None, decimate)
        for i in range(nbursts):
            bursts.append({
                "burst_i": i,
                "tsbs": np.asarray(tsbs[i, :], float)[sl],
                "vr":   np.asarray(vr[i, :], float)[sl],
                "xp":   np.asarray(xp[i, :], float)[sl],
            })
        return bursts

    tsbs = as_1d(tsbs).astype(float)
    vr   = as_1d(vr).astype(float)
    xp   = as_1d(xp).astype(float)

    slices = split_by_time_reset(tsbs, min_len=min_len)
    for i, (a, b) in enumerate(slices):
        sl = slice(a, b, decimate)
        bursts.append({"burst_i": i, "tsbs": tsbs[sl], "vr": vr[sl], "xp": xp[sl]})
    return bursts

def find_step_time(tsbs, x_proxy, smooth_pts=11):
    t = np.asarray(tsbs, float)
    p = np.asarray(x_proxy, float)
    if len(t) < 30:
        return np.nan

    k = int(smooth_pts)
    if k % 2 == 0:
        k += 1
    kernel = np.ones(k) / k
    p_s = np.convolve(p, kernel, mode="same")

    dp = np.gradient(p_s, t)
    idx = int(np.argmax(np.abs(dp)))  # step up OR down
    return float(t[idx])

def interp_nan_safe(x_new, x_old, y_old):
    x_old = np.asarray(x_old)
    y_old = np.asarray(y_old)
    good = np.isfinite(x_old) & np.isfinite(y_old)
    if np.count_nonzero(good) < 2:
        return np.full_like(x_new, np.nan, dtype=float)
    return np.interp(x_new, x_old[good], y_old[good], left=np.nan, right=np.nan)

# ==================== MAIN ====================
def main():
    initialise_or_create_database_at(db_path)

    sum_v = np.zeros_like(tgrid, dtype=float)
    cnt_v = np.zeros_like(tgrid, dtype=float)

    used = 0
    skipped = 0

    for rid in range(run_start, run_end + 1):
        try:
            ds = load_by_id(rid)
            bursts = extract_bursts(ds, min_len=20, decimate=DECIMATE)
        except Exception as e:
            print(f"[run {rid}] extract FAILED: {e}")
            continue

        if len(bursts) == 0:
            skipped += 1
            continue

        # scegli burst più vicino a 4 s
        burst_starts = np.array([b["burst_i"] * BURST_DURATION for b in bursts], float)
        bi = int(np.argmin(np.abs(burst_starts - t_center)))
        b = bursts[bi]

        tsbs = b["tsbs"]
        vr   = b["vr"]
        xp   = b["xp"]

        t0_rel = find_step_time(tsbs, xp, smooth_pts=SMOOTH_PTS)
        if not np.isfinite(t0_rel):
            skipped += 1
            continue

        tshift = tsbs - t0_rel

        mask = (tshift >= tmin_shift) & (tshift <= tmax_shift)
        if np.count_nonzero(mask) < 10:
            skipped += 1
            continue

        ts = tshift[mask]
        vs = vr[mask]
        order = np.argsort(ts)
        ts, vs = ts[order], vs[order]

        vi = interp_nan_safe(tgrid, ts, vs)

        good = np.isfinite(vi)
        if np.count_nonzero(good) < 10:
            skipped += 1
            continue

        sum_v[good] += vi[good]
        cnt_v[good] += 1
        used += 1

        if (rid - run_start) % 10 == 0:
            print(f"processed run {rid} | used={used} skipped={skipped}", flush=True)

    v_avg = sum_v / np.maximum(cnt_v, 1)

    print(f"DONE. used={used}, skipped={skipped}, median count per point = "
          f"{np.median(cnt_v[cnt_v>0]) if np.any(cnt_v>0) else 0}")

    # ===== ONE FINAL PLOT =====
    plt.figure(figsize=(8, 4))
    plt.plot(1e3 * tgrid, v_avg, linewidth=2)
    plt.axvline(0, ls="--")
    plt.xlabel("t' (ms) aligned (t0 from x_proxy step)")
    plt.ylabel("⟨v_r⟩ (avg)")
    plt.title(f"Average v_r (runs {run_start}-{run_end}), used={used}, skipped={skipped}")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()





