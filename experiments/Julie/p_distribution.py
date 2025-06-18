import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from qcodes.dataset import load_by_id as qc_load_by_id
import pandas as pd


RUN_IDS          = [2170, 2171, 2172, 2173,2174, 2175, 2176]   
TIME_WINDOW      = (0, 100)             # seconds to keep from each trace
POINTS_PER_AVG   = 5                    # box-average width
DROP_THRESHOLD   = 1.6e-5                 # V   (defines a drop / rise)
DURATION_CUTOFF  = 0.03                 # s   (ignore longer “high” intervals)
BIN_US           = 2000.0               # μs  (histogram bin width)


def extract_1d(run_id: int,
               data_name: str = "v_r",
               setpoint_name: str = "time_param"
               ) -> tuple[np.ndarray, np.ndarray]:

    ds   = qc_load_by_id(run_id)
    sp   = ds.get_parameter_data(data_name)[data_name][setpoint_name]
    data = ds.to_pandas_dataframe_dict()[data_name]
    return np.asarray(sp).flatten(), np.asarray(data).flatten()

def high_state_durations(v_r: np.ndarray,
                         t: np.ndarray,
                         drop_thr: float,
                         pts_per_avg: int,
                         dur_cut: float
                         ) -> np.ndarray:
 

    n      = len(v_r) // pts_per_avg
    v_avg  = v_r[:n*pts_per_avg].reshape(n, pts_per_avg).mean(axis=1)
    t_avg  = t[:n*pts_per_avg].reshape(n, pts_per_avg).mean(axis=1)

    dv     = np.diff(v_avg)
    start  = np.where(dv <= -drop_thr)[0] + 1
    end    = np.where(dv >=  drop_thr)[0] + 1

    paired_end = [end[end > s][0] for s in start if (end > s).any()]
    start      = start[:len(paired_end)]
    end        = np.asarray(paired_end)

    high_idx_start = end[:-1]
    high_idx_end   = start[1:]
    dur            = t_avg[high_idx_end] - t_avg[high_idx_start]
    return dur[dur <= dur_cut]


all_durations = []

for rid in RUN_IDS:
    t_abs, v_r  = extract_1d(rid)
    t_rel       = t_abs - t_abs[0]
    sel         = (t_rel >= TIME_WINDOW[0]) & (t_rel <= TIME_WINDOW[1])
    t_rel, v_r  = t_rel[sel], v_r[sel]

    dur = high_state_durations(v_r, t_rel,
                               DROP_THRESHOLD,
                               POINTS_PER_AVG,
                               DURATION_CUTOFF)
    all_durations.append(dur)


durations = np.concatenate(all_durations)

if durations.size == 0:
    raise RuntimeError("No high-state events detected in the selected runs.")


bin_s   = BIN_US * 1e-6
edges   = np.arange(0.0, durations.max() + bin_s, bin_s)
counts, _ = np.histogram(durations, bins=edges)
centers = edges[:-1] + 0.5 * bin_s

mask    = counts > 0          
slope, intercept = np.polyfit(centers[mask], np.log(counts[mask]), 1)
tau     = -1 / slope
A       = np.exp(intercept)
y_fit   = A * np.exp(-centers / tau)

plt.figure(figsize=(6,4))
plt.bar(centers, counts, width=bin_s, align="center", alpha=0.6, label="all runs")
plt.plot(centers, y_fit, "r--", lw=2,
         label=rf"$N(x)= {A:.1f}\,e^{{-x/{tau:.3f}\,\mathrm{{s}}}}$")
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.xlabel("High-state duration (s)")
plt.ylabel("Number of events")
plt.title(f"Histogram of {len(RUN_IDS)} runs " +
          f"({sum(len(d) for d in all_durations)} events)")
plt.legend()
plt.tight_layout()
plt.show()