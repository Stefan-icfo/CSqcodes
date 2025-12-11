import os
import numpy as np
import time
from datetime import datetime
from zhinst.toolkit import Session
from IPython import get_ipython

# ---------------- USER SETTINGS ----------------
DEVICE_ID = "DEV20039"
DEMOD_CH = 3          # demod index in Python (0-based)
ORDER = 3
TC_START = 100e-9
TC_STOP  = 10.0e-6
TC_STEP  = 150e-9

N_REPEATS = 10
SCRIPT_PATH = "experiments/zurich_experiments/takedemodtimetrace_fidelity4.py"

SLEEP_AFTER_SET = 0.4   # s
SLEEP_BETWEEN_RUNS = 0.05# s
# -----------------------------------------------

ip = get_ipython()
if ip is None:
    raise RuntimeError("Run this in Jupyter/IPython (needed for %run).")

# Connect once
session = Session("localhost")
device = session.connect_device(DEVICE_ID)

# Configure only demod 3 once
device.demods[DEMOD_CH].enable(True)
device.demods[DEMOD_CH].order(ORDER)

# Build TC sweep (inclusive)
tcs = np.arange(TC_START, TC_STOP + 0.5 * TC_STEP, TC_STEP)

print(
    f"Sweeping demod {DEMOD_CH} TC from {TC_START*1e6:.3f} to {TC_STOP*1e6:.3f} µs "
    f"step {TC_STEP*1e9:.0f} ns  ->  {len(tcs)} points"
)
print(f"Each point runs: {N_REPEATS} × %run {SCRIPT_PATH}")

log = []

for k, tc in enumerate(tcs, start=1):
    # 1) Set TC on instrument
    device.demods[DEMOD_CH].timeconstant(float(tc))
    time.sleep(SLEEP_AFTER_SET)

    # 2) Read back TC (what the instrument actually uses)
    tc_read = float(device.demods[DEMOD_CH].timeconstant())

    print("\n" + "=" * 80)
    print(f"[{k}/{len(tcs)}] TC = {tc_read*1e6:.3f} µs   |   {datetime.now().isoformat(timespec='seconds')}")
    print("=" * 80)

    # 3) Run your script N times, passing TC via environment variable
    for r in range(1, N_REPEATS + 1):
        print(f"  Run {r}/{N_REPEATS} (TC={tc_read*1e6:.3f} µs)")
        os.environ["TC_NOW_S"] = str(tc_read)   # <- script reads this to build dataset name
        ip.run_line_magic("run", SCRIPT_PATH)
        time.sleep(SLEEP_BETWEEN_RUNS)

    log.append((datetime.now().isoformat(timespec="seconds"), tc_read))

print("\nDONE. Summary (timestamp, TC_us):")
for ts, tc_read in log:
    print(ts, f"{tc_read*1e6:.3f}")


