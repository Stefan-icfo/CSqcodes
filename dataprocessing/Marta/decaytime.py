import os
import numpy as np
import matplotlib.pyplot as plt
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ==================== USER SETTINGS ====================

# Path to your QCoDeS database
db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v32_19_11_25.db"

# Run ID that contains the demod time trace
run_id = 592

# =======================================================


def plot_demod_timetrace(db_path, run_id):
    """
    Load a QCoDeS dataset and plot the demodulated time trace v_r(t)
    for the selected run.

    This function:
      1. attaches to the database,
      2. loads the dataset,
      3. uses get_parameter_data('v_r') to obtain:
         - the setpoint array (time),
         - the dependent data (v_r),
      4. plots v_r vs time.
    """

    # Attach or create the database
    initialise_or_create_database_at(db_path)

    # Load the dataset by run ID
    ds = load_by_id(run_id)
    print(f"Loaded dataset with run_id = {run_id}")
    print(f"Dataset snapshot name: {ds.name}")

    # ----------------------------------------------------
    # Get parameter data ONLY for 'v_r'
    # This returns a nested dictionary structure of the form:
    #   {'v_r': {'time_param': array(...), 'v_r': array(...)}}
    # or similar (names can vary slightly).
    # ----------------------------------------------------
    param_data = ds.get_parameter_data('v_r')
    print("Keys in param_data:", list(param_data.keys()))

    # The top-level key is usually the dependent parameter name, e.g. 'v_r'
    dep_block_name = list(param_data.keys())[0]
    dep_block = param_data[dep_block_name]

    print("Keys inside dep_block:", list(dep_block.keys()))

    # Identify which entry is the dependent data (v_r) and which is the setpoint (time)
    # Conventionally, one key contains "v_r" and another is the setpoint (e.g. 'time_param')
    dep_name_candidates = [k for k in dep_block.keys() if 'v_r' in k.lower()]
    if not dep_name_candidates:
        raise RuntimeError("Could not identify the dependent v_r key in dep_block.")
    dep_name = dep_name_candidates[0]

    # The setpoint is any other key than the dependent
    setpoint_names = [k for k in dep_block.keys() if k != dep_name]
    if not setpoint_names:
        raise RuntimeError("Could not find a setpoint array for v_r.")
    setpoint_name = setpoint_names[0]

    print(f"Using setpoint key: {setpoint_name}")
    print(f"Using dependent key: {dep_name}")

    # Extract the raw numpy arrays.
    # QCoDeS often returns lists of arrays; we flatten everything to 1D.
    t_raw = np.asarray(dep_block[setpoint_name]).ravel()
    v_r_raw = np.asarray(dep_block[dep_name]).ravel()

    print(f"Length of time array: {t_raw.size}")
    print(f"Length of v_r array:  {v_r_raw.size}")
    if t_raw.size == 0 or v_r_raw.size == 0:
        raise RuntimeError("Empty arrays for time or v_r; something is wrong with the dataset.")

    # Optionally, print basic ranges
    print(f"time range: {t_raw[0]:.3e} s  ->  {t_raw[-1]:.3e} s")
    print(f"v_r range:  {np.min(v_r_raw):.3e} V  ->  {np.max(v_r_raw):.3e} V")

    # ----------------------------------------------------
    # Plot v_r vs time
    # ----------------------------------------------------
    plt.figure(figsize=(8, 4))
    plt.plot(t_raw * 1e6, v_r_raw * 1e3, 'o-', markersize=2, linewidth=0.8)
    plt.xlabel('time (µs)')
    plt.ylabel('v_r (mV)')
    plt.title(f'Demod time trace – run {run_id}')
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def main():
    plot_demod_timetrace(db_path, run_id)


if __name__ == "__main__":
    main()







