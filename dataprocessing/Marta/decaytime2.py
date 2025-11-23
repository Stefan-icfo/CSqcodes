import numpy as np
from qcodes.dataset import load_by_id, initialise_or_create_database_at

# ==================== USER SETTINGS ====================

db_path = r"C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v32_19_11_25.db"
run_id = 592

# output txt file
out_path = r"Z:\Users\Marta\demod_trace_run592.txt"

# =======================================================

def load_demod_trace(db_path, run_id):
    """
    Load time (t) and v_r from QCoDeS dataset using get_parameter_data('v_r')
    and flatten everything to 1D arrays.
    """
    initialise_or_create_database_at(db_path)
    ds = load_by_id(run_id)
    print(f"Loaded dataset with run_id = {run_id}")
    print(f"Dataset snapshot name: {ds.name}")

    param_data = ds.get_parameter_data('v_r')
    print("Keys in param_data:", list(param_data.keys()))

    dep_block_name = list(param_data.keys())[0]
    dep_block = param_data[dep_block_name]

    print("Keys inside dep_block:", list(dep_block.keys()))

    # find v_r key
    dep_name_candidates = [k for k in dep_block.keys() if 'v_r' in k.lower()]
    if not dep_name_candidates:
        raise RuntimeError("Could not identify the dependent v_r key in dep_block.")
    dep_name = dep_name_candidates[0]

    # setpoint (time) key = anything different from v_r
    setpoint_names = [k for k in dep_block.keys() if k != dep_name]
    if not setpoint_names:
        raise RuntimeError("Could not find a setpoint array for v_r.")
    setpoint_name = setpoint_names[0]

    print(f"Using setpoint key: {setpoint_name}")
    print(f"Using dependent key: {dep_name}")

    t_raw = np.asarray(dep_block[setpoint_name]).ravel()
    v_r_raw = np.asarray(dep_block[dep_name]).ravel()

    if t_raw.size == 0 or v_r_raw.size == 0:
        raise RuntimeError("Empty arrays for time or v_r; something is wrong with the dataset.")

    print(f"Length of time array: {t_raw.size}")
    print(f"Length of v_r array:  {v_r_raw.size}")

    return t_raw, v_r_raw


def save_trace_to_txt(t, v_r, out_path):
    """
    Save time and v_r to a txt file as two columns:
    column 1: time [s]
    column 2: v_r [V]
    """
    data = np.column_stack((t, v_r))
    header = "time_s\tv_r_V"
    np.savetxt(out_path, data, header=header)
    print(f"Saved {len(t)} points to:\n{out_path}")


def main():
    t, v_r = load_demod_trace(db_path, run_id)
    save_trace_to_txt(t, v_r, out_path)


if __name__ == "__main__":
    main()










