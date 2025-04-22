
import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
from qcodes.dataset import load_by_id as qc_load_by_id  

# Fully automated code
# Set the run ids AND the value to calculate the detuning 


# List of main run IDs (for the 1D data extraction)
run_ids = [936, 988, 1040, 1092]
# Corresponding run IDs for temperature metadata (each temperature measurement is assumed to be run_id - 1)
# those run ids are also used to extract the charge sensor gate voltage
run_ids_temp = [run_id - 1 for run_id in run_ids]

def get_metadata(meas_id):
    """
    Retrieve and return the metadata for a single measurement ID.
    """
    dataset = qc_load_by_id(meas_id)
    return dataset.metadata

def extract_1d(run_id, data_1d_name="Resistance", setpoint_name='single_gate', plot=False):
    """
    Extract 1D data (x and y arrays) from the dataset corresponding to run_id.
    """
    dataset = qc_load_by_id(run_id)
    data_x = dataset.get_parameter_data(data_1d_name)
    setpoints_raw = data_x[data_1d_name][setpoint_name]
    setpoints_np = np.array(setpoints_raw)
    
    pdf_temp = dataset.to_pandas_dataframe_dict()
    data1d_raw = pdf_temp[data_1d_name]
    data1d_np = np.array(data1d_raw)
    
    if plot:
        plt.plot(setpoints_np, data1d_np)
        plt.title(f"Measurement {run_id}")
        plt.ylabel("Avg PSD W/Hz")
        plt.xlabel("Frequency MHz")
        plt.show()
        
    return setpoints_np.flatten(), data1d_np.flatten()

meta_info = []

for run_temp_id in run_ids_temp:
    metadata = get_metadata(run_temp_id)
    temp = metadata.get("Temperature")
    v_gate = metadata.get("qdac_ch06_dc_constant_V")
    
    if temp is None or v_gate is None:
        raise KeyError(f"Missing metadata in measurement {run_temp_id} (Temperature or Gate Voltage).")
    
    meta_info.append({
        "temperature": temp,
        "v_gate_ch06": v_gate
    })


plt.figure(figsize=(10, 6))
for run_id, meta in zip(run_ids, meta_info):
    temp = meta["temperature"]
    v_gate = meta["v_gate_ch06"]
    # set the detuning 
    detuning = v_gate+1.303
    freq, psd = extract_1d(run_id, data_1d_name="avg_avg_psd_nodrive", setpoint_name='freq_param', plot=False)
    label_str = f"T = {round(temp * 1000, 2)} mK, Δ  = {round(detuning * 1000, 2)} mV"
    plt.plot(freq, psd, label=label_str)

plt.xlabel("Frequency (MHz)")
plt.ylabel("Avg PSD (W/Hz)")
plt.title("PSD vs Frequency for Each Measurement")
plt.legend()
plt.tight_layout()
plt.show()


