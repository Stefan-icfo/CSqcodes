import numpy as np
import matplotlib.pyplot as plt
import qcodes as qc
import os
qc.config["core"]["db_location"]="C:"+"\\"+"Users"+"\\"+"LAB-nanooptomechanic"+"\\"+"Documents"+"\\"+"MartaStefan"+"\\"+"CSqcodes"+"\\"+"Data"+"\\"+"Raw_data"+"\\"+'CD11_D7_C1_part2.db'
def extract_1d(run_id, data_1d_name="v_r", setpoint_name="zurich_oscs_freq", plot=True):
    experiments = qc.experiments()
    dataset = qc.load_by_id(run_id)

    interdeps = dataset.description.interdeps
    data_x = dataset.get_parameter_data(data_1d_name)
    setpoints_raw = data_x[data_1d_name][setpoint_name]
    setpoints_np = np.array(setpoints_raw)

    pdf_temp = dataset.to_pandas_dataframe_dict()
    data1d_raw = pdf_temp[data_1d_name]
    data1d_np = np.array(data1d_raw)
    
    if plot:
        plt.plot(setpoints_np, data1d_np)
        plt.title(f"Measurement {run_id}")
        plt.xlabel(setpoint_name)
        plt.ylabel(data_1d_name)
        plt.grid()
        plt.show()

    return setpoints_np.flatten(), data1d_np.flatten()

def save_extracted_data(run_id, output_folder, data_1d_name="v_r", setpoint_name="zurich_oscs_freq"):
    # Extract the data
    time, data = extract_1d(run_id, data_1d_name=data_1d_name, setpoint_name=setpoint_name, plot=False)

    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)
    
    # Save the data as a text file
    output_file_path = os.path.join(output_folder, f"run_{run_id}.txt")
    with open(output_file_path, 'w') as f:
        f.write(f"# {setpoint_name}, {data_1d_name}\n")  # Add header
        for t, d in zip(time, data):
            f.write(f"{t:.10f}, {d:.10f}\n")

    print(f"Data extracted and saved to '{output_file_path}'")

# Main function
def main():
    run_id = 1134 # ID of the measurement
    output_folder = r"\\files\groups\NanoOptoMechanics\Users\Marta\autoccrelationforRoger"
    save_extracted_data(run_id, output_folder)

if __name__ == "__main__":
    main()
