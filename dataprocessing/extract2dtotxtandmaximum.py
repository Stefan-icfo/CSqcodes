import math
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
import numpy as np
import os

# =======================
# Configurazione percorso
# =======================

run_id = 508  # <-- Cambia qui se vuoi analizzare altri run

dB_location = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part2.db"
output_folder = "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD11_D7_C1_part2"

# ========================
# Funzione principale
# ========================

def extract_2d(run_id, 
               data_2d_name="G",
               setpoints1_name='delta',
               setpoints2_name='gateV',
               plot=True,
               save_txt=True,
               dB_location=dB_location):

    print(f"Running extract_2d for run {run_id}...")

    qc.config["core"]["db_location"] = dB_location
    dataset = qc.load_by_id(run_id)

    # Dati grezzi
    pdf_temp = dataset.to_pandas_dataframe_dict()
    data2d_raw = pdf_temp[data_2d_name]
    data2d_np = np.array(data2d_raw)

    # Interdipendenze per gli assi
    interdeps = dataset.description.interdeps
    param_spec = interdeps.non_dependencies[0]  
    param_name = param_spec.name
    data_xy = dataset.get_parameter_data(param_spec)

    # Estrazione setpoint
    setpoints1_np = np.array(data_xy[param_name][setpoints1_name])
    setpoints2_np = np.array(data_xy[param_name][setpoints2_name])

    setpoints1 = np.unique(setpoints1_np)
    setpoints2 = np.unique(setpoints2_np)

    # Ricostruzione mappa 2D
    data_2d = np.zeros((len(setpoints1), len(setpoints2)))
    for i in range(len(data2d_np)): 
        row = np.where(setpoints1 == setpoints1_np[i])[0][0]  
        col = np.where(setpoints2 == setpoints2_np[i])[0][0]  
        data_2d[row, col] = data2d_np[i]

    # Plot
    if plot:
        plt.figure(figsize=(6, 5))
        plt.pcolor(setpoints2, setpoints1, data_2d, shading='auto', cmap='plasma')
        plt.title(f"Measurement {run_id}")
        plt.colorbar(label=data_2d_name)
        plt.xlabel('Gate voltage (V)')
        plt.ylabel('Delta (mV)')
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, f"measurement_{run_id}.png"), dpi=300, bbox_inches='tight')
        plt.show()

    # ============================
    # Salva file completo
    # ============================

    if save_txt:
        data_txt_path = os.path.join(output_folder, f"CD11_D7_C1_run{run_id}.txt")
        with open(data_txt_path, 'w', encoding='utf-8') as f:
            f.write("gateV(V)\tdelta(mV)\tG(μS)\n")
            for i in range(len(setpoints1)):
                for j in range(len(setpoints2)):
                    f.write(f"{setpoints2[j]:.8f}\t{setpoints1[i]:.8f}\t{data_2d[i, j]:.8f}\n")
        print(f"✅ File dati completo salvato in: {data_txt_path}")

    # ============================
    # Max G per ogni delta (valore)
    # ============================

    delta_vs_peakgate = []
    for i, delta in enumerate(setpoints1):
        g_row = data_2d[i, :]
        max_idx = np.argmax(g_row)
        gate_max = setpoints2[max_idx]
        delta_vs_peakgate.append((delta, gate_max))

    peak_txt_path = os.path.join(output_folder, f"maxpeakvsdetuning_{run_id}.txt")
    with open(peak_txt_path, 'w', encoding='utf-8') as f:
        f.write("delta(mV)\tgateV_at_max_G(V)\n")
        for delta, gate_peak in delta_vs_peakgate:
            f.write(f"{delta:.8f}\t{gate_peak:.8f}\n")
    print(f"✅ File picchi salvato in: {peak_txt_path}")

    # ============================
    # Max G per ogni delta (indice)
    # ============================

    delta_vs_index = []
    for i, delta in enumerate(setpoints1):
        g_row = data_2d[i, :]
        max_idx = int(np.argmax(g_row))  # posizione intera
        delta_vs_index.append((delta, max_idx))

    index_txt_path = os.path.join(output_folder, f"massimovsdetuning_{run_id}.txt")
    with open(index_txt_path, 'w', encoding='utf-8') as f:
        f.write("delta(mV)\tposizione_max\n")
        for delta, idx in delta_vs_index:
            f.write(f"{delta:.8f}\t{idx}\n")
    print(f"✅ File posizione massima salvato in: {index_txt_path}")

    return setpoints1, setpoints2, data_2d

# ============================
# Esegui quando chiamato
# ============================

if __name__ == "__main__":
    extract_2d(run_id=run_id)


