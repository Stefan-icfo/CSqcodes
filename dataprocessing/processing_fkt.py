import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
from scipy.ndimage import gaussian_filter
from database import *
import qcodes as qc
import time


import sqlite3

from qcodes.dataset.data_set import load_by_id

def rename(run_id, new_name, db_path=None):
    if db_path==None:
        db_path=qc.config["core"]["db_location"]

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    try:
        cursor.execute("UPDATE runs SET name = ? WHERE run_id = ?", (new_name, run_id))
        conn.commit()
        print(f"✅ Renamed run {run_id} to '{new_name}'")
    except Exception as e:
        print(f"❌ Failed to rename run {run_id}: {e}")
    finally:
        conn.close()


def compact_database_into_new_file(original_path, output_path):
    import sqlite3
    print(f"🧼 Compacting {original_path} → {output_path}...")
    conn = sqlite3.connect(original_path)
    try:
        conn.execute(f"VACUUM INTO '{output_path}';")
        print("✅ VACUUM INTO successful — new file is compacted.")
    except Exception as e:
        print(f"❌ VACUUM INTO failed: {e}")
    finally:
        conn.close()


def vacuum_database(db_path):
    print(f"\n🧼 Running VACUUM on {db_path}...")
    conn = sqlite3.connect(db_path)
    try:
        conn.execute("VACUUM;")
        print("✅ VACUUM complete — database file compacted.")
    except sqlite3.OperationalError as e:
        print(f"❌ VACUUM failed: {e}")
    finally:
        conn.close()


def delete_runs_by_id_list(run_ids, db_path):
    """
    Delete multiple QCoDeS runs from a database by run_id, including their results-* tables.

    Args:
        run_ids (list[int]): List of run_id integers to delete.
        db_path (str): Full path to the QCoDeS SQLite database file.
    """
    print(f"🧹 Connecting to {db_path}")
    time.sleep(2)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Get all table names in advance
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    all_tables = [row[0] for row in cursor.fetchall()]

    deleted_tables = 0

    for run_id in run_ids:
        print(f"\n🗑 Deleting run_id {run_id}")
        time.sleep(0.2)

        # Delete from common tables
        for table in all_tables:
            try:
                cursor.execute(f'PRAGMA table_info("{table}")')
                columns = [row[1] for row in cursor.fetchall()]
                if "run_id" in columns:
                    cursor.execute(f'DELETE FROM "{table}" WHERE run_id = ?', (run_id,))
                    print(f"  ✔ Deleted from {table}")
            except Exception as e:
                print(f"  ⚠ Could not clean {table}: {e}")

        # Drop any results-<run_id>-* tables
        for table in all_tables:
            if table.startswith(f"results-{run_id}-"):
                try:
                    cursor.execute(f'DROP TABLE "{table}"')
                    print(f"  🗑 Dropped table {table}")
                    deleted_tables += 1
                except Exception as e:
                    print(f"  ⚠ Could not drop table {table}: {e}")

    conn.commit()
    conn.close()
    print(f"\n✅ Deleted runs and dropped {deleted_tables} results-* tables.")




def Charging_lines(gate2, gate4, data_2d, run_id=0, plot=False): 

    # Threshold value was calculated by taking the mean + constant * standard deviation of all data points
    # This is to generalize the recognition of the charging lines to all possible charge stability diagrams in the database  
    threshold_value = np.mean(data_2d) + 3.8 * np.std(data_2d)


    coordinates = peak_local_max(
        data_2d,
        min_distance=2,              
        threshold_abs=threshold_value,
        exclude_border=False
    )

    print(f"Found {len(coordinates)} Coloumb peaks in the data.")


    peak_points = [(gate2[row], gate4[col]) for (row, col) in coordinates]

    if peak_points:
        g2_vals = [p[0] for p in peak_points]
        g4_vals = [p[1] for p in peak_points]
        print(f"\nGate 2 range of peaks: {min(g2_vals):.4f} V to {max(g2_vals):.4f} V")
        print(f"Gate 4 range of peaks: {min(g4_vals):.4f} V to {max(g4_vals):.4f} V")


    else:
        print("No peaks found")

    if plot: 
        plt.figure(figsize=(6, 5))
        plt.pcolor(gate4, gate2, data_2d)
        plt.title(f"Measurement {run_id} with Charging Peaks")
        plt.colorbar(label='signal_shift_Vxn_deriv')
        plt.xlabel('Gate 4 (V)')
        plt.ylabel('Gate 2 (V)')


        peak_g4 = [p[1] for p in peak_points]  # gate 4 is x axis
        peak_g2 = [p[0] for p in peak_points]  # gate 2 is y-axis
        plt.scatter(peak_g4, peak_g2, color='red', marker='x')

        plt.savefig(f"measurement_{run_id}_peaks.png", dpi=300, bbox_inches='tight')
        plt.show()

    # Outputs the ranges of gate 2 and gate 4 values for which the the charging lines are recognized
    return min(g2_vals), max(g2_vals), min(g4_vals), max(g4_vals)



def ICT_points(gate2, gate4, data_2d, run_id=0, threshold_std=4.2, plot=False):   
    
    # Apply Gaussian filter to reduce noise
    smoothed_data = gaussian_filter(data_2d, sigma=1) 
    # Arbiitary constant 
    threshold_value = np.mean(smoothed_data) - threshold_std * np.std(smoothed_data)

    # we are looking for the local minimum, hence we invert the data and threshold value 
    coordinates1 = peak_local_max(
        -smoothed_data,
        min_distance=0,              
        threshold_abs=-threshold_value,
        exclude_border=False
    )

    print(f"Found {len(coordinates1)} local minimum in the data.")


    min_points = [(gate2[row], gate4[col]) for (row, col) in coordinates1]

    
    if len(min_points) > 1:
        g2_vals1 = [p[0] for p in min_points]
        g4_vals1 = [p[1] for p in min_points]
        print(f"\nGate 2 range: {min(g2_vals1):.4f} V to {max(g2_vals1):.4f} V")
        print(f"Gate 4 range: {min(g4_vals1):.4f} V to {max(g4_vals1):.4f} V")
    
    else:
        print("No ICT found.")

    if plot: 
        min_g4 = [p[1] for p in min_points]  # gate 4 is x axis
        min_g2 = [p[0] for p in min_points]  # gate 2 is y-axis

        plt.figure(figsize=(6, 5))
        plt.pcolor(gate4, gate2, data_2d)
        plt.title(f"Measurement {run_id} with ICT length definied")
        plt.colorbar(label='signal_shift_Vxn_deriv')
        plt.xlabel('Gate 4 (V)')
        plt.ylabel('Gate 2 (V)')


        plt.scatter(min(min_g4), min(min_g2), color='yellow', marker='x')
        plt.scatter(max(min_g4), max(min_g2), color='yellow', marker='x')

        plt.show()
    # Possible to return a more relevant calculation
    return min(min_g2), max(min_g2), min(min_g4), max(min_g4)

def ICT_width(gate2, gate4, data_2d,threshold_std=4.2):   
    
    # Apply Gaussian filter to reduce noise
    smoothed_data = gaussian_filter(data_2d, sigma=1)


    threshold_value = np.mean(smoothed_data) - threshold_std * np.std(smoothed_data)

    # we are looking for the local minimum, hence we invert the data 
    coordinates1 = peak_local_max(
        -smoothed_data,
        min_distance=0,              
        threshold_abs=-threshold_value,
        exclude_border=False
    )

    print(f"Found {len(coordinates1)} local minimum in the data.")


    min_points = [(gate2[row], gate4[col]) for (row, col) in coordinates1]

    if len(min_points) < 2:
        print("Not enough points to measure width.")
    else:
        
        x_vals = [p[1] for p in min_points]
        y_vals = [p[0] for p in min_points]
        
        # Fit a line y = m*x + b using np.polyfit
        m, b = np.polyfit(x_vals, y_vals, 1)
        
        # Compute perpendicular distances from each point to that line
        distances = []
        for (y_i, x_i) in min_points:
            dist = abs(m * x_i - y_i + b) / np.sqrt(m**2 + 1)
            distances.append(dist)
        
    
        width = np.mean(distances)

    
        #print(f"Width (range of distances) = {width_range:.4e}")
        print(f"Mean Distances = {width:.4e}")
        # to be decided which calculation is more accurate 
        return width





