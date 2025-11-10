import math
import matplotlib.pyplot as plt
import pandas as pd
import qcodes as qc
import numpy as np
from matplotlib.colors import LogNorm
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
# Database location
#qc.config["core"]["db_location"] = ".\Data\Raw_data\CD12_B5_F4v11"
#run_id=1946 #g2 compensated linesweep
#run_id=3 #g2 compensated linesweep
threshold = 150e-6
constant_slope=-200e-6
outer_gate_ch=2
avg_num=5
crosscap=-0.019
asymmetry=1  # Position between electron numbers, 1: no assymetry

qc.config["core"]["db_location"] = ".\Data\Raw_data\CD12_B5_F4v18_171025"
run_id=3 #g2 compensated linesweep
threshold = 150e-6
constant_slope=-200e-6
outer_gate_ch=2
avg_num=5
crosscap=-0.019
asymmetry=3  # Position between electron numbers, 1: no assymetry

#qc.config["core"]["db_location"] = '.\\Data\\Raw_data\\CD12_B5_F4v19_211025.db'
#run_id=35 #g2 compensated linesweep
#threshold = 200e-6
#constant_slope=0e-6
#outer_gate_ch=3
#avg_num=5
#crosscap=-0.012
#asymmetry=1  # Position between electron numbers, 1: no assymetry

#qc.config["core"]["db_location"] = r"C:\Users\sforstner\Desktop\Triton database\CD12_B5_F4v13.db"#for lever arm determinations
#run_id=98#all5g
#run_id=102#g1
#run_id=104#g2
#run_id=106#g3
#run_id=108#g4
#run_id=110#g5
#threshold = 150e-6
#constant_slope=-250e-6#for 5g together
#constant_slope=-100e-6#g1
#constant_slope=-70e-6#g2
#constant_slope=-50e-6#g3
#constant_slope=-40e-6#g4g5

#outer_gate_ch=2
#avg_num=9


#qc.config["core"]["db_location"] = r"C:\Users\sforstner\Desktop\Triton database\CD12_B5_F4v14.db"#g2 compensated linesweep






qc.config["core"]["db_location"] = ".\Data\Raw_data\CD12_B5_F4v18_171025.db"
print("tryna open db at "+qc.config["core"]["db_location"])
run_id=866#35 #g2 compensated linesweep
threshold = 150e-6
constant_slope=0e-6
outer_gate_ch=2
avg_num=5
crosscap=-0.018
asymmetry=1  # Position between electron numbers, 1: no assymetry

qc.config["core"]["db_location"] = ".\Data\Raw_data\CD12_B5_F4v18_171025.db"#repeat 1946 linesweep for x-th time on 311025
run_id=874 #g2 compensated linesweep
threshold = 150e-6
constant_slope=-200e-6
outer_gate_ch=2
avg_num=5
crosscap=-0.019
asymmetry=1  # Position between electron numbers, 1: no assymetry


#qc.config["core"]["db_location"] = '.\\Data\\Raw_data\\CD12_B5_F4v22_29_10_25.db'
#print("tryna open db at "+qc.config["core"]["db_location"])
#run_id=33#tensioned, g1=0.8
#run_id=35#tensioned, g1=0.84
#run_id=37#tensioned, g1=0.78
#threshold = 150e-6
#constant_slope=0e-6
#outer_gate_ch=2
#avg_num=5
#crosscap=-0.018
#asymmetry=1  # Position between electron numbers, 1: no assymetry
#qc.config["core"]["db_location"] = '.\\Data\\Raw_data\\CD12_B5_F4v19_211025.db'
#print("tryna open db at "+qc.config["core"]["db_location"])
#run_id=35#35 #g2 compensated linesweep
#threshold = 150e-6
#constant_slope=0e-6
#outer_gate_ch=2
#avg_num=5
#crosscap=-0.018
#asymmetry=1  # Position between electron numbers, 1: no assymetry


qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v25_02_11_25.db"
)

run_id=331 #g2 compensated linesweep
threshold = 150e-6
constant_slope=-200e-6
outer_gate_ch=2
avg_num=5
crosscap=-0.019
asymmetry=1  # Position between electron numbers, 1: no assymetry


qc.config["core"]["db_location"] = (
    "C:\\Users\\LAB-nanooptomechanic\\Documents\\MartaStefan\\CSqcodes\\Data\\Raw_data\\CD12_B5_F4v22_29_10_25.db"
)

run_id=223 #g2 compensated linesweep
threshold = 150e-6
constant_slope=-200e-6
outer_gate_ch=2
avg_num=5
crosscap=-0.019
asymmetry=1  # Position between electron numbers, 1: no assymetry




cs_gate_V, outer_gate_V, G_data=extract_2d(run_id,
               data_2d_name="G",
               setpoints1_name='QDAC_ch06_dc_constant_V',  # cs
               setpoints2_name=f'QDAC_ch0{outer_gate_ch}_dc_constant_V',  # gate 4
               plot=True, log=False, progress_report=False)

_, _, sens_data=extract_2d(run_id,
               data_2d_name="I_sens",
               setpoints1_name='QDAC_ch06_dc_constant_V',  # cs
               setpoints2_name=f'QDAC_ch0{outer_gate_ch}_dc_constant_V',  # gate 4
               plot=False, log=False, progress_report=False)

max_G_Vcs = cs_gate_V[np.argmax(G_data, axis=1)]



G_data_smoothed = np.apply_along_axis(
    lambda x: centered_moving_average(x, n=avg_num), 
    axis=1,  # Apply along cs_gate direction
    arr=G_data
)

max_G_Vcs_Gavg = cs_gate_V[np.argmax(G_data_smoothed, axis=1)]

left_max_sens_list, right_max_sens_list = [], []
for i in range(len(max_G_Vcs_Gavg)):
    Gmax_id = np.argmax(G_data_smoothed[i])

    left_slice  = sens_data[i, :Gmax_id]          # row i, left of peak
    right_slice = sens_data[i, Gmax_id+1:]        # row i, right of peak

    left_max  = np.max(left_slice)  if left_slice.size  else np.nan
    right_max = np.max(right_slice) if right_slice.size else np.nan

    left_max_sens_list.append(left_max)
    right_max_sens_list.append(right_max)




# Plot the cs_gate voltage of maximum G vs outer_gate voltage
plt.plot(outer_gate_V, max_G_Vcs_Gavg, '-',linewidth=1)
plt.xlabel(f'Outer Gate Voltage (V) - ch0{outer_gate_ch}')
plt.ylabel('CS Gate Voltage at Max avgG (V)')
plt.title('Position of Maximum Conductance')
plt.grid(True, alpha=0.3)
plt.show()

max_G_Vcs_Gavg_Virtual=max_G_Vcs_Gavg-outer_gate_V*crosscap

plt.plot(outer_gate_V, max_G_Vcs_Gavg_Virtual, '-',linewidth=1)
plt.xlabel(f'Outer Gate Voltage (V) - ch0{outer_gate_ch}')
plt.ylabel('corrected CS Gate Voltage at Max G (V)')
plt.title('corrected cs gate')
plt.grid(True, alpha=0.3)
plt.show()

derivative = np.diff(max_G_Vcs_Gavg_Virtual)  # This gives you n-1 elements
derivative = np.insert(derivative, 0, 0)  # Insert 0 at the beginning

plt.plot(outer_gate_V, derivative, '-',linewidth=1)
plt.xlabel(f'Outer Gate Voltage (V) - ch0{outer_gate_ch}')
plt.ylabel('CS Gate Voltage at Max G (V)')
plt.title('derivative')
plt.grid(True, alpha=0.3)
plt.show()


derivative_unsloped=derivative

plt.plot(outer_gate_V, derivative_unsloped, '-',linewidth=1)
plt.xlabel(f'Outer Gate Voltage (V) - ch0{outer_gate_ch}')
plt.ylabel('CS Gate Voltage at Max G (V)')
plt.title('derivative_unsloped')
plt.grid(True, alpha=0.3)
plt.show()



# Mask of values above threshold
above_threshold = abs(derivative_unsloped) > threshold

# Check if adjacent values are above threshold
adjacent_above = np.zeros_like(above_threshold, dtype=bool)
adjacent_above[1:] |= above_threshold[:-1]   # Previous neighbor
adjacent_above[:-1] |= above_threshold[1:]   # Next neighbor

# Keep value if it's above threshold OR has a neighbor above threshold
keep_mask = above_threshold | adjacent_above

# Set to zero where we don't want to keep
derivative_filtered = derivative_unsloped.copy()
derivative_filtered[~keep_mask] = 0

plt.plot(outer_gate_V, derivative_filtered, '-',linewidth=1)
plt.xlabel(f'Outer Gate Voltage (V) - ch0{outer_gate_ch}')
plt.ylabel('CS Gate Voltage at Max G (V)')
plt.title('derivative_filtered')
#plt.grid(True, alpha=0.3)
plt.show()

# Step 1: Cut off values where outer_gate_V < 0
derivative_positive_gate = derivative_filtered.copy()
derivative_positive_gate[outer_gate_V < 0] = 0

# Plot only the positive outer gate region
mask = outer_gate_V >= 0  # Changed from < to >=
plt.plot(outer_gate_V[mask], derivative_positive_gate[mask], '-', linewidth=1)
plt.xlabel(f'Outer Gate Voltage (V) - ch0{outer_gate_ch}')
plt.ylabel('Derivative (unsloped)')
plt.title('Derivative filtered (outer gate >= 0 only)')
#plt.grid(True, alpha=0.3)
plt.show()

# Step 2: Find and sum clusters
derivative_clustered = np.zeros_like(derivative_positive_gate)

# Find where non-zero values are
non_zero_mask = derivative_positive_gate != 0
non_zero_indices = np.where(non_zero_mask)[0]

print(f"Number of non-zero indices: {len(non_zero_indices)}")

if len(non_zero_indices) > 0:
    # Group consecutive indices into clusters
    clusters = []
    current_cluster = [non_zero_indices[0]]
    
    for i in range(1, len(non_zero_indices)):
        if non_zero_indices[i] == non_zero_indices[i-1] + 1:
            # Consecutive - add to current cluster
            current_cluster.append(non_zero_indices[i])
        else:
            # Gap - save current cluster and start new one
            clusters.append(current_cluster)
            current_cluster = [non_zero_indices[i]]
    
    # Don't forget the last cluster
    clusters.append(current_cluster)
    
    print(f"Number of clusters found: {len(clusters)}")
    
    # For each cluster, sum and place at peak
    for cluster_indices in clusters:
        cluster_values = derivative_positive_gate[cluster_indices]
        cluster_sum = np.sum(cluster_values)
        
        # Find maximum absolute value position
        max_pos = cluster_indices[np.argmax(np.abs(cluster_values))]
        
        derivative_clustered[max_pos] = cluster_sum
        
        print(f"Cluster at indices {cluster_indices[0]}-{cluster_indices[-1]}: sum={cluster_sum:.6f}, peak at index {max_pos}")

# Plot
plt.figure(figsize=(10, 6))
plt.stem(outer_gate_V, derivative_clustered, linefmt='C0-', markerfmt='C0o', basefmt='k-')
plt.xlabel(f'Outer Gate Voltage (V) - ch0{outer_gate_ch}')
plt.ylabel('Summed Derivative at Peak')
plt.title('Clustered Derivative (summed at peaks)')
#plt.grid(True, alpha=0.3)
plt.show()

# Extract non-zero values (cluster peaks) and their positions
cluster_positions = np.where(derivative_clustered != 0)[0]
cluster_heights = derivative_clustered[cluster_positions]
cluster_gate_voltages = outer_gate_V[cluster_positions]

# Enumerate clusters
cluster_numbers = np.arange(1, len(cluster_heights) + 1)

# Calculate distances between consecutive clusters
distances = np.diff(cluster_gate_voltages)
distances_full = np.insert(distances, 0, 0)

# Create figure with two y-axes
fig, ax1 = plt.subplots(figsize=(12, 6))

# Bar width and positions
bar_width = 0.35
x_pos = np.arange(len(cluster_numbers))

# First y-axis: Jump sizes (blue bars)
color1 = 'tab:blue'
ax1.set_xlabel('Electron Number')
ax1.set_ylabel('CS voltage jump size (mV)', color=color1)
bars1 = ax1.bar(x_pos - bar_width/2, cluster_heights * 1000, bar_width,
                label='CS voltage jump size (mV)', color=color1, alpha=0.8)
ax1.tick_params(axis='y', labelcolor=color1)
ax1.set_xticks(x_pos)
ax1.set_xticklabels(cluster_numbers)

# Second y-axis: Distances (red bars)
ax2 = ax1.twinx()
color2 = 'tab:red'
ax2.set_ylabel('Distance to previous jump (V)', color=color2)
bars2 = ax2.bar(x_pos + bar_width/2, distances_full, bar_width,
                label='Distance to previous jump', color=color2, alpha=0.8)
ax2.tick_params(axis='y', labelcolor=color2)

# Add legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.title(f'Filling Electrons - Electrostatic Linesweep Analysis - run_id {run_id} ')
#ax1.grid(True, alpha=0.3, axis='y')
fig.tight_layout()
plt.show()

# Print detailed info
for i in range(len(cluster_numbers)):
    if i == 0:
        print(f"Electron {cluster_numbers[i]}: jump size={cluster_heights[i]:.6f}, at V={cluster_gate_voltages[i]:.3f} V")
    else:
        print(f"Electron {cluster_numbers[i]}: jump size={cluster_heights[i]:.6f}, at V={cluster_gate_voltages[i]:.3f} V, distance from previous={distances[i-1]:.4f} V")


# Plot only distances with left y-axis
fig, ax = plt.subplots(figsize=(12, 6))

bar_width = 0.4
x_pos = np.arange(len(cluster_numbers))

# Plot distance bars
bars = ax.bar(x_pos, distances_full, bar_width, 
              color='tab:red', alpha=0.8, label='Distance to previous jump')

# Add gate voltage labels on top of each bar
for i, (x, height, voltage) in enumerate(zip(x_pos, distances_full, cluster_gate_voltages)):
    ax.text(x, height, f'{voltage:.2f}V', 
            ha='center', va='bottom', fontsize=8, rotation=0)

ax.set_xlabel('Electron Number')
ax.set_ylabel('Distance to previous jump (V)')
ax.set_title(f'Distance Between Electron Additions - run_id {run_id}')
ax.set_xticks(x_pos)
ax.set_xticklabels(cluster_numbers)
#ax.grid(True, alpha=0.3, axis='y')
ax.legend()

fig.tight_layout()
plt.show()

# Calculate midpoints between consecutive electron additions
midpoints = (cluster_gate_voltages[:-1] + asymmetry*cluster_gate_voltages[1:]) / (1+asymmetry)
midpoint_numbers = cluster_numbers[:-1] + 0.5  # Position between electron numbers

##find sensivitities and maxG around midpoints
idx = np.abs(outer_gate_V[:, None] - midpoints).argmin(axis=0)

idx = np.asarray(idx, dtype=int)
maxG_at_mid_nearest  = np.asarray(max_G_Vcs_Gavg)[idx]
leftsenspoints  = np.asarray(left_max_sens_list)[idx]
rightsenspoints = np.asarray(right_max_sens_list)[idx]


midpoints_rounded = [float(f'{x:.6g}') for x in midpoints]
print(f"midpoints = {midpoints_rounded}")
# Plot midpoints
fig, ax = plt.subplots(figsize=(12, 6))

ax.plot(midpoint_numbers, midpoints, 'o-', color='tab:green', 
        markersize=8, linewidth=2, label='Midpoint voltage')

# Optionally add the actual electron positions as reference
ax.plot(cluster_numbers, cluster_gate_voltages, 's', color='tab:blue', 
        markersize=6, alpha=0.5, label='Electron addition voltage')

ax.set_xlabel('Electron Number')
ax.set_ylabel('Outer Gate Voltage (V)')
ax.set_title(f'Midpoints Between Electron Additions - run_id {run_id}')
#ax.grid(True, alpha=0.3)
ax.legend()

fig.tight_layout()
plt.show()

plt.plot(cluster_numbers[:-1] + 0.5,maxG_at_mid_nearest)
plt.title("maxG")
plt.show()

plt.plot(cluster_numbers[:-1] + 0.5,leftsenspoints,'g*',label='left')
plt.plot(cluster_numbers[:-1] + 0.5,rightsenspoints,'r*',label='right')
plt.title("sensitivity")
plt.legend()
plt.show()


