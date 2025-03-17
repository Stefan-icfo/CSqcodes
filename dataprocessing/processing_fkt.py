import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
from scipy.ndimage import gaussian_filter



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





