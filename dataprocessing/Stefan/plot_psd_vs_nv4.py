#v4: 150 and 140M modes from two different databases

import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator


def analyze_fit_quality_l(x, y, popt):
    """
    This tells you how good your fit actually is
    """
    y_fit = lorentzian_fkt(x, *popt)
    residuals = y - y_fit
    
   
    

    return residuals


def analyze_fit_quality_g(x, y, popt):
    """
    This tells you how good your fit actually is
    """
    y_fit = gaussian_fkt(x, *popt)
    residuals = y - y_fit
    
   
    

    return residuals


#qc.config["core"]["db_location"] = r"C:\Users\sforstner\Desktop\Triton database\CD12_B5_F4v11.db"
qc.config["core"]["db_location"] ='.\\Data\\Raw_data\\CD12_B5_F4v11.db'
qc.config["core"]["db_location"] ='.\\Data\\Raw_data\\CD12_B5_F4v19_211025.db'#database for first dataset (150M)
print("Opening DB:", qc.config["core"]["db_location"])


e_nr=True
fit=True
Gamma_guess=2e3
background_id=654







run_ids1=[1171,1188,1205,1221,1238,1255,1271]#in db19_2110

electron_nrs1=[1,2,3,4,5,6,7]



run_ids2=[527,540,553,566,579,595,609,647,661,673,682,691,700,709,721,747,756]#in dbv1171025

electron_nrs2=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

run_ids2=[553,566,579,595,609,647,661,673,682,691,700,709,721,747,756]#in dbv1171025

electron_nrs2=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]




areas=[]
frequencies=[]
g2_voltages=[]
sensitivities=[]
 

areas1=[]
frequencies1=[]
g2_voltages1=[]
sensitivities1=[]
linewidths1l=[]
linewidths1g=[]
L_errors1=[]
G_errors1=[]

for run_id in run_ids1:#all
    metadata_temp=get_metadata(run_id-1,print_it=False,return_data=True)
    area=metadata_temp['integral_over_substracted_psd']
    g2_voltage=metadata_temp['qdac_ch02_dc_constant_V']
    frequency=metadata_temp['center_freq']
    sensitivity=metadata_temp['I_sens_sit']
    print(f"area {area:.4g}")
    areas1.append(area)
    g2_voltages1.append(g2_voltage)
    frequencies1.append(frequency)
    sensitivities1.append(sensitivity)
    if fit:
            freq_axis,spectrum=extract_1d(run_id, data_1d_name = "avg_avg_psd_nodrive_substracted", setpoint_name = 'freq_param',  plot = False,return_exp_name=False)
    
            

            initial_guess=[frequency,Gamma_guess,max(centered_moving_average(spectrum,n=5)),0]
            #freq_span=max(compressed_freq_array)-min(compressed_freq_array)
            # Define bounds for the parameters
            #lower_bounds = [min(compressed_freq_array)+0.25*freq_span, 0, max(avg_avg_psd_nodrive/1.5, 0]  # Replace with appropriate lower bounds
            #upper_bounds = [max(compressed_freq_array)-0.25*freq_span, 1e3, np.inf, np.inf]  # Gamma is bounded to <= 1e3
            #plt.plot(freq_axis,spectrum)
            #gaussian,g_area=gaussian_fkt_w_area(freq_axis,initial_guess[0],initial_guess[1],initial_guess[2],initial_guess[3])
            #plt.plot(freq_axis,gaussian)
            #plt.show()
            try:
                    popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, freq_axis, centered_moving_average(spectrum,n=5), p0=initial_guess)
                    lorentzian,l_area=lorentzian_fkt_w_area(freq_axis,popt[0],popt[1],popt[2],popt[3])
                   
                    lorentzian=lorentzian_fkt(freq_axis,popt[0],popt[1],popt[2],popt[3])

                    mask = np.abs(freq_axis - popt[0]) <= 3e3
                    residuals_l = analyze_fit_quality_l(freq_axis[mask], spectrum[mask], popt)
                    
                    
                    
                    area_error_l=sum(abs(residuals_l))
                    #area_error_l=np.sqrt(sum(abs(residuals_l)**2))
                    L_errors1.append(area_error_l)
                    gamma_l=abs(popt[1])
                    linewidths1l.append(gamma_l)

                    popt, pcov = scp.optimize.curve_fit(gaussian_fkt, freq_axis, centered_moving_average(spectrum,n=5), p0=initial_guess)
                    gaussian,g_area=gaussian_fkt_w_area(freq_axis,popt[0],popt[1],popt[2],popt[3])
                    residuals_g = analyze_fit_quality_g(freq_axis[mask], spectrum[mask], popt)
                    area_error_g=sum(abs(residuals_g))
                    G_errors1.append(area_error_g)
                    gamma_g=abs(popt[1])
                    linewidths1g.append(gamma_g)
                   # plt.plot(freq_axis,centered_moving_average(spectrum,n=5))
                   # plt.plot(freq_axis,lorentzian,"g")
                   # plt.plot(freq_axis[mask],residuals_l,"g")
                   # plt.plot(freq_axis,gaussian,"r-")
                   # plt.title(run_id)
                   # plt.show()
            except:
                    popt=initial_guess
                    print("fitting error")


    

qc.config["core"]["db_location"] ='.\\Data\\Raw_data\\CD12_B5_F4v18_171025.db'#database for second dataset (140M)
g2_voltages2=[]
areas2=[]
sensitivities2=[]
frequencies2=[]
linewidths2l=[]
linewidths2g=[]
L_errors2=[]
G_errors2=[]
for run_id in run_ids2:
        metadata_temp=get_metadata(run_id-1,print_it=False,return_data=True)
        area=metadata_temp['integral_over_substracted_psd']
        g2_voltage=metadata_temp['qdac_ch02_dc_constant_V']
        frequency=metadata_temp['center_freq']
        sensitivity=metadata_temp['I_sens_sit']
        print(f"area {area:.4g}")
        areas2.append(area)
        g2_voltages2.append(g2_voltage)
        frequencies2.append(frequency)
        sensitivities2.append(sensitivity)
        if fit:
            freq_axis,spectrum=extract_1d(run_id, data_1d_name = "avg_avg_psd_nodrive_substracted", setpoint_name = 'freq_param',  plot = False,return_exp_name=False)
    
            

            initial_guess=[frequency,Gamma_guess,max(centered_moving_average(spectrum,n=5)),0]
            #freq_span=max(compressed_freq_array)-min(compressed_freq_array)
            # Define bounds for the parameters
            #lower_bounds = [min(compressed_freq_array)+0.25*freq_span, 0, max(avg_avg_psd_nodrive/1.5, 0]  # Replace with appropriate lower bounds
            #upper_bounds = [max(compressed_freq_array)-0.25*freq_span, 1e3, np.inf, np.inf]  # Gamma is bounded to <= 1e3
            #plt.plot(freq_axis,spectrum)
            #gaussian,g_area=gaussian_fkt_w_area(freq_axis,initial_guess[0],initial_guess[1],initial_guess[2],initial_guess[3])
            #plt.plot(freq_axis,gaussian)
            #plt.show()
            try:
                    popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, freq_axis, centered_moving_average(spectrum,n=5), p0=initial_guess)
                    lorentzian,l_area=lorentzian_fkt_w_area(freq_axis,popt[0],popt[1],popt[2],popt[3])
                   
                    lorentzian=lorentzian_fkt(freq_axis,popt[0],popt[1],popt[2],popt[3])

                    mask = np.abs(freq_axis - popt[0]) <= 3e3
                    residuals_l = analyze_fit_quality_l(freq_axis[mask], spectrum[mask], popt)
                    
                    
                    
                    area_error_l=sum(abs(residuals_l))
                    #area_error_l=np.sqrt(sum(abs(residuals_l)**2))
                    L_errors2.append(area_error_l)
                    gamma_l=abs(popt[1])
                    linewidths2l.append(gamma_l)

                    popt, pcov = scp.optimize.curve_fit(gaussian_fkt, freq_axis, centered_moving_average(spectrum,n=5), p0=initial_guess)
                    gaussian,g_area=gaussian_fkt_w_area(freq_axis,popt[0],popt[1],popt[2],popt[3])
                    residuals_g = analyze_fit_quality_g(freq_axis[mask], spectrum[mask], popt)
                    area_error_g=sum(abs(residuals_g))
                    G_errors2.append(area_error_g)
                    gamma_g=abs(popt[1])
                    linewidths2g.append(gamma_g)
                    #plt.plot(freq_axis,centered_moving_average(spectrum,n=5))
                    #plt.plot(freq_axis,lorentzian,"g")
                    #plt.plot(freq_axis[mask],residuals_l,"g")
                    #plt.plot(freq_axis,gaussian,"r-")
                    #plt.title(run_id)
                    #plt.show()
            except:
                    popt=initial_guess
                    print("fitting error")
            

    
 

plt.plot(g2_voltages1,areas1,'g*')
plt.plot(g2_voltages2,areas2,'r*')
plt.title("areas_vs_g2_voltages with best sensitivity ref linesweep 1946 (repeated)",fontsize=14)
plt.xlabel("g2voltage")
plt.ylabel("PSD area")
for i, (x, y,run_id) in enumerate(zip(g2_voltages1, areas1,run_ids1)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
for i, (x, y,run_id) in enumerate(zip(g2_voltages2, areas2,run_ids2)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
plt.show()

if e_nr==True:
    #plt.plot(electron_nrs1,areas1,'g*')
    plt.errorbar(electron_nrs2, areas2, yerr=L_errors2, fmt='g*', capsize=5, label='Areas with Lerrors')
    plt.errorbar(electron_nrs1, areas1, yerr=L_errors1, fmt='g*', capsize=5, label='Areas with Lerrors')
    #plt.plot(electron_nrs2,L_errors,'r*')
    plt.title("areas_vs_e_nr with best sensitivity ref linesweep 1946", fontsize=14)
    plt.xlabel("nr electrons")
    plt.ylabel("PSD area")
    plt.xlim(left=0)
    plt.ylim(bottom=0) 
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for i, (x, y,run_id) in enumerate(zip(electron_nrs2, areas2,run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
    plt.show()

plt.plot(g2_voltages1,frequencies1,'g*')
plt.plot(g2_voltages2,frequencies2,'r*')
plt.title("frequencies_vs_e_nr with best sensitivity ref linesweep 1946 (repeated)")
plt.xlabel("g2voltage")
plt.ylabel("frequency")

for i, (x, y,run_id) in enumerate(zip(g2_voltages2, frequencies2,run_ids2)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
for i, (x, y,run_id) in enumerate(zip(g2_voltages1, frequencies1,run_ids1)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
plt.show()

if e_nr==True:
    plt.plot(electron_nrs1,frequencies1,'g*')
    plt.plot(electron_nrs2,frequencies2,'r*')
    plt.title("frequencies_vs_e_nr with best sensitivity ref linesweep 1946(repeated)")
    plt.xlabel("nr electrons")
    plt.ylabel("frequency")
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for i, (x, y,run_id) in enumerate(zip(electron_nrs2, frequencies,run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    for i, (x, y,run_id) in enumerate(zip(electron_nrs1, frequencies,run_ids1)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
    plt.show()
 
if e_nr==True:
    plt.plot(electron_nrs1,sensitivities1,'g*')
    plt.plot(electron_nrs2,sensitivities2,'r*')
    plt.title("sensitivities")
    plt.xlabel("nr electrons")
    plt.ylabel("sensitivities")
    plt.xlim(left=0)
    plt.ylim(bottom=0) 
    for i, (x, y,run_id2) in enumerate(zip(electron_nrs2, areas2,run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    for i, (x, y,run_id2) in enumerate(zip(electron_nrs1, areas1,run_ids1)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
 
    plt.show()


    #plt.plot(electron_nrs1, np.sqrt(areas1)*np.array(sensitivities1), 'g*')
    #plt.plot(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), 'r*')




    
    #################ERRORBARS##################

    errors_scaled2 = np.array(sensitivities2) * np.array(G_errors2) / (2 * np.sqrt(areas2))
    errors_scaled1 = np.array(sensitivities1) * np.array(G_errors1) / (2 * np.sqrt(areas1))

    plt.errorbar(electron_nrs1, np.sqrt(areas1)*np.array(sensitivities1), yerr=errors_scaled1,fmt='g*', label='Dataset 1')
    plt.errorbar(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), yerr=errors_scaled2, fmt='r*', capsize=5, label='Dataset 2')










    plt.title("scaled sqrt areas, both gaussian error")
    plt.xlabel("nr electrons")
    plt.ylabel("scaled psd")
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    # Calculate linear regression through origin first dataset
    x1 = np.array(electron_nrs1)
    y1 = np.sqrt(areas1) * np.array(sensitivities1)
    # For regression through origin: slope = sum(x*y) / sum(x^2)
    slope = np.sum(x1 * y1) / np.sum(x1**2)
    # Create regression line
    x1_line = np.linspace(0, max(electron_nrs1), 100)
    y1_line = slope * x1_line
    # Plot the regression line
    plt.plot(x1_line, y1_line, 'g-', label=f'y = {slope:.4e}x', linewidth=2)
    # Add labels
    for i, (x, y, run_id) in enumerate(zip(electron_nrs1, np.sqrt(areas1)*np.array(sensitivities1), run_ids1)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    
    # Calculate linear regression through origin second dataset
    x2 = np.array(electron_nrs2)
    y2 = np.sqrt(areas2) * np.array(sensitivities2)
    # For regression through origin: slope = sum(x*y) / sum(x^2)
    slope = np.sum(x2 * y2) / np.sum(x2**2)
    # Create regression line
    x2_line = np.linspace(0, max(electron_nrs2), 100)
    y2_line = slope * x2_line
    # Plot the regression line
    plt.plot(x2_line, y2_line, 'r-', label=f'y = {slope:.4e}x', linewidth=2)
    # Add labels
    for i, (x, y, run_id) in enumerate(zip(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')

    plt.legend()
    plt.show()


    plt.errorbar(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), yerr=errors_scaled2, fmt='g*', capsize=5, label='Dataset 2')
    plt.errorbar(electron_nrs1, np.sqrt(areas1)*np.array(sensitivities1), yerr=errors_scaled1, fmt='g*', capsize=5, label='Dataset 1')

    plt.title("scaled sqrt areas w Lortenzian error")
    plt.xlabel("nr electrons")
    plt.ylabel("scaled psd")
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    
    # Calculate linear regression through origin first dataset
    x2 = np.array(electron_nrs2)
    y2 = np.sqrt(areas2) * np.array(sensitivities2)
    # For regression through origin: slope = sum(x*y) / sum(x^2)
    slope = np.sum(x2 * y2) / np.sum(x2**2)
    # Create regression line
    x2_line = np.linspace(0, max(electron_nrs2), 100)
    y2_line = slope * x2_line
    # Plot the regression line
    plt.plot(x2_line, y2_line, 'r-', label=f'y = {slope:.4e}x', linewidth=2)
    # Add labels
    for i, (x, y, run_id) in enumerate(zip(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')

    plt.legend()
    plt.show()



    errors_scaledG = np.array(sensitivities2) * np.array(G_errors2) / (2 * np.sqrt(areas2))

    plt.errorbar(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), yerr=errors_scaledG, fmt='r*', capsize=5, label='Dataset 2')

    plt.title("scaled sqrt areas w Gaussian error")
    plt.xlabel("nr electrons")
    plt.ylabel("scaled psd")
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    
    # Calculate linear regression through origin first dataset
    x2 = np.array(electron_nrs2)
    y2 = np.sqrt(areas2) * np.array(sensitivities2)
    # For regression through origin: slope = sum(x*y) / sum(x^2)
    slope = np.sum(x2 * y2) / np.sum(x2**2)
    # Create regression line
    x2_line = np.linspace(0, max(electron_nrs2), 100)
    y2_line = slope * x2_line
    # Plot the regression line
    plt.plot(x2_line, y2_line, 'r-', label=f'y = {slope:.4e}x', linewidth=2)
    # Add labels
    for i, (x, y, run_id) in enumerate(zip(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')

    plt.legend()
    plt.show()


    
    
    idx = run_ids2.index(661)
    linewidths2l.pop(idx)
    linewidths2g.pop(idx)
    electron_nrs2.pop(idx)


    #idx = run_ids2.index(540)
    #linewidths2.pop(idx)
    #electron_nrs2.pop(idx)

    plt.plot(electron_nrs2, linewidths2l, 'r*')
    plt.plot(electron_nrs1, linewidths1l, 'g*')

    plt.title("linewidths Lorentzian fit")
    plt.xlabel("nr electrons")
    plt.ylabel("linewidths")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.show()


    plt.plot(electron_nrs2, linewidths2g, 'r*')
    plt.plot(electron_nrs1, linewidths1g, 'g*')

    plt.title("linewidths Gaussian fit")
    plt.xlabel("nr electrons")
    plt.ylabel("linewidths")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.show()

  

  