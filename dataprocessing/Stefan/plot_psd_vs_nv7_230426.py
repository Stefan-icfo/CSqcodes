#v6: 141125; tensioned and untensioned 140 M; for now have not yet removed the 150M points

import matplotlib.pyplot as plt
from dataprocessing.extract_fkts import *
from utils.CS_utils import *
import qcodes as qc
from database import *
import re
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
import copy
import pickle
import hashlib
import os


def analyze_fit_quality_l(x, y, popt):
    y_fit = lorentzian_fkt(x, *popt)
    residuals = y - y_fit
    return residuals


def analyze_fit_quality_g(x, y, popt):
    y_fit = gaussian_fkt(x, *popt)
    residuals = y - y_fit
    return residuals


def psd_area_to_vrms(area, R=50):
    return np.sqrt(area * R)  # PSD is in W/Hz, area in W → V_rms = sqrt(P*R)


def make_cache_key(run_ids, fit):
    key = str(run_ids) + str(fit)
    return hashlib.md5(key.encode()).hexdigest()[:8]


qc.config["core"]["db_location"] = r"D:\databases CD12_B5_F4\CD12_B5_F4v19_211025.db"
print("Opening DB:", qc.config["core"]["db_location"])


e_nr = True
fit = True
Gamma_guess = 2e3
background_id = 654


#manual_spectra_sum=np.array([4.3e-16,1.1e-15,1.3e-15,2.11e-15,5.7e-15])/170*190#these are the values from dBv27
#e_nrs_manual=[1,1,2,3,4]
#manual_sens=[0.73e-11,0.73e-11,0.99e-11,1.17e-11,1.29e-11]

manual_spectra_sum = [4.7e-16, 1.3e-15, 2.11e-15, 5.7e-15]  # dBv30, tensioned configuration
e_nrs_manual = [1, 2, 3, 4]
manual_sens = [1.1e-11, 0.99e-11, 1.17e-11, 1.29e-11]
manual_frequencies = [142.830e6, 142.924e6, 142.616e6, 142.0665e6]  # midpoint of calculation window, uncertainty 1-2kHz
g2V_manual = []  # to be filled
manual_errors = [2*3.4e-16, 2*6.6e-16, 7.7e-16, 2*2.1e-17]


run_ids1 = [1171, 1188, 1205, 1221, 1238, 1255, 1271]  # in db19_2110
electron_nrs1 = [1, 2, 3, 4, 5, 6, 7]

#run_ids2=[527,540,553,566,579,595,609,647,661,673,682,691,700,709,721,747,756]#in dbv18 71025
#electron_nrs2=[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

#run_ids2=[540,553,566,579,595,609,647,661,673,682,691,700,709,721,747,756]#in dbv18 71025
#electron_nrs2=[5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

run_ids2 = [553, 566, 579, 595, 609, 647, 661, 673, 682, 691, 700, 709, 721, 747, 756]  # in dbv18 171025
electron_nrs2 = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]


n_ph_th = 26
background_noise_psd = 5.5e-15


# ─── Dataset 1 ───────────────────────────────────────────────────────────────

cache_file1 = f"cache_ds1_{make_cache_key(run_ids1, fit)}.pkl"

if os.path.exists(cache_file1):
    print(f"Loading dataset 1 from cache ({cache_file1})...")
    with open(cache_file1, 'rb') as f:
        d = pickle.load(f)
    areas1         = d['areas1']
    g2_voltages1   = d['g2_voltages1']
    frequencies1   = d['frequencies1']
    sensitivities1 = d['sensitivities1']
    n_imp1         = d['n_imp1']
    linewidths1l   = d['linewidths1l']
    linewidths1g   = d['linewidths1g']
    L_errors1      = d['L_errors1']
    G_errors1      = d['G_errors1']
    vrms_fft1      = d['vrms_fft1']

else:
    areas1 = []
    frequencies1 = []
    g2_voltages1 = []
    sensitivities1 = []
    linewidths1l = []
    linewidths1g = []
    L_errors1 = []
    G_errors1 = []
    n_imp1 = []
    vrms_fft1 = []

    for run_id, e_nr_ in zip(run_ids1, electron_nrs1):
        metadata_temp = get_metadata(run_id-1, print_it=False, return_data=True)
        area = metadata_temp['integral_over_substracted_psd']
        g2_voltage = metadata_temp['qdac_ch02_dc_constant_V']
        frequency = metadata_temp['center_freq']
        sensitivity = metadata_temp['I_sens_sit']
        max_psd = metadata_temp['max_avg_avg_psd_']
        SNR_temp = (max_psd - background_noise_psd) / background_noise_psd
        print(f"area {area:.4g}")
        areas1.append(area)
        g2_voltages1.append(g2_voltage)
        frequencies1.append(frequency)
        sensitivities1.append(sensitivity)
        n_imp1.append(n_ph_th / SNR_temp)

        freq_axis_fft, v_fft = extract_1d(run_id, data_1d_name="V_fft_avg_avg", setpoint_name='freq_param', plot=False, return_exp_name=False)
        bg_mask_fft = np.abs(freq_axis_fft - frequency) > 10e3
        v_fft_sub = v_fft - np.mean(v_fft[bg_mask_fft])
        signal_mask_fft = np.abs(freq_axis_fft - frequency) <= 3e3
        vrms_fft1.append(np.sqrt(np.sum(v_fft_sub[signal_mask_fft]**2)))

        if fit:
            freq_axis, spectrum = extract_1d(run_id, data_1d_name="avg_avg_psd_nodrive_substracted", setpoint_name='freq_param', plot=False, return_exp_name=False)
            freq_axis, spectrum_unsubstracted = extract_1d(run_id, data_1d_name="avg_avg_psd_nodrive", setpoint_name='freq_param', plot=False, return_exp_name=False)

            initial_guess = [frequency, Gamma_guess, max(centered_moving_average(spectrum, n=5)), 0]
            try:
                popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, freq_axis, centered_moving_average(spectrum, n=5), p0=initial_guess)
                lorentzian, l_area = lorentzian_fkt_w_area(freq_axis, popt[0], popt[1], popt[2], popt[3])
                lorentzian = lorentzian_fkt(freq_axis, popt[0], popt[1], popt[2], popt[3])
                mask = np.abs(freq_axis - popt[0]) <= 3e3
                residuals_l = analyze_fit_quality_l(freq_axis[mask], spectrum[mask], popt)
                area_error_l = sum(abs(residuals_l))
                L_errors1.append(area_error_l)
                linewidths1l.append(abs(popt[1]))

                popt, pcov = scp.optimize.curve_fit(gaussian_fkt, freq_axis, centered_moving_average(spectrum, n=5), p0=initial_guess)
                gaussian, g_area = gaussian_fkt_w_area(freq_axis, popt[0], popt[1], popt[2], popt[3])
                residuals_g = analyze_fit_quality_g(freq_axis[mask], spectrum[mask], popt)
                area_error_g = sum(abs(residuals_g))
                G_errors1.append(area_error_g)
                linewidths1g.append(abs(popt[1]))
            except:
                popt = initial_guess
                print("fitting error")

    with open(cache_file1, 'wb') as f:
        pickle.dump({'areas1': areas1, 'g2_voltages1': g2_voltages1, 'frequencies1': frequencies1,
                     'sensitivities1': sensitivities1, 'n_imp1': n_imp1, 'linewidths1l': linewidths1l,
                     'linewidths1g': linewidths1g, 'L_errors1': L_errors1, 'G_errors1': G_errors1,
                     'vrms_fft1': vrms_fft1}, f)
    print(f"Dataset 1 cached to {cache_file1}")


# ─── Dataset 2 ───────────────────────────────────────────────────────────────

cache_file2 = f"cache_ds2_{make_cache_key(run_ids2, fit)}.pkl"

if os.path.exists(cache_file2):
    print(f"Loading dataset 2 from cache ({cache_file2})...")
    with open(cache_file2, 'rb') as f:
        d = pickle.load(f)
    areas2         = d['areas2']
    g2_voltages2   = d['g2_voltages2']
    frequencies2   = d['frequencies2']
    sensitivities2 = d['sensitivities2']
    n_imp2         = d['n_imp2']
    linewidths2l   = d['linewidths2l']
    linewidths2g   = d['linewidths2g']
    L_errors2      = d['L_errors2']
    G_errors2      = d['G_errors2']
    vrms_fft2      = d['vrms_fft2']

else:
    qc.config["core"]["db_location"] = r"C:\Users\sforstner\Desktop\Triton database\CD12_B5_F4v18_171025.db"  # 140MHz data

    g2_voltages2 = []
    areas2 = []
    sensitivities2 = []
    frequencies2 = []
    linewidths2l = []
    linewidths2g = []
    L_errors2 = []
    G_errors2 = []
    n_imp2 = []
    vrms_fft2 = []

    for run_id, e_nr_ in zip(run_ids2, electron_nrs2):
        metadata_temp = get_metadata(run_id-1, print_it=False, return_data=True)
        area = metadata_temp['integral_over_substracted_psd']
        g2_voltage = metadata_temp['qdac_ch02_dc_constant_V']
        frequency = metadata_temp['center_freq']
        sensitivity = metadata_temp['I_sens_sit']
        max_psd = metadata_temp['max_avg_avg_psd_']
        SNR_temp = (max_psd - background_noise_psd) / background_noise_psd
        print(f"area {area:.4g}")
        areas2.append(area)
        g2_voltages2.append(g2_voltage)
        frequencies2.append(frequency)
        sensitivities2.append(sensitivity)
        n_imp2.append(n_ph_th / SNR_temp)

        freq_axis_fft, v_fft = extract_1d(run_id, data_1d_name="V_fft_avg_avg", setpoint_name='freq_param', plot=False, return_exp_name=False)
        bg_mask_fft = np.abs(freq_axis_fft - frequency) > 10e3
        v_fft_sub = v_fft - np.mean(v_fft[bg_mask_fft])
        signal_mask_fft = np.abs(freq_axis_fft - frequency) <= 3e3
        vrms_fft2.append(np.sqrt(np.sum(v_fft_sub[signal_mask_fft]**2)))

        if fit:
            freq_axis, spectrum = extract_1d(run_id, data_1d_name="avg_avg_psd_nodrive_substracted", setpoint_name='freq_param', plot=False, return_exp_name=False)
            freq_axis, spectrum_unsubstracted = extract_1d(run_id, data_1d_name="avg_avg_psd_nodrive", setpoint_name='freq_param', plot=False, return_exp_name=False)

            initial_guess = [frequency, Gamma_guess, max(centered_moving_average(spectrum, n=5)), 0]
            try:
                popt, pcov = scp.optimize.curve_fit(lorentzian_fkt, freq_axis, centered_moving_average(spectrum, n=5), p0=initial_guess)
                lorentzian, l_area = lorentzian_fkt_w_area(freq_axis, popt[0], popt[1], popt[2], popt[3])
                lorentzian = lorentzian_fkt(freq_axis, popt[0], popt[1], popt[2], popt[3])
                mask = np.abs(freq_axis - popt[0]) <= 3e3
                residuals_l = analyze_fit_quality_l(freq_axis[mask], spectrum[mask], popt)
                area_error_l = sum(abs(residuals_l))
                L_errors2.append(area_error_l)
                linewidths2l.append(abs(popt[1]))

                popt, pcov = scp.optimize.curve_fit(gaussian_fkt, freq_axis, centered_moving_average(spectrum, n=5), p0=initial_guess)
                gaussian, g_area = gaussian_fkt_w_area(freq_axis, popt[0], popt[1], popt[2], popt[3])
                residuals_g = analyze_fit_quality_g(freq_axis[mask], spectrum[mask], popt)
                area_error_g = sum(abs(residuals_g))
                G_errors2.append(area_error_g)
                linewidths2g.append(abs(popt[1]))
            except:
                popt = initial_guess
                print("fitting error")

    with open(cache_file2, 'wb') as f:
        pickle.dump({'areas2': areas2, 'g2_voltages2': g2_voltages2, 'frequencies2': frequencies2,
                     'sensitivities2': sensitivities2, 'n_imp2': n_imp2, 'linewidths2l': linewidths2l,
                     'linewidths2g': linewidths2g, 'L_errors2': L_errors2, 'G_errors2': G_errors2,
                     'vrms_fft2': vrms_fft2}, f)
    print(f"Dataset 2 cached to {cache_file2}")


# ─── Plotting ────────────────────────────────────────────────────────────────

#plt.plot(g2_voltages1,areas1,'g*')
#plt.plot(e_nrs_manual,manual_spectra_sum,'r*')
plt.plot(g2_voltages2, areas2, 'r*')
plt.title("plot1: areas_vs_g2_voltages untensioned", fontsize=14)
plt.xlabel("g2voltage")
plt.ylabel("PSD area")
for i, (x, y, run_id) in enumerate(zip(g2_voltages2, areas2, run_ids2)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
plt.show()

if e_nr == True:
    #plt.plot(electron_nrs1,areas1,'g*')
    plt.errorbar(electron_nrs2, areas2, yerr=G_errors2, fmt='g*', capsize=5, label='raw areas with G_errors, untensioned')
    #add here G_errors for the first few datapoints!!!
    #plt.errorbar(electron_nrs1, areas1, yerr=L_errors1, fmt='g*', capsize=5, label='Areas with Lerrors')
    plt.errorbar(e_nrs_manual, manual_spectra_sum, yerr=manual_errors, fmt='r*', capsize=5, label='raw areas with G_errors, tensioned')
    plt.title("plot2: areas_vs_e_nr with best sensitivity ref linesweep 1946", fontsize=14)
    plt.xlabel("nr electrons")
    plt.ylabel("PSD area")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for i, (x, y, run_id) in enumerate(zip(electron_nrs2, areas2, run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    plt.legend()
    plt.show()

#plt.plot(g2_voltages1,frequencies1,'g*')
plt.plot(g2_voltages2, frequencies2, 'r*')
plt.title("frequencies_vs_e_nr with best sensitivity ref linesweep 1946 (repeated)")
plt.xlabel("g2voltage")
plt.ylabel("frequency")
for i, (x, y, run_id) in enumerate(zip(g2_voltages2, frequencies2, run_ids2)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
plt.show()

if e_nr == True:
    #plt.plot(electron_nrs1,frequencies1,'g*')
    plt.plot(electron_nrs2, frequencies2, 'g*')
    plt.plot(e_nrs_manual, manual_frequencies, 'r*')
    plt.title("plot3: frequencies_vs_e_nr with best sensitivity ref linesweep 1946(repeated)")
    plt.xlabel("nr electrons")
    plt.ylabel("frequency")
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for i, (x, y, run_id) in enumerate(zip(electron_nrs2, frequencies2, run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    for i, (x, y, run_id) in enumerate(zip(electron_nrs1, frequencies1, run_ids1)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    plt.show()

if e_nr == True:
    #plt.plot(electron_nrs1,sensitivities1,'g*')
    plt.plot(electron_nrs2, sensitivities2, 'g*')
    plt.plot(e_nrs_manual, manual_sens, 'r*')
    plt.title("plot4: sensitivities")
    plt.xlabel("nr electrons")
    plt.ylabel("sensitivities")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    for i, (x, y, run_id) in enumerate(zip(electron_nrs2, sensitivities2, run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    #for i, (x, y,run_id) in enumerate(zip(electron_nrs1, areas1,run_ids1)):
    #    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    plt.show()


    #plt.plot(electron_nrs1, np.sqrt(areas1)*np.array(sensitivities1), 'g*')
    #plt.plot(electron_nrs2, np.sqrt(areas2)*np.array(sensitivities2), 'r*')


    #################ERRORBARS##################

    errors_scaled2 = 1/np.array(sensitivities2) * np.array(G_errors2) / (2 * np.sqrt(areas2))
    errors_scaled1 = 1/np.array(sensitivities1) * np.array(G_errors1) / (2 * np.sqrt(areas1))
    manual_errors_scaled = 1/np.array(manual_sens) * np.array(manual_errors) / (2 * np.sqrt(manual_spectra_sum))
    regression_area = copy.copy(manual_spectra_sum)
    regression_area.extend(areas2[0:12])
    e_nrs_reg = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    reg_sens = copy.copy(manual_sens)
    reg_sens.extend(sensitivities2[0:12])

    # Remove datapoint at 12 electrons from regression data
    idx_12 = e_nrs_reg.index(12)
    e_nrs_reg.pop(idx_12)
    regression_area.pop(idx_12)
    reg_sens.pop(idx_12)

    # Remove datapoint at 12 electrons from plot data
    electron_nrs2_plot = copy.copy(electron_nrs2)
    areas2_plot = copy.copy(areas2)
    sensitivities2_plot = copy.copy(sensitivities2)
    errors_scaled2_plot = list(errors_scaled2)
    linewidths2l_plot = copy.copy(linewidths2l)

    idx_12_plot = electron_nrs2_plot.index(12)
    electron_nrs2_plot.pop(idx_12_plot)
    areas2_plot.pop(idx_12_plot)
    sensitivities2_plot.pop(idx_12_plot)
    errors_scaled2_plot.pop(idx_12_plot)
    linewidths2l_plot.pop(idx_12_plot)

    # Create figure with two y-axes
    fig, ax1 = plt.subplots()

    # Left y-axis: scaled sqrt areas (all green, no distinction) - divided by 1e3
    ax1.errorbar(electron_nrs2_plot, np.sqrt(areas2_plot)/np.array(sensitivities2_plot)/1e3, yerr=np.array(errors_scaled2_plot)/1e3, fmt='g*', capsize=5)
    ax1.errorbar(e_nrs_manual, np.sqrt(manual_spectra_sum)/np.array(manual_sens)/1e3, yerr=np.array(manual_errors_scaled)/1e3, fmt='g*', capsize=5)

    # Calculate linear regression through origin
    x1 = np.array(e_nrs_reg)
    y1 = np.sqrt(regression_area) / np.array(reg_sens)
    slope = np.sum(x1 * y1) / np.sum(x1**2)
    x1_line = np.linspace(0, max(e_nrs_reg), 100)
    y1_line = slope * x1_line / 1e3  # also divide by 1e3
    ax1.plot(x1_line, y1_line, 'g-', linewidth=2)

    #ax1.text(0.05, 0.95, f'y = {slope:.4e}x', transform=ax1.transAxes, fontsize=10, verticalalignment='top', color='g')

    ax1.set_xlabel("nr electrons")
    ax1.set_ylabel("scaled psd sqrt(Svv)*dVg/dG", color='g')
    ax1.tick_params(axis='y', labelcolor='g')
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

    # Right y-axis: Lorentzian linewidths (also with 12 electron point removed)
    ax2 = ax1.twinx()
    ax2.plot(electron_nrs2_plot, linewidths2l_plot, 'b^', label='linewidth (Lorentzian)')
    ax2.set_ylabel("linewidth (Hz)", color='b')
    ax2.tick_params(axis='y', labelcolor='b')
    ax2.set_ylim(bottom=0, top=2000)

    plt.title("plot5: scaled sqrt areas, gaussian error, regression full data")
    plt.tight_layout()
    plt.show()

    plt.plot(electron_nrs1, n_imp1, '*b')
    plt.plot(electron_nrs2, n_imp2, '*g')
    plt.title("plot6: n_imp vs n_e", fontsize=14)
    plt.xlabel("nr electrons")
    plt.ylabel("n_imp")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.legend()
    plt.show()

    idx = run_ids2.index(661)
    linewidths2l.pop(idx)
    linewidths2g.pop(idx)
    electron_nrs2.pop(idx)

    #idx = run_ids2.index(540)
    #linewidths2.pop(idx)
    #electron_nrs2.pop(idx)

    plt.plot(electron_nrs2, linewidths2l, 'g*')
    #plt.plot(electron_nrs1, linewidths1l, 'g*')
    plt.title("linewidths Lorentzian fit")
    plt.xlabel("nr electrons")
    plt.ylabel("linewidths")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()

    plt.plot(electron_nrs2, linewidths2g, 'g*')
    #plt.plot(electron_nrs1, linewidths1g, 'g*')
    plt.title("linewidths Gaussian fit")
    plt.xlabel("nr electrons")
    plt.ylabel("linewidths")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()

    #############here add calculation for coupling strength##############

    # ─── V_rms plot ──────────────────────────────────────────────────────────

    vrms1 = np.array([psd_area_to_vrms(a) for a in areas1])
    vrms2 = np.array([psd_area_to_vrms(a) for a in areas2_plot])
    vrms_manual = np.array([psd_area_to_vrms(a) for a in manual_spectra_sum])

    vrms_fft2_plot = list(vrms_fft2)
    vrms_fft2_plot.pop(idx_12_plot)

    plt.plot(electron_nrs2_plot, vrms2, 'g*', label='PSD-derived, untensioned (140 MHz)')
    plt.plot(e_nrs_manual, vrms_manual, 'r*', label='PSD-derived, tensioned')
    plt.plot(electron_nrs2_plot, vrms_fft2_plot, 'g^', label='FFT-derived, untensioned (140 MHz)')
    plt.plot(electron_nrs1, vrms_fft1, 'b^', label='FFT-derived, untensioned (150 MHz)')
    plt.title("V_rms vs nr electrons", fontsize=14)
    plt.xlabel("nr electrons")
    plt.ylabel("V_rms (V)")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    ax = plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for i, (x, y, run_id) in enumerate(zip(electron_nrs2_plot, vrms2, run_ids2)):
        plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
    plt.legend()
    plt.show()


#calculate delta_epsilon from (dVLI/dVg)^-1; get conversion factor of I_sens to dVg/dVLI - first, can extrapolate from any one value
#  lever arm=0.25
# slope (from db18 run 685): 4.29 V_LI/V_g so 17.1 V_LI/eV corresponds to ~13pA sensitivity
#so conversion factor is dVLIdepsdIsens 1.31 V_LI/eV/pA

#convert Area to voltage area (psd to V)
# convert sens to pA
#use n_ph_th, can you add a little function to convert the areas (eg area=metadata_temp['integral_over_substracted_psd']) from psd to voltage, assuming 50 ohm?


sensitivities1_pA = np.array(sensitivities1) * 1e12
sensitivities2_pA = np.array(sensitivities2) * 1e12
manual_sens_pA = np.array(manual_sens) * 1e12


dVLIdepsdIsens=1.31 #V_LI/eV/pA

dVLIdeps1=sensitivities1_pA*dVLIdepsdIsens
dVLIdeps2=np.array(sensitivities2_plot)*1e12*dVLIdepsdIsens
manual_dVLIdeps=manual_sens_pA*dVLIdepsdIsens

# ─── PSD-derived ─────────────────────────────────────────────────────────────
deps1_eV = vrms1 / dVLIdeps1
deps2_eV = vrms2 / dVLIdeps2
deps_manual_eV = vrms_manual / manual_dVLIdeps

eV_to_Hz = 1.602e-19 / 6.626e-34  # E = h*f
deps1_Hz = deps1_eV * eV_to_Hz
deps2_Hz = deps2_eV * eV_to_Hz
deps_manual_Hz = deps_manual_eV * eV_to_Hz

g2_Hz = deps2_Hz / np.sqrt(2 * n_ph_th)
g_manual_Hz = deps_manual_Hz / np.sqrt(2 * n_ph_th)

# ─── FFT-derived ──────────────────────────────────────────────────────────────
vrms_fft2_plot_arr = np.array(vrms_fft2_plot)
vrms_fft1_arr = np.array(vrms_fft1)

deps2_fft_eV = vrms_fft2_plot_arr / dVLIdeps2
deps1_fft_eV = vrms_fft1_arr / dVLIdeps1
deps2_fft_Hz = deps2_fft_eV * eV_to_Hz
deps1_fft_Hz = deps1_fft_eV * eV_to_Hz

g2_fft_Hz = deps2_fft_Hz / np.sqrt(2 * n_ph_th)
g1_fft_Hz = deps1_fft_Hz / np.sqrt(2 * n_ph_th)

# ─── delta_epsilon (MHz) ─────────────────────────────────────────────────────
plt.plot(electron_nrs2_plot, deps2_Hz / 1e6, 'g*', label='PSD, untensioned (140 MHz)')
plt.plot(e_nrs_manual, deps_manual_Hz / 1e6, 'r*', label='PSD, tensioned')
plt.plot(electron_nrs1, deps1_Hz / 1e6, 'b*', label='PSD, untensioned (150 MHz)')
plt.plot(electron_nrs2_plot, deps2_fft_Hz / 1e6, 'g^', label='FFT, untensioned (140 MHz)')
plt.plot(electron_nrs1, deps1_fft_Hz / 1e6, 'b^', label='FFT, untensioned (150 MHz)')
plt.title("delta_epsilon vs nr electrons", fontsize=14)
plt.xlabel("nr electrons")
plt.ylabel("delta_epsilon (MHz)")
plt.xlim(left=0)
plt.ylim(bottom=0)
ax = plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
for i, (x, y, run_id) in enumerate(zip(electron_nrs2_plot, deps2_Hz / 1e6, run_ids2)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
plt.legend()
plt.show()

# ─── delta_epsilon (eV) ──────────────────────────────────────────────────────
plt.plot(electron_nrs2_plot, deps2_eV, 'g*', label='PSD, untensioned (140 MHz)')
plt.plot(e_nrs_manual, deps_manual_eV, 'r*', label='PSD, tensioned')
plt.plot(electron_nrs1, deps1_eV, 'b*', label='PSD, untensioned (150 MHz)')
plt.plot(electron_nrs2_plot, deps2_fft_eV, 'g^', label='FFT, untensioned (140 MHz)')
plt.plot(electron_nrs1, deps1_fft_eV, 'b^', label='FFT, untensioned (150 MHz)')
plt.title("delta_epsilon vs nr electrons", fontsize=14)
plt.xlabel("nr electrons")
plt.ylabel("delta_epsilon (eV)")
plt.xlim(left=0)
plt.ylim(bottom=0)
ax = plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
for i, (x, y, run_id) in enumerate(zip(electron_nrs2_plot, deps2_eV, run_ids2)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
plt.legend()
plt.show()

# ─── g simple (PSD + FFT) ────────────────────────────────────────────────────
g1_Hz = deps1_Hz / np.sqrt(2 * n_ph_th)

plt.plot(electron_nrs2_plot, g2_Hz / 1e6, 'g*', label='PSD, untensioned (140 MHz)')
plt.plot(e_nrs_manual, g_manual_Hz / 1e6, 'r*', label='PSD, tensioned')
plt.plot(electron_nrs1, g1_Hz / 1e6, 'b*', label='PSD, untensioned (150 MHz)')
plt.plot(electron_nrs2_plot, g2_fft_Hz / 1e6, 'g^', label='FFT, untensioned (140 MHz)')
plt.plot(electron_nrs1, g1_fft_Hz / 1e6, 'b^', label='FFT, untensioned (150 MHz)')
plt.title("coupling strength g vs nr electrons", fontsize=14)
plt.xlabel("nr electrons")
plt.ylabel("g/2\u03c0 (MHz)")
plt.xlim(left=0)
plt.ylim(bottom=0)
ax = plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
for i, (x, y, run_id) in enumerate(zip(electron_nrs2_plot, g2_Hz / 1e6, run_ids2)):
    plt.text(x, y, f"{run_id}", fontsize=8, ha='left', va='bottom')
plt.legend()
plt.show()

# ─── g + linewidths: PSD-derived ─────────────────────────────────────────────
g_errors2_Hz = np.array(errors_scaled2_plot) * np.sqrt(50) / (1e12*dVLIdepsdIsens) * eV_to_Hz / np.sqrt(2*n_ph_th)
g_errors_manual_Hz = manual_errors_scaled * np.sqrt(50) / (1e12*dVLIdepsdIsens) * eV_to_Hz / np.sqrt(2*n_ph_th)

x1 = np.array(e_nrs_reg)
y1_g = (np.sqrt(np.array(regression_area) * 50) / (np.array(reg_sens)*1e12*dVLIdepsdIsens)
        * eV_to_Hz / np.sqrt(2*n_ph_th))
slope_g = np.sum(x1 * y1_g) / np.sum(x1**2)
x1_line = np.linspace(0, max(e_nrs_reg), 100)

fig, ax1 = plt.subplots()
ax1.errorbar(electron_nrs2_plot, g2_Hz / 1e6, yerr=g_errors2_Hz / 1e6, fmt='g*', capsize=5, label='untensioned (140 MHz)')
ax1.errorbar(e_nrs_manual, g_manual_Hz / 1e6, yerr=g_errors_manual_Hz / 1e6, fmt='g*', capsize=5, label='tensioned')
ax1.plot(x1_line, slope_g * x1_line / 1e6, 'g-', linewidth=2)
ax1.set_xlabel("nr electrons")
ax1.set_ylabel("g/2\u03c0 (MHz)", color='g')
ax1.tick_params(axis='y', labelcolor='g')
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
#ax1.legend(loc='upper left')
ax2 = ax1.twinx()
ax2.plot(electron_nrs2_plot, np.array(linewidths2l_plot) / 1e3, 'b^')
ax2.set_ylabel("linewidth (kHz)", color='b')
ax2.tick_params(axis='y', labelcolor='b')
ax2.set_ylim(bottom=0, top=2)
#plt.title("g/2\u03c0 PSD-derived + linewidths")
plt.title("coupling rate vs electron nr")
plt.tight_layout()
plt.show()

# ─── g + linewidths: FFT-derived ─────────────────────────────────────────────
x1_fft = np.array(electron_nrs2_plot)
slope_g_fft = np.sum(x1_fft * g2_fft_Hz) / np.sum(x1_fft**2)
x1_line_fft = np.linspace(0, max(electron_nrs2_plot), 100)

fig, ax1 = plt.subplots()
ax1.plot(electron_nrs2_plot, g2_fft_Hz / 1e6, 'g^', label='untensioned (140 MHz)')
ax1.plot(e_nrs_manual, g_manual_Hz / 1e6, 'r*', label='tensioned')
ax1.plot(x1_line_fft, slope_g_fft * x1_line_fft / 1e6, 'g-', linewidth=2)
ax1.set_xlabel("nr electrons")
ax1.set_ylabel("g/2\u03c0 (MHz)", color='g')
ax1.tick_params(axis='y', labelcolor='g')
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.legend(loc='upper left')
ax2 = ax1.twinx()
ax2.plot(electron_nrs2_plot, np.array(linewidths2l_plot) / 1e3, 'b^')
ax2.set_ylabel("linewidth (kHz)", color='b')
ax2.tick_params(axis='y', labelcolor='b')
ax2.set_ylim(bottom=0, top=2)
plt.title("g/2\u03c0 FFT-derived + linewidths")
plt.tight_layout()
plt.show()
