import matplotlib.pyplot as plt

# Integration times in nanoseconds
integration_times_ns = [29.09, 60.9, 102.6, 150.0, 3560.0]  # 3.56 microseconds = 3560 ns

# Fidelities in percent (from your screenshots)
fidelities = [100.0, 100.0, 100.0, 100.0, 100.0]

# SNR values (from your screenshots)
snr_values = [77.67, 149.09, 209.32, 226.17, 1256.89]

# -------- Figure 1: Fidelity vs integration time --------
plt.figure()
plt.plot(integration_times_ns, fidelities, marker='o')
plt.xlabel("Integration time (ns)")
plt.ylabel("Fidelity (%)")
plt.title("Fidelity vs Integration Time")


# -------- Figure 2: SNR vs integration time --------
plt.figure()
plt.plot(integration_times_ns, snr_values, marker='o')
plt.xlabel("Integration time (ns)")
plt.ylabel("SNR")
plt.title("SNR vs Integration Time")


# Show both figures (two separate windows)
plt.show()

