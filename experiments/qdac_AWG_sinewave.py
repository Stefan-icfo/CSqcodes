import qcodes as qc
from instruments import qdac
import time


# Access Channel 6
channel = qdac.ch06

# Define amplitude and frequency for the sine wave
amplitude = 1e-3  # peaktopeak
frequency = 1e3    # kHz

# Configure the sine wave on Channel 8
sine_wave_context = channel.sine_wave(
    frequency_Hz=frequency,  # Set frequency to 5 kHz
    span_V=amplitude,        # Set amplitude span to 10 mV
    offset_V=0,            # No offset, centered around 0V
    repetitions=-1,          # Run indefinitely (-1 for infinite repetitions)
)

# Start the sine wave
time.sleep(1)
sine_wave_context.start()

# Let it run for a while (e.g., 10 seconds)
time.sleep(1)

# To stop the sine wave, use the abort method
sine_wave_context.abort()


