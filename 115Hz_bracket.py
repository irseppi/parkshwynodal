import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import datetime
import obspy
from scipy.signal import spectrogram

ht = datetime.datetime.utcfromtimestamp(1550089044)
mins = ht.minute
secs = ht.second

tim = 120

p = "/scratch/naalexeev/NODAL/2019-02-13T20:00:00.000000Z.2019-02-13T21:00:00.000000Z.1004.mseed"
tr = obspy.read(p)

tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
data = tr[2][0:-1]
fs = int(tr[2].stats.sampling_rate)
title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						

# Compute spectrogram
f, t, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 


frequency = 115  # Frequency of interest (115 Hz)

# Extract the relevant frequency bin (closest to 115 Hz)
freq_index = np.argmin(np.abs(f - frequency))

# Plot the histogram of power values in that frequency bin
plt.hist(10*np.log10(Sxx[freq_index]), bins=50, color='skyblue', edgecolor='black')

plt.xlabel('Power')
plt.ylabel('Count')
plt.title(f'Histogram of Power at {frequency} Hz')
plt.grid(True)
plt.show()

