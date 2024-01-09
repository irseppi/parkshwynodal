import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
import obspy
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import os
import obspy
from obspy.core import UTCDateTime
import datetime
from numpy.fft import fft, ifft

#20190213,528271370,1550087918,16025,280,1183,AT73

ht = datetime.datetime.utcfromtimestamp(1550005905)
mins = ht.minute
secs = ht.second
h = ht.hour
station = str(1011)
tim = 120	
h_u = str(h+1)
h = str(h)
	
n = "/scratch/naalexeev/NODAL/2019-02-12T"+str(h)+":00:00.000000Z.2019-02-12T"+h_u+":00:00.000000Z."+station+".mseed"
tr = obspy.read(n)

tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
data = tr[2][0:-1]
fs = int(tr[2].stats.sampling_rate)
title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
t = tr[2].times()

# Time array
t = np.arange(len(data)) / fs
g = fs*240

# Compute spectrogram
frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9) 

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))     

ax1.plot(t, data, 'k', linewidth=0.5)
ax1.set_title(title)
ax1.axvline(x=tim, c = 'r', ls = '--')
ax1.margins(x=0)
#vmin = np.min(np.abs(10 * np.log10(Sxx)))
vmax = np.abs(np.max(10 * np.log10(Sxx)))
vmin = - vmax * (1/5)
# Plot spectrogram
cax = ax2.pcolormesh(times, frequencies, np.abs(10 * np.log10(Sxx)), shading='gouraud', cmap='hsv_r', vmin=vmin, vmax=vmax)
ax2.set_xlabel('Time [s]')

# Find the center of the trace
center_index = len(data) // 2
center_time = t[center_index]

ax2.axvline(x=center_time, c = 'r', ls = '--')
ax2.set_ylabel('Frequency (Hz)')
ax2.margins(x=0)
ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

plt.colorbar(mappable=cax, cax=ax3)
ax3.set_ylabel('Relative Amplitude (dB)')


plt.figure()
# Compute and plot amplitude spectrum for the center of the trace
window = np.hanning(fs)  # one second window
center_data = data[int(center_index):int(center_index+fs)] 
freqs = np.fft.rfftfreq(fs, 1/fs)
fft = 10*np.log10(np.fft.rfft(center_data))
plt.plot(freqs, np.abs(fft))

plt.xlim(2,int(fs/2))


# Find the index of the middle frequency
middle_index = len(times) // 2

# Extract the middle line of the spectrogram
middle_column = Sxx[:, middle_index]

# Plot the middle line of the spectrogram
plt.figure(figsize=(10,6))
plt.grid()
plt.plot(frequencies, np.abs(10 * np.log10(middle_column)))
plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))

plt.xlim(2,int(fs/2))

plt.xlabel('Freq [Hz]')
plt.ylabel('Amplitude [dB]')
plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))
plt.show()
plt.close()
	
