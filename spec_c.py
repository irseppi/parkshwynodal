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
from pathlib import Path

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()
#20190225_530342801_1551066051_1022
#20190214_528485724_1550172833_1272
#20190214_528473229_1550168070_1173
#20190214_528407493_1550165577_1283
#20190213_528293430_1550089044_1004

flight_num = [530342801,528485724,528473229,528407493,528293430]
time = [1551066051,1550172833,1550168070,1550165577,1550089044]
sta = [1022,1272,1173,1283,1004]
day = [25,14,14,14,13]
for n in range(0,5):
	ht = datetime.datetime.utcfromtimestamp(time[n])
	mins = ht.minute
	secs = ht.second
	h = ht.hour
	station = str(sta[n])
	tim = 120	
	h_u = str(h+1)
	if h < 23:			
		day2 = str(day[n])
		if h < 10:
			h_u = '0'+str(h+1)
			h = '0'+str(h)
		else:
			h_u = str(h+1)
			h = str(h)
	else:
		h_u = '00'
		day2 = str(day[n]+1)

			
			
	n = "/scratch/naalexeev/NODAL/2019-02-"+str(day[n])+"T"+str(h)+":00:00.000000Z.2019-02-"+day2+"T"+h_u+":00:00.000000Z."+station+".mseed"
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
	vmin = np.min(10 * np.log10(Sxx))
	vmax = np.max(10 * np.log10(Sxx))

	# Plot spectrogram
	cax = ax2.pcolormesh(times, frequencies, 10 * np.log10(Sxx), shading='gouraud', cmap='hsv', vmin=vmin, vmax=vmax)
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
	BASE_DIR = "/scratch/irseppi/nodal_data/plane_info/5plane_spec/"
	make_base_dir(BASE_DIR)
	fig.savefig('/scratch/irseppi/nodal_data/plane_info/5plane_spec/201902'+str(day[n])+ '_'+station+'_'+str(time[n])+'_'+flight_num[n]+'.png')
	plt.close()
	# Find the index of the middle frequency
	middle_index = len(times) // 2

	# Extract the middle line of the spectrogram
	middle_column = Sxx[:, middle_index]

	# Plot the middle line of the spectrogram
	plt.figure(figsize=(10,6))
	plt.grid()
	plt.plot(frequencies, 10 * np.log10(middle_column))
	plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))

	plt.xlim(2,int(fs/2))

	plt.xlabel('Freq [Hz]')
	plt.ylabel('Amplitude [dB]')
	plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))
	BASE_DIR = "/scratch/irseppi/nodal_data/plane_info/5spec/"
	make_base_dir(BASE_DIR)
	fig.savefig('/scratch/irseppi/nodal_data/plane_info/5spec/201902'+str(day[n])+'_'+station+'-'+str(time[n])+'_'+flight_num[n]+'.png')
	plt.close()

