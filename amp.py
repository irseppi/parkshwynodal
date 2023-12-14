import pandas as pd
import matplotlib.pyplot as plt
import os

import datetime
import pytz

import math
import matplotlib.dates as mdates

# Import all dependencies
import numpy as np
import datetime
from matplotlib import dates, pyplot as plt
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from scipy.signal import spectrogram
import obspy
from pathlib import Path

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()


text = open('input/all_station_crossing_db.txt', 'r')
window_duration = 10  # spectrogram window duration [s]
freq_lims = (0.5,250)  # frequency limits for output spectrogram. If `None`, the limits will be adaptive
v_percent_lims = (20,97.5)  # colorbar limits
for line in text.readlines():
	val = line.split(',')
	ht = datetime.datetime.utcfromtimestamp(int(val[2]))
	if ht.month == 2 and ht.day >= 11 or ht.month == 3:
		day_of_year = str((ht - datetime.datetime(2019, 1, 1)).days + 1)
		if val[5].isdigit() == False:
			try:
				
				n = "/aec/wf/2019/0"+day_of_year+"/"+str(val[5])+".HHZ.20190"+day_of_year+"000000+"
				print(n)
				tr = obspy.read(n)

				h = ht.hour
				mins = ht.minute
				secs = ht.second

				tm = 300
				tr[0].trim(tr[0].stats.starttime + (h *60 *60) + (mins * 60) + secs - tm, tr[0].stats.starttime + (h *60 *60) +  (mins * 60) + secs + tm)
				tr = tr[0]
				# Define reference value for spectrogram
				REFERENCE_VALUE = 1  # Velocity, in m/s

				# Extract trace information for FFT
				sampling_rate = tr.stats.sampling_rate
				samples_per_segment = int(window_duration * sampling_rate)

				# Compute spectrogram (Note that overlap is 90% of samples_per_segment)
				sample_frequencies, segment_times, spec = spectrogram(tr.data, sampling_rate, window='hann', scaling='density',
								                      nperseg=samples_per_segment, noverlap=samples_per_segment * .9)

				# Convert spectrogram matrix to decibels for plotting
				spec_db = 10 * np.log10(abs(spec) / (REFERENCE_VALUE ** 2))

				# Convert trace times to matplotlib dates
				trace_time_matplotlib = tr.stats.starttime.matplotlib_date + (segment_times / dates.SEC_PER_DAY)

				# Plot simple
				fig, ax = plt.subplots(3,1)
				# Plot trace
				ax[0].plot(tr.times(), tr.data, 'k-')
				# Plot spectrogram
				spec_db_plot = spec_db[np.flatnonzero((sample_frequencies > freq_lims[0]) & (sample_frequencies < freq_lims[1])), :]
				ax[1].imshow(spec_db, extent=[trace_time_matplotlib[0], trace_time_matplotlib[-1], sample_frequencies[0],
								sample_frequencies[-1]],
					       vmin=np.percentile(spec_db_plot, v_percent_lims[0]),
					       vmax=np.percentile(spec_db_plot, v_percent_lims[1]),
					       origin='lower', aspect='auto', interpolation='None')
				# Plot every 20th column of spectrogram to see how it changes with time
				#ax[1].axvline(val[2], linestyle='--')
				#ax[2].plot(np.arange(freq_lims[0]+0.05, freq_lims[1]-0.05, 0.1), spec_db_plot[:, i])
				for i in range(0, np.shape(spec_db_plot)[1], 3):
					ax[1].axvline(trace_time_matplotlib[i], linestyle='--')
					ax[2].plot(np.arange(freq_lims[0]+0.05, freq_lims[1]-0.05, 0.1), spec_db_plot[:, i])

				#fig.show()
				
				BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/allspec/'+ str(val[0]) + '/'+str(val[1])+'/'+str(val[5])+'/'
				make_base_dir(BASE_DIR)
				plt.savefig('/scratch/irseppi/nodal_data/plane_info/allspec/'+ str(val[0]) + '/'+str(val[1])+'/'+str(val[5])+'/'+str(val[1])+'_' + str(val[2]) + '.png')
				plt.close()
			except:
				print(day_of_year+' '+val[5])
				continue
	text.close()

