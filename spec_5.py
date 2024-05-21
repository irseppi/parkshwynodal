import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import obspy
from scipy.signal import spectrogram
from scipy import signal
from obspy.geodetics import gps2dist_azimuth
from datetime import datetime, timedelta
from prelude import make_base_dir, dist_less, calc_time

# Load the seismometer location data
seismo_data = pd.read_csv('input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
file_names = ['/scratch/irseppi/nodal_data/flightradar24/20190225_positions/20190225_530342801.csv','/scratch/irseppi/nodal_data/flightradar24/20190214_positions/20190214_528485724.csv','/scratch/irseppi/nodal_data/flightradar24/20190214_positions/20190214_528473220.csv','/scratch/irseppi/nodal_data/flightradar24/20190214_positions/20190214_528407493.csv','/scratch/irseppi/nodal_data/flightradar24/20190213_positions/20190213_528293430.csv']

flight_num = [530342801,528485724,528473220,528407493,528293430,527937367,529741194,529776675,529179112,530165646]
time = [1551066051,1550172833,1550168070,1550165577,1550089044,1549912188,1550773710,1550787637,1550511447,1550974151]
sta = [1022,1272,1173,1283,1004,"CCB","F6TP","F4TN","F3TN","F7TV"]
day = [25,14,14,14,13,11,21,21,18,24]


for i in range(0,6):
	flight_data = pd.read_csv(file_names[i], sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	timestamp = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']
	head = flight_data['heading']

	ht = datetime.utcfromtimestamp(time[i])
	h = ht.hour
	month = ht.month
	day = ht.day
	mins = ht.minute
	secs = ht.second
	
	month2 = str(month)
	if month == 3 and day < 10:
		day1 = '0'+str(day)
	else:
		day1 = str(day)
	h_u = str(h+1)
	if h < 23:
			
		day2 = day1
		if h < 10:
			h_u = '0'+str(h+1)
			h = '0'+str(h)
		else:
			h_u = str(h+1)
			h = str(h)
	else:
		h_u = '00'
		if month == '02' and day == '28':
			month2 = '03'
			day2 = '01'
		else:
			day2 = str(day + 1)
				
	
	tim = 120	

		
	n = "/scratch/naalexeev/NODAL/2019-0"+ str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month2)+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(sta[i])+".mseed"

	if os.path.isfile(n):
		print('made it')
		tr = obspy.read(n)
		tr[2].trim(tr[2].stats.starttime + (mins*60) + secs - tim, tr[2].stats.starttime + (mins*60) + secs + tim)
		data = tr[2][0:-1]
		fs = int(tr[2].stats.sampling_rate)
		title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'	
	
		# Time array
		t = np.arange(len(data)) / fs
		g = fs*240

		# Compute spectrogram
		frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 
		a, b = Sxx.shape
		
		MDF = np.zeros((a,b))
		for row in range(len(Sxx)):
			m = len(Sxx[row])
			p = sorted(Sxx[row])
			median = p[int(m/2)]

			for col in range(m):
				MDF[row][col] = median

		fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))
		
		# Find the center of the trace
		center_index = len(data) // 2
		center_time = t[center_index]

		ax1.plot(t, data, 'k', linewidth=0.5)
		ax1.set_title(title)
		#ax1.axvline(x=time, c = 'r', ls = '--')
		#ax1.axvline(x=tm, c = 'r', ls = '--')
		ax1.margins(x=0)
		
		spec = (10 * np.log10(Sxx)) - (10 * np.log10(MDF))

		# Find the index of the middle frequency
		middle_index = len(times) // 2

		# Extract the middle line of the spectrogram
		middle_column = spec[:, middle_index]

		vmin = np.min(middle_column)
		vmax = np.max(middle_column)
		# Plot spectrogram
		cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='hsv', vmin = vmin, vmax = vmax)



		ax2.axvline(x=center_time, c = 'k', ls = '--')
		ax2.set_ylabel('Frequency (Hz)')
		ax2.margins(x=0)
		ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

		plt.colorbar(mappable=cax, cax=ax3)
		ax3.set_ylabel('Relative Amplitude (dB)')

		#make_base_dir('/scratch/irseppi/nodal_data/Plane_map_spec/')
		#fig.savefig('/scratch/irseppi/nodal_data/P_map_spec/spec_'+str(ctime)+'_'+str(csta)+'_'+str(flight_num[i])+'.png')
		
		plt.close()
		
		peaks, _ = signal.find_peaks(middle_column, prominence=20) #, distance = 10) 
		np.diff(peaks)
		fig = plt.figure(figsize=(10,6))
		plt.grid()
		
		plt.plot(frequencies, middle_column, c='c')
		plt.plot(peaks, middle_column[peaks], "x")
		for g in range(len(peaks)):
			plt.text(peaks[g], middle_column[peaks[g]], peaks[g])
		plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))

		plt.xlim(0,int(fs/2))
		plt.ylim(vmin,vmax*1.1)
		plt.xlabel('Freq [Hz]')
		plt.ylabel('Amplitude [dB]')
		plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))
		plt.show()
		#make_base_dir('/scratch/irseppi/nodal_data/P_map_spec/')
		#fig.savefig('/scratch/irseppi/nodal_data/P_map_spec/fft_'+str(ctime)+'_'+str(csta)+'_'+str(flight_num[i])+'.png')
		plt.close()

		spec2 = 10 * np.log10(MDF)

		middle_column2 = spec2[:, middle_index]
		vmin = np.min(middle_column2)
		vmax = np.max(middle_column2)

		fig = plt.figure(figsize=(1.5,3))
		plt.margins(x=0)
		plt.plot(middle_column2,frequencies, c='c')
		plt.ylim(0,int(fs/2))
		plt.xlim(vmax*1.1,vmin)
		plt.tick_params(left = False, right = False , labelleft = False , 
	labelbottom = False, bottom = False)

		#BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5spec2/201902'+str(day[n])+'/'+str(flight_num[n])+'/'+sta[i]+'/'
		#make_base_dir(BASE_DIR)
		#fig.savefig('/scratch/irseppi/nodal_data/plane_info/5spec2/201902'+str(day[n])+'/'+str(flight_num[n])+'/'+sta[i]+'/'+sta[i]+'_'+str(time[n])+'.png', bbox_inches='tight',pad_inches =0)
		#plt.close()
	
