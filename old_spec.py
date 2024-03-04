import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy import signal
import obspy
import math
from obspy.core import UTCDateTime
import datetime
from obspy.geodetics import gps2dist_azimuth
from prelude import make_base_dir, distance
import datetime
from obspy import UTCDateTime
from prelude import make_base_dir, distance, closest_encounter
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station = seismo_data['Station']
flight_num = [530342801,528485724,528473220,528407493,528293430]
time = [1551066051,1550172833,1550168070,1550165577,1550089044]
sta = [1022,1272,1173,1283,1004]
day = [25,14,14,14,13]
for n in range(0,5):
	ht = datetime.datetime.utcfromtimestamp(time[n])
	mins = ht.minute
	secs = ht.second
	h = ht.hour
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
	flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/201902'+str(day[n])+'_positions/201902'+str(day[n])+'_'+str(flight_num[n])+'.csv', sep=",")

	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	tm = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']
	head = flight_data['heading']
	for line in range(len(tm)):
		if str(tm[line]) == str(time[n]):
			speed = flight_data['speed'][line]
			alt = flight_data['altitude'][line]
			for y in range(len(station)):
				if str(station[y]) == str(sta[n]):
					dist = distance(seismo_latitudes[y], seismo_longitudes[y], flight_latitudes[line], flight_longitudes[line])	

					p = "/scratch/naalexeev/NODAL/2019-02-"+str(day[n])+"T"+str(h)+":00:00.000000Z.2019-02-"+day2+"T"+h_u+":00:00.000000Z."+station[y]+".mseed"
					tr = obspy.read(p)
					tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
					data = tr[2][0:-1]
					fs = int(tr[2].stats.sampling_rate)
					title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
					t = tr[2].times()
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

					ax1.plot(t, data, 'k', linewidth=0.5)
					ax1.set_title(title)
					ax1.axvline(x=tim, c = 'r', ls = '--')
					ax1.margins(x=0)

					spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))

					# Find the index of the middle frequency
					middle_index = len(times) // 2
					middle_column = spec[:, middle_index]
					vmin = 0  
					vmax = np.max(middle_column) 
					
					# Plot spectrogram
					cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
					ax2.set_xlabel('Time [s]')
				
					freq1 = []
					time1 = []

					l = closest_encounter(flight_latitudes, flight_longitudes,line, tm, seismo_latitudes[y], seismo_longitudes[y])
					start_time = tr[2].stats.starttime + (mins * 60) + secs - tim

					# Convert ObsPy time to timestamp
					start_time = UTCDateTime(start_time)
					start_time= start_time.datetime.timestamp()

					end_time = tr[2].stats.starttime + (mins * 60) + secs + tim
					end_time  = UTCDateTime(end_time)
					end_time = end_time.datetime.timestamp()
					x_values = []
					y_values = []
					for gg in range(len(times)):
						column = spec[50:250, gg]
						
						tnew = times[gg]
					
						peaks, _ = signal.find_peaks(column, prominence=10)

						np.diff(peaks)
			
						largest_peak_index = np.argmax(column[peaks])

						# Scatter plot
						ax2.scatter(tnew, peaks[largest_peak_index]+50, marker='x', color='k', linewidth=0.5)
					
						freq1.append(peaks[largest_peak_index]+50)
						time1.append(tnew)
					
					#ax2.axvline(x=center_time, c = 'k', ls = '--')
					ax2.set_ylabel('Frequency (Hz)')
					ax2.margins(x=0)
					ax2.set_xlim(0, 240)
					ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

					plt.colorbar(mappable=cax, cax=ax3)
					ax3.set_ylabel('Relative Amplitude (dB)')
					plt.show()
				
					plt.figure()
					# Spectrogram 
					plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
				

					coords = []

					def onclick(event):
						global coords
						coords.append((event.xdata, event.ydata))
						plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
						plt.draw() 
						print('Clicked:', event.xdata, event.ydata)  
					cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

					plt.show(block=True)
					# Convert the list of coordinates to a numpy array
					coords_array = np.array(coords)
					print(coords_array)
					v0 = speed*0.514444
					c = 343
					for jj in range(len(coords_array)):
						print(coords_array[jj][0], coords_array[jj][1])
						tflight = coords_array[jj][0]-(np.sqrt(coords_array[jj][0]**2-(1-v0**2/c**2)*(coords_array[jj][0]**2-l**2/c**2)))/(1-v0**2/c**2)
						f0 = coords_array[jj][1]*(1+(v0/c)*(v0*tflight/(np.sqrt(l**2+(v0*tflight)**2))))
						print('f0 = ', f0)
						
					# Save the coordinates to a text file
					#np.savetxt('home/irseppi/REPOSITORIES/parkshwynodal/coords_'+str(flight_num[n])+'.txt', coords_array)
				
					plt.figure()
					plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
					plt.plot(coords_array[:,0], coords_array[:,1])
					plt.show()
			
