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
from prelude import make_base_dir, distance, calculate_distance, closest_encounter, calc_time

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
			speed_knots = flight_data['speed'][line]
			speed_mps = speed_knots * 0.514444
			alt_ft = flight_data['altitude'][line]
			alt_m = alt_ft * 0.3048
			for y in range(len(station)):
				if str(station[y]) == str(sta[n]):
					dist_km = distance(seismo_latitudes[y], seismo_longitudes[y], flight_latitudes[line], flight_longitudes[line])	
					dist_m = dist_km * 1000
					
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
					fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6)) #, gridspec_kw={'height_ratios': [3, 1]})     

					ax1.plot(t, data, 'k', linewidth=0.5)
					ax1.set_title(title)
					ax1.axvline(x=tim, c = 'c', ls = '--')
					ax1.margins(x=0)

					spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))

					# Find the index of the middle frequency
					middle_index = len(times) // 2
					middle_column = spec[:, middle_index]
					vmin = 0  
					vmax = np.max(middle_column) 
					if n == 0:
						tprime0 = 112
						fnot = [93, 115, 153, 172, 228]
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 68
						l = 2135

					if n == 1:
						fnot = [37, 56, 73, 110, 146, 165, 182, 218, 238, 256, 275]
						tprime0 = 107
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 100
						l = 2700

					if n == 2:
						fnot = [79,131,261]
						tprime0 = 93
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 139
						l = 4650

					if n == 3:
						fnot = [36,73,121,136,144]
						tprime0 = 116
						tpr = np.arange(80, 170, 1)
						c = 343
						v0 = 142
						l = 2450

					if n == 4:
						fnot = [13,27,40,47,54,60,67,74,80,87,90,94,101,108,114,121,127,134,148,160,177,189,202,223,239,247,270]
						tprime0 = 140
						tpr = np.arange(40, 230, 1)
						c = 343
						v0 = 64
						l = 580

					for f0 in fnot:
						ft = []
						for tprime in tpr:
							ft0p = f0*1/(1+(v0/c)*(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
								
							ft.append(ft0p)
						ax2.plot(tpr, ft, 'g', linewidth=0.5)
							
					# Plot spectrogram
					cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
					ax2.set_xlabel('Time [s]')
					v0 = speed * 0.000514444
					ax2.axvline(x=tim, c = 'c', ls = '--', label='Wave generated (t0): '+str(tim)+' sec')
					ax2.axvline(x=tprime0, c = 'g', ls = '--', label='Estimated t0: '+str(tprime0)+' sec')
					print(calc_time(tim,dist_m,alt_m))
					tarrive = calc_time(tim,dist_m,alt_m)
					ax2.axvline(x=calc_time(tim,dist_m,alt_m), c = 'r', ls = '--',label='Wave arrvial: '+str(int(tarrive-120))+' sec after t0')
					ax2.legend(loc='upper right',fontsize = 'x-small')
					ax2.set_ylabel('Frequency (Hz)')
					ax2.set_title("Forward Model: t'= "+str(tprime0)+' sec, v0 = '+str(int(speed_mps))+' m/s, l = '+str(int(dist_m))+' m, \n' + 'f0 = '+str(fnot)+' Hz', fontsize='x-small')
					ax2.margins(x=0)
					ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

					plt.colorbar(mappable=cax, cax=ax3)
					ax3.set_ylabel('Relative Amplitude (dB)')

					ax2.margins(x=0)
					ax2.set_xlim(0, 240)
					
					# Plot overlay
					spec2 = 10 * np.log10(MDF)
					middle_column2 = spec2[:, middle_index]
					vmin = np.min(middle_column2)
					vmax = np.max(middle_column2)
				
					# Create ax4 and plot on the same y-axis as ax2
					ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2) #, width=vmax*1.1-vmin, height=int(fs/2))
					ax4.plot(middle_column2, frequencies, c='orange')  
					ax4.set_ylim(0, int(fs/2))
					ax4.set_xlim(vmax*1.1, vmin) #, width=vmax*1.1-vmin, height=int(fs/2))
					ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
					ax4.grid(axis='y')
					plt.show()
					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-02-'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'/'
					make_base_dir(BASE_DIR)
					fig.savefig('/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-02-'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'/'+str(time[n])+'_'+str(flight_num[n])+'.png')
					plt.close()


					'''
										l = closest_encounter(flight_latitudes, flight_longitudes,line, tm, seismo_latitudes[y], seismo_longitudes[y])
										start_time = tr[2].stats.starttime + (mins * 60) + secs - tim


										tpr = np.arange(40, 230, 1)
										middle_index = len(times) // 2

										for tnew in tpr:
											column = spec[50:250, tnew]
											peaks, _ = signal.find_peaks(column, prominence=10)
											sorted_peaks = sorted(peaks, key=lambda x: column[x], reverse=True)[:12]  # Select the 12 largest peaks
											
											# Scatter plot
											ax2.scatter([tnew] * len(sorted_peaks), sorted_peaks + 50, marker='x', color='k', linewidth=0.5)
										plt.show()
					'''
					# Find the center of the trace
		
					#plt.figure()
					# Spectrogram 
					#plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
				

					# '''
					# 					coords = []

					# 					def onclick(event):
					# 						global coords
					# 						coords.append((event.xdata, event.ydata))
					# 						plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
					# 						plt.draw() 
					# 						print('Clicked:', event.xdata, event.ydata)  
					# 					cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

					# 					plt.show(block=True)
					# 					# Convert the list of coordinates to a numpy array
					# 					coords_array = np.array(coords)
					# 					print(coords_array)
					# 					v0 = speed*0.514444
					# 					c = 343
					# 					for jj in range(len(coords_array)):
					# 						print(coords_array[jj][0], coords_array[jj][1])
					# 						tflight = coords_array[jj][0]-(np.sqrt(coords_array[jj][0]**2-(1-v0**2/c**2)*(coords_array[jj][0]**2-l**2/c**2)))/(1-v0**2/c**2)
					# 						f0 = coords_array[jj][1]*(1+(v0/c)*(v0*tflight/(np.sqrt(l**2+(v0*tflight)**2))))
					# 						print('f0 = ', f0)
											
					# 					# Save the coordinates to a text file
					# 					#np.savetxt('home/irseppi/REPOSITORIES/parkshwynodal/coords_'+str(flight_num[n])+'.txt', coords_array)
									
					# 					plt.figure()
					# 					plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
					# 					plt.plot(coords_array[:,0], coords_array[:,1])
					# 					plt.show()
					# '''
