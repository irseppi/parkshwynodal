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
from scipy.stats import mode
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station = seismo_data['Station']

flight_num = [530342801,528485724,528473220,528407493,528293430,527937367,529741194,529776675,529179112,530165646,531605202,531715679,529805251,529948401,530122923]
time = [1551066051,1550172833,1550168070,1550165577,1550089044,1549912188,1550773710,1550787637,1550511447,1550974151,1551662362,1551736354,1550803701,1550867033,1550950429]
sta = [1022,1272,1173,1283,1004,"CCB","F6TP","F4TN","F3TN","F7TV",1010,1021,1006,1109,1298]
day = [25,14,14,14,13,11,21,21,18,24,4,4,22,22,23]
month = [2,2,2,2,2,2,2,2,2,2,3,3,2,2,2]
forward_model = False

for n in range(0,15):
	ht = datetime.datetime.utcfromtimestamp(time[n])
	
	mins = ht.minute
	secs = ht.second
	h = ht.hour
	tim = 120	
	h_u = h+1
	if h < 23:			
		day2 = day[n]
		if h < 10:
			h_u = '0'+str(h+1)
			h = '0'+str(h)
		else:
			h_u = h+1
			h = h
	else:
		h_u = '00'
		day2 = str(day[n]+1)
	if len(str(day[n])) == 1:
		day[n] = '0'+str(day[n])
		day2 = day[n]
	flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/20190'+str(month[n])+str(day[n])+'_positions/20190'+str(month[n])+str(day[n])+'_'+str(flight_num[n])+'.csv', sep=",")

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
					if isinstance(sta[n], str):
						day_of_year = str((ht - datetime.datetime(2019, 1, 1)).days + 1)
	
						p = "/aec/wf/2019/0"+day_of_year+"/"+str(sta[n])+".*Z.20190"+day_of_year+"000000+"
						tr = obspy.read(p)
						
						tr[0].trim(tr[0].stats.starttime +(int(h) *60 *60) + (mins * 60) + secs - tim, tr[0].stats.starttime +(int(h) *60 *60) + (mins * 60) + secs + tim)
						data = tr[0][0:-1]
						fs = int(tr[0].stats.sampling_rate)
						title    = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'						
						t                  = tr[0].times()
					else:
						p = "/scratch/naalexeev/NODAL/2019-0"+str(month[n])+"-"+str(day[n])+"T"+str(h)+":00:00.000000Z.2019-0"+str(month[n])+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(station[y])+".mseed"
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
					
					ty = False
					if ty == True:
						if isinstance(sta[n], str):
							spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))
						else:
							spec = np.zeros((a,b))
							for col in range(0,b):
								p = sorted(Sxx[:, col])
								median = p[int(len(p)/2)]

								for row in range(len(Sxx)):

									spec[row][col] = 10 * np.log10(Sxx[row][col]) - ((10 * np.log10(MDF[row][col])) + ((10*np.log10(median))))


					else:
						spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))

					# Find the index of the middle frequency
					middle_index = len(times) // 2
					middle_column = spec[:, middle_index]
					vmin = 0  
					vmax = np.max(middle_column) 

					fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     

					ax1.plot(t, data, 'k', linewidth=0.5)
					ax1.set_title(title)
					
					ax1.margins(x=0)

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
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 142
						l = 2450

					if n == 4:
						fnot = [13,27,40,47,54,60,67,74,80,87,90,94,101,108,114,121,127,134,148,160,177,189,202,223,239,247,270]
						tprime0 = 140
						tpr = np.arange(40, 230, 1)
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 67
						l = 580

					if n == 5:
						fnot = [12.5,17.5]
						tprime0 = 123
						tpr = np.arange(105, 140, 1)
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 112
						l = 1150

					if n == 6:
						fnot = [18,36]
						tprime0 = 133
						tpr = np.arange(100, 200, 1)
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 92
						l = 2400

					if n == 7:
						fnot = [26, 49]		
						tprime0 = 122
						tpr = np.arange(50, 200, 1)
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 126
						l = 3000

					if n == 8:
						fnot = [27,57.7,87.7]
						tprime0 = 100
						tpr = np.arange(50, 250, 1)
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 67
						l = 2300

					if n == 9:
						fnot = [26]
						tprime0 = 114
						tpr = np.arange(60, 170, 1)
						tpr = np.arange(0, 241, 1)
						c = 343
						v0 = 144
						l = 1900

					if n == 10:
						fnot = [38, 57, 76, 96, 116, 135, 154, 173, 192, 211, 231]
						tprime0 = 112
						tpr = np.arange(0, 241, 1)
						v0 = 53 
						l = 831

					if n == 11:
						fnot = [19,40,59,79,100,120,140,160,180,200,221,241,261]
						tprime0 = 118
						tpr = np.arange(0, 241, 1)
						v0 = 59
						l = 479

					if n == 12:
						fnot = [14,32,43,48,64,80,86,96,112,129,145,158,161,180,194,202,210,227,243,260,277]
						tprime0 = 110
						tpr = np.arange(0, 241, 1)
						v0 = 89
						l = 1307
					c = 343
					# Plot forward model
					if forward_model == True:
						for f0 in fnot:
							ft = []
							for tprime in tpr:
								#l = np.sqrt(dist_h**2 + alt**2)
								ft0p = f0*1/(1+(v0/c)*(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
									
								ft.append(ft0p)
							ax2.plot(tpr, ft, 'g', linewidth=0.5)
							ax2.set_title("Forward Model: t'= "+str(tprime0)+' sec, v0 = '+str(v0)+' m/s, l = '+str(l)+' m, \n' + 'f0 = '+str(fnot)+' Hz', fontsize='x-small')
					# Plot spectrogram
					cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
					ax2.set_xlabel('Time (s)')
					_,_,dist_km, tmid = closest_encounter(flight_latitudes, flight_longitudes,line, tm, seismo_latitudes[y], seismo_longitudes[y])
					dist_m = dist_km * 1000
					tarrive = tim + (time[n] - calc_time(tmid,dist_m,alt_m))
					tarrive_est = calc_time(tprime0,dist_m,alt_m)
					print(tmid, tarrive)

					ax2.axvline(x=tarrive, c = 'r', ls = '--',label='Wave arrvial: '+str(np.round(tarrive,2))+' s')
					if n < 10:
						ax2.axvline(x=tprime0, c = 'g', ls = '--', label='Estimated arrival: '+str(tprime0)+' s')
					ax2.legend(loc='upper right',fontsize = 'x-small')
					ax2.set_ylabel('Frequency (Hz)')

					
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
					ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2) 
					ax4.plot(middle_column2, frequencies, c='orange')  
					ax4.set_ylim(0, int(fs/2))
					ax4.set_xlim(vmax*1.1, vmin) 
					ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
					ax4.grid(axis='y')
					
					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'
					make_base_dir(BASE_DIR)
					fig.savefig('/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(time[n])+'_'+str(flight_num[n])+'.png')
					plt.close()


					closest_index = np.argmin(np.abs(tarrive - times))
					arrive_time = spec[:,closest_index]
					for i in range(len(arrive_time)):
						if arrive_time[i] < 0:
							arrive_time[i] = 0
					vmin = np.min(arrive_time) 
					vmax = np.max(arrive_time) 
					peaks, _ = signal.find_peaks(arrive_time, prominence=10, distance = 10, height = 5, width=1) #for later change parameters for jets and permenant stations
					np.diff(peaks)
					fig = plt.figure(figsize=(10,6))
					plt.grid()
					
					plt.plot(frequencies, arrive_time, c='c')
					plt.plot(peaks, arrive_time[peaks], "x")
					for g in range(len(peaks)):
						plt.text(peaks[g], arrive_time[peaks[g]], peaks[g], fontsize=15)


					plt.xlim(0, int(fs/2))
					plt.xticks(fontsize=12)
					plt.yticks(fontsize=12)
					plt.ylim(0,vmax*1.1)
					plt.xlabel('Frequency (Hz)', fontsize=17)
					plt.ylabel('Amplitude Spectrum at t = {:.2f} s (dB)'.format(tarrive), fontsize=17)
					
					
					make_base_dir('/scratch/irseppi/nodal_data/plane_info/5spec/20190'+str(month[n])+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/')
					
					fig.savefig('/scratch/irseppi/nodal_data/plane_info/5spec/20190'+str(month[n])+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(sta[n])+'_' + str(time[n]) + '.png')
					plt.close()

