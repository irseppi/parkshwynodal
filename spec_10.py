import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
import obspy
import datetime
from prelude import make_base_dir, distance, closest_encounter, calc_time

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station = seismo_data['Station']
flight_num = [530342801,528485724,528473220,528407493,528293430,527937367,529741194,529776675,529179112,530165646]
time = [1551066051,1550172833,1550168070,1550165577,1550089044,1549912188,1550773710,1550787637,1550511447,1550974151]
sta = [1022,1272,1173,1283,1004,"CCB","F6TP","F4TN","F3TN","F7TV"]
day = [25,14,14,14,13,11,21,21,18,24]

for n in range(0,10):
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
						p = "/scratch/naalexeev/NODAL/2019-02-"+str(day[n])+"T"+str(h)+":00:00.000000Z.2019-02-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(station[y])+".mseed"
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
					ax1.margins(x=0)

					spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))

					# Find the index of the middle frequency
					middle_index = len(times) // 2
					middle_column = spec[:, middle_index]
					vmin = 0  
					vmax = np.max(middle_column)
		
					if n == 0:
						tprime0 = 112
						
					if n == 1:
						tprime0 = 107
						
					if n == 2:
						tprime0 = 93
					if n == 3:
						tprime0 = 116
						
					if n == 4:
						tprime0 = 140

					if n == 5:
						tprime0 = 123

					if n == 6:
						tprime0 = 133

					if n == 7:	
						tprime0 = 122

					if n == 8:
						tprime0 = 100
					if n == 9:
						tprime0 = 114

					# Plot spectrogram
					cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
					ax2.set_xlabel('Time [s]')
					dist_m, tmid = closest_encounter(flight_latitudes, flight_longitudes,line, tm, seismo_latitudes[y], seismo_longitudes[y])
					tarrive = tim + (time[n] - calc_time(tmid,dist_m,alt_m))
					tarrive_est = calc_time(tprime0,dist_m,alt_m)
					print(tmid, tarrive)

					ax2.axvline(x=tarrive, c = 'r', ls = '--',label='Wave arrvial: '+str(np.round(tarrive,2))+'sec')
					ax2.axvline(x=tprime0, c = 'g', ls = '--', label='Estimated arrival: '+str(tprime0)+' sec')
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
					ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2) #, width=vmax*1.1-vmin, height=int(fs/2))
					ax4.plot(middle_column2, frequencies, c='orange')  
					ax4.set_ylim(0, int(fs/2))
					ax4.set_xlim(vmax*1.1, vmin) #, width=vmax*1.1-vmin, height=int(fs/2))
					ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
					ax4.grid(axis='y')
					plt.show()
					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-02-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'
					make_base_dir(BASE_DIR)
					fig.savefig('/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-02-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(time[n])+'_'+str(flight_num[n])+'.png')
					plt.close()