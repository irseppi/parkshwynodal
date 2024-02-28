import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy import signal
from obspy.geodetics import gps2dist_azimuth
import obspy
import math
import scipy.io
import os
from obspy.core import UTCDateTime
import datetime
from prelude import make_base_dir, dist_less, load_flights, distance


def closest_encounter(flight_latitudes, flight_longitudes, timestamp, altitude, speed, head, seismo_latitudes, seismo_longitudes, stations):
	closest_distance = float('inf')
	closest_time = 0
	closest_altitude = None
	closest_speed = None
	closest_head = None
	closest_lat = 0 
	closest_lon = 0
	closest_seislat = 0
	closest_seislon = 0
	closest_station = 0
	closest_distance2 = float('inf')
	closest_time2 = None
	closest_altitude2 = None
	closest_speed2 = None
	closest_head2 = None
	closest_lat2 = 0 
	closest_lon2 = 0

	for n in range(len(flight_latitudes)):
		for t in range(len(seismo_latitudes)):
			_, _, distance = gps2dist_azimuth(flight_latitudes[n], flight_longitudes[n], seismo_latitudes[t], seismo_longitudes[t])

			if float(distance) < float(closest_distance):
				if int(stations[t]) == int(closest_station):
					if closest_distance2 < closest_distance:
						closest_distance2 = closest_distance
						closest_time2 = closest_time
						closest_altitude2 = closest_altitude
						closest_speed2 = closest_speed
						closest_head2 = closest_head
						closest_lat2 = closest_lat
						closest_lon2 = closest_lon


				closest_distance = distance
				closest_time = timestamp[n]
				closest_altitude = altitude[n] * 0.3048
				closest_speed = speed[n] * 0.514
				closest_head = head[n]
				closest_lat = flight_latitudes[n] 
				closest_lon = flight_longitudes[n]

				closest_seislat = seismo_latitudes[t] 
				closest_seislon = seismo_longitudes[t]
				closest_station = stations[t]

			elif int(stations[t]) == int(closest_station):
				if float(distance) < float(closest_distance2):
					closes_distance2 = distance
					closest_time2 = timestamp[n]
					closest_altitude2 = altitude[n] * 0.3048
					closest_speed2 = speed[n] * 0.514
					closest_head2 = head[n]
					closest_lat2 = flight_latitudes[n] 
					closest_lon2 = flight_longitudes[n]
	if float(closest_distance2) == float('inf'):
		average_speed = closest_speed
		heading_direction = closest_head
		time_of_closest_approach = closest_time
		avg_alt = closest_altitude
		closest_seis = closest_seislat
		closest_seislon = closest_seislon
	else:
		line_vector = (closest_lat2 - closest_lat, closest_lon2 - closest_lon)
		station_vector = (closest_seislat - closest_lat, closest_seislon - closest_lon)

		# Calculate the projection of the station vector onto the line vector
		dot_product = line_vector[0]*station_vector[0] + line_vector[1]*station_vector[1]
		line_magnitude_squared = line_vector[0]**2 + line_vector[1]**2
		projection_length_ratio = dot_product / line_magnitude_squared

		# Calculate the coordinates of the closest point on the line to the station
		closest_point_on_line = (closest_lat + projection_length_ratio * line_vector[0], closest_lon + projection_length_ratio * line_vector[1])

		# Calculate the distance from the closest point on the line to the station
		closest_distance = distance(closest_point_on_line[0],closest_point_on_line[1], closest_seislat, closest_seislon)*1000

		# Calculate the average speed and heading direction between the two closest points
		time_difference = (closest_time2 - closest_time).total_seconds()
		distance = distance(closest_lat, closest_lon, closest_lat2, closest_lon2)*1000
		average_speed = distance / time_difference
		#print(average_speed)
		#average_speed = (closest_speed + closest_speed2)/2
		print(average_speed)
		heading_direction = (closest_head + closest_head2) / 2
		avg_alt = (closest_altitude + closest_altitude2) / 2

		# Calculate the time of the closest approach
		time_of_closest_approach = datetime.fromtimestamp(closest_time) + datetime.timedelta(seconds=time_difference*projection_length_ratio)
	
	return time_of_closest_approach, closest_distance, average_speed, heading_direction, avg_alt, closest_seislat, closest_seislon, closest_station


def calculate_wave_arrival(closest_time, closest_distance, closest_altitude, aircraft_speed, aircraft_heading, seismo_latitudes, seismo_longitudes):
	speed_of_sound = 343  # speed of sound in m/s at sea level
	aircraft_speed_mps = aircraft_speed 

	# calculate the component of the aircraft's speed in the direction of the seismometer
	#relative_heading = radians(aircraft_heading)
	#relative_speed = aircraft_speed_mps * cos(relative_heading)

	# calculate the time it takes for the sound to travel from the aircraft to the seismometer
	time_for_sound_to_travel = np.sqrt((closest_distance)**2 + (closest_altitude)**2) / (speed_of_sound) # + relative_speed)

	wave_arrival_time = float(closest_time + time_for_sound_to_travel)
	print(closest_time, wave_arrival_time)
	return wave_arrival_time


def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km
def calc_time(t, l, vo):
	#t is epoch time at time wave is generated by aircraft sec
	#l is the shortest distance to between the station and aircraft km
	#vo is the aircraft velocity km/sec
	#c is the speed of aucostic wave propogation km/sec
	c = 0.343 #km/sec
	to=t+(math.sqrt(l**2+(vo*t)**2))/c
	return to
#20190225_530342801_1551066051_1022
#20190214_528485724_1550172833_1272
#20190214_528473220_1550168070_1173
#20190214_528407493_1550165577_1283
#20190213_528293430_1550089044_1004
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
					#spec = 10 * np.log10(Sxx)
					#spec = 10 * np.log10(MDF)

					middle_index = len(times) // 2
					middle_column = spec[:, middle_index]
					vmin = 0  #np.min(middle_column)
					vmax = np.max(middle_column) # - np.min(middle_column)
					# Plot spectrogram
					cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='hot_r', vmin=vmin, vmax=vmax)
					ax2.set_xlabel('Time [s]')
					# Find the center of the trace
					center_index = len(data) // 2
					center_time = t[center_index]
					
					l = np.sqrt(dist**2+(alt*0.0003048)**2) #km
					vo = speed*0.000514444
					xloc = calc_time(center_time, l, vo)
					ax2.axvline(x=center_time, c = 'r', ls = '--')
					ax2.axvline(x=xloc, c = 'k', ls = '--')
					ax2.set_ylabel('Frequency (Hz)')
					ax2.margins(x=0)
					ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])
					plt.colorbar(mappable=cax, cax=ax3)
					ax3.set_ylabel('Relative Amplitude (dB)')
					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-02-'+str(day[n])+ '/'+str(flight_num[n])+'/'+station[y]+'/'
					make_base_dir(BASE_DIR)
					fig.savefig('/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-02-'+str(day[n])+ '/'+str(flight_num[n])+'/'+station[y]+'/'+str(time[n])+'_'+str(flight_num[n])+'.png')
					plt.close()
					# Find the index of the middle frequency
					middle_index = len(times) // 2
					# Extract the middle line of the spectrogram
					middle_column = spec[:, middle_index]
					
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
					#plt.ylim(vmin,vmax)
					plt.xlabel('Freq [Hz]')
					plt.ylabel('Amplitude [dB]')
					plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))
					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5spec/201902'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'/'
					make_base_dir(BASE_DIR)
					fig.savefig('/scratch/irseppi/nodal_data/plane_info/5spec/201902'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'/'+station[y]+'_'+str(time[n])+'.png')
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
					
					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5spec2/201902'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'/'
					make_base_dir(BASE_DIR)
					fig.savefig('/scratch/irseppi/nodal_data/plane_info/5spec2/201902'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'/'+station[y]+'_'+str(time[n])+'.png', bbox_inches='tight',pad_inches =0)
					plt.close()