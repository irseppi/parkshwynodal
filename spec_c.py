import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import scipy.io
import os
import obspy
from scipy.signal import spectrogram
from scipy import signal
from obspy.geodetics import gps2dist_azimuth
from datetime import datetime
from pathlib import Path
from math import radians, sin, cos, sqrt, atan2
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


# Load the seismometer location data
seismo_data = pd.read_csv('input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']

flight_files = load_flights(2,3,15,20)[0]				
filenames = load_flights(2,3,15,20)[1]
for i, flight_file in enumerate(flight_files):
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	timestamp = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']
	head = flight_data['heading']

	if dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes) == True:
		ctime, cdist, cspeed, chead, calt,  cseislat, cseislon, csta = closest_encounter(flight_latitudes, flight_longitudes, timestamp, alt, speed, head, seismo_latitudes, seismo_longitudes, stations)	
		tm = calculate_wave_arrival(ctime, cdist, calt, cspeed, chead, cseislat, cseislon)

		ht = datetime.utcfromtimestamp(tm)
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

		
		n = "/scratch/naalexeev/NODAL/2019-0"+ str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month2)+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(csta)+".mseed"

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
			ax1.axvline(x=center_time, c = 'r', ls = '--')
			ax1.margins(x=0)
			
			spec = (10 * np.log10(Sxx)) - (10 * np.log10(MDF))

			# Find the index of the middle frequency
			middle_index = len(times) // 2

			# Extract the middle line of the spectrogram
			middle_column = spec[:, middle_index]

			vmin = 0
			vmax = 1
			if np.min(middle_column) < 0:
				vmax = np.abs(np.min(middle_column))+np.max(middle_column)
				spec = np.abs(np.min(middle_column)) + np.array(spec)
				middle_column = np.abs(np.min(middle_column)) + np.array(middle_column)
			else:
				vmin = 0
				vmax = np.max(middle_column)
			# Plot spectrogram
			cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='hot_r', vmin=vmin, vmax=vmax) #'hsv' w/ vmin = np.min(middle_column) 
			ax2.set_xlabel('Time [s]')



			ax2.axvline(x=center_time, c = 'k', ls = '--')
			ax2.set_ylabel('Frequency (Hz)')
			ax2.margins(x=0)
			ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

			plt.colorbar(mappable=cax, cax=ax3)
			ax3.set_ylabel('Relative Amplitude (dB)')

			make_base_dir('/scratch/irseppi/nodal_data/Plane_map_spec/')
			fig.savefig('/scratch/irseppi/nodal_data/Plane_map_spec/spec_'+str(tm)+'_'+str(csta)+'_'+filenames[i]+'.png')
			
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

			make_base_dir('/scratch/irseppi/nodal_data/Plane_map_spec/')
			fig.savefig('/scratch/irseppi/nodal_data/Plane_map_spec/fft_'+str(tm)+'_'+str(csta)+'_'+filenames[i]+'.png')
			plt.close()
	else:
		continue
print(i/len(flight_files), '% Done')	
