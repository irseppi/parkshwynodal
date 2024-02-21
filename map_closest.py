import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib.patches as mpatch
from matplotlib.patches import Rectangle
import os
from obspy.geodetics import gps2dist_azimuth
from datetime import datetime
from pathlib import Path
from math import radians, sin, cos, sqrt, atan2
from prelude import make_base_dir, dist_less, load_flights, distance
from datetime import datetime, timedelta

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
			idistance = gps2dist_azimuth(flight_latitudes[n], flight_longitudes[n], seismo_latitudes[t], seismo_longitudes[t])[0]

			if float(idistance) < float(closest_distance):
				closest_distance = idistance
				closest_time = timestamp[n]
				closest_altitude = altitude[n] * 0.3048
				closest_speed = speed[n] * 0.514
				closest_head = head[n]
				closest_lat = flight_latitudes[n] 
				closest_lon = flight_longitudes[n]

				closest_seislat = seismo_latitudes[t] 
				closest_seislon = seismo_longitudes[t]
				closest_station = stations[t]
				try:
					if gps2dist_azimuth(flight_latitudes[n-1], flight_longitudes[n-1], seismo_latitudes[t], seismo_longitudes[t])[0] < gps2dist_azimuth(flight_latitudes[n+1], flight_longitudes[n+1], seismo_latitudes[t], seismo_longitudes[t])[0]:
						closes_distance2 = gps2dist_azimuth(flight_latitudes[n-1], flight_longitudes[n-1], seismo_latitudes[t], seismo_longitudes[t])[0]
						closest_time2 = timestamp[n-1]
						closest_altitude2 = altitude[n-1] * 0.3048
						closest_speed2 = speed[n-1] * 0.514
						closest_head2 = head[n-1]
						closest_lat2 = flight_latitudes[n-1] 
						closest_lon2 = flight_longitudes[n-1]
					elif gps2dist_azimuth(flight_latitudes[n-1], flight_longitudes[n-1], seismo_latitudes[t], seismo_longitudes[t])[0] > gps2dist_azimuth(flight_latitudes[n+1], flight_longitudes[n+1], seismo_latitudes[t], seismo_longitudes[t])[0]:
						closes_distance2 = gps2dist_azimuth(flight_latitudes[n+1], flight_longitudes[n+1], seismo_latitudes[t], seismo_longitudes[t])[0]
						closest_time2 = timestamp[n+1]
						closest_altitude2 = altitude[n+1] * 0.3048
						closest_speed2 = speed[n+1] * 0.514
						closest_head2 = head[n+1]
						closest_lat2 = flight_latitudes[n+1] 
						closest_lon2 = flight_longitudes[n+1]
				except:
						average_speed = closest_speed
						heading_direction = closest_head
						time_of_closest_approach = closest_time
						avg_alt = closest_altitude
						closest_lat = closest_lat
						closest_lon = closest_lon
			
	if closest_time2 is not None:
		line_vector = (closest_lat2 - closest_lat, closest_lon2 - closest_lon)
		station_vector = (closest_seislat - closest_lat, closest_seislon - closest_lon)
		# Calculate the projection of the station vector onto the line vector
		dot_product = line_vector[0]*station_vector[0] + line_vector[1]*station_vector[1]
		line_magnitude_squared = line_vector[0]**2 + line_vector[1]**2
		projection_length_ratio = dot_product / line_magnitude_squared

		# Calculate the coordinates of the closest point on the line to the station
		closest_point_on_line = (closest_lat + projection_length_ratio * line_vector[0], closest_lon + projection_length_ratio * line_vector[1])

		# Calculate the distance from the closest point on the line to the station
		closest_distance = gps2dist_azimuth(closest_point_on_line[0],closest_point_on_line[1], closest_seislat, closest_seislon)[0]

		# Calculate the average speed and heading direction between the two closest points
		time_difference = (datetime.fromtimestamp(closest_time2) - datetime.fromtimestamp(closest_time)).total_seconds()
		distance = gps2dist_azimuth(closest_lat, closest_lon, closest_lat2, closest_lon2)[0]
		average_speed = float(distance) / float(time_difference)
		#print(average_speed)
		average_speed = (closest_speed + closest_speed2)/2
	
		heading_direction = (closest_head + closest_head2) / 2
		avg_alt = (closest_altitude + closest_altitude2) / 2

		# Calculate the time of the closest approach
		time_of_closest_approach = datetime.fromtimestamp(closest_time) + timedelta(seconds=time_difference*projection_length_ratio)

	return time_of_closest_approach, closest_distance, average_speed, heading_direction, avg_alt, closest_seislat, closest_seislon, closest_station, closest_point_on_line


def calculate_wave_arrival(closest_time, closest_distance, closest_altitude, aircraft_speed, aircraft_heading, seismo_latitudes, seismo_longitudes):
	speed_of_sound = 343  # speed of sound in m/s at sea level
	aircraft_speed_mps = aircraft_speed 

	# calculate the component of the aircraft's speed in the direction of the seismometer
	relative_heading = radians(aircraft_heading)
	relative_speed = aircraft_speed_mps * cos(relative_heading)

	# calculate the time it takes for the sound to travel from the aircraft to the seismometer
	time_for_sound_to_travel = np.sqrt((closest_distance)**2 + (closest_altitude)**2) / (speed_of_sound + relative_speed)

	wave_arrival_time = float(closest_time + time_for_sound_to_travel)
	print(closest_time, wave_arrival_time)
	return wave_arrival_time


# Load the seismometer location data
seismo_data = pd.read_csv('input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3
f = 1.0/np.cos(62*np.pi/180)

file_names = ['/scratch/irseppi/nodal_data/flightradar24/20190225_positions/20190225_530342801.csv','/scratch/irseppi/nodal_data/flightradar24/20190214_positions/20190214_528485724.csv','/scratch/irseppi/nodal_data/flightradar24/20190214_positions/20190214_528473220.csv','/scratch/irseppi/nodal_data/flightradar24/20190214_positions/20190214_528407493.csv','/scratch/irseppi/nodal_data/flightradar24/20190213_positions/20190213_528293430.csv']

flight_num = [530342801,528485724,528473220,528407493,528293430]
time = [1551066051,1550172833,1550168070,1550165577,1550089044]
sta = [1022,1272,1173,1283,1004]
day = [25,14,14,14,13]

for i in range(0,6):
	flight_data = pd.read_csv(file_names[i], sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	timestamp = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']
	head = flight_data['heading']

	if dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes) == True:
		ctime, cdist, cspeed, chead, calt,  cseislat, cseislon, csta, cpoint = closest_encounter(flight_latitudes, flight_longitudes, timestamp, alt, speed, head, seismo_latitudes, seismo_longitudes, stations)	
		#tm = calculate_wave_arrival(ctime, cdist, calt, cspeed, chead, cseislat, cseislon)

		dist = cdist
		ht = ctime #datetime.utcfromtimestamp(ctime)
		# Create a figure with two subplots side by side
		fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

		axs[0].set_aspect(f)
		axs[1].set_aspect(f)

		axs[0].scatter(seismo_longitudes, seismo_latitudes, c='red', s = 3, label='seismometers')
		axs[0].plot(flight_longitudes, flight_latitudes, '-o', c='c', lw=1, ms = 1, label='flight path')

		# Set labels and title
		axs[0].set_xlim(min_lon, max_lon)
		axs[0].set_ylim(min_lat, max_lat)
		axs[0].tick_params(axis='both', which='major', labelsize=13)


		y =[cpoint[0],  cseislat]
		x = [cpoint[1], cseislon]
		yy = sum(y)/len(y)
		xx = sum(x)/len(x)

		minl = xx - 0.025
		maxl = xx + 0.025
		minla = yy - 0.015
		maxla = yy + 0.015
		direction = np.deg2rad(chead)

		rect = Rectangle((minl, minla), 0.1, 0.06, ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
		axs[0].add_patch(rect)
		axs[1].plot(x,y, '--', c='orange')
		# Draw the zoomed in map on the second subplot
		axs[1].scatter(seismo_longitudes, seismo_latitudes, c='red')
		#axs[1].quiver(flight_longitudes[l], flight_latitudes[l], np.cos(direction), np.sin(direction), angles='xy') #, scale_units='xy', scale=0.002)
		axs[1].plot(flight_longitudes, flight_latitudes, c='c',linestyle ='dotted')
		axs[1].set_xlim(minl, maxl)
		axs[1].set_ylim(minla, maxla)
		axs[1].text(cseislon, cseislat, csta, fontsize=11, fontweight='bold')
		axs[1].tick_params(axis='both', which='major', labelsize=13)
		axs[1].text(cpoint[1], cpoint[0], ht, fontsize=11, fontweight='bold')
		axs[1].scatter(cpoint[1], cpoint[0], c='lawngreen')
		axs[1].scatter(cseislon, cseislat, c='pink')


		axs[1].text(xx,yy, str(round(dist, 2)/1000)+'km', fontsize=10, fontweight='bold')

		# Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
		con = mpatch.ConnectionPatch(xyA=(minl, minla), xyB=(maxl, minla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
		fig.add_artist(con)
		con = mpatch.ConnectionPatch(xyA=(minl, maxla), xyB=(maxl, maxla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
		fig.add_artist(con)
		      
		BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/map_5/'
		make_base_dir(BASE_DIR)
		plt.savefig('/scratch/irseppi/nodal_data/plane_info/map_5/'+str(flight_num[i])+'.png')
		plt.close()
	else:
		continue


