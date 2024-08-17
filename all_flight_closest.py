import pandas as pd
from obspy.geodetics import gps2dist_azimuth
import numpy as np
from prelude import load_flights, closest_encounter, distance

# Load flight files and filenames
flight_files,filenames = load_flights(2, 4, 11, 27)

# Open output file for writing
output = open('all_station_crossing_db_updated111.txt','w')

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_sta = seismo_data['Station']

# Iterate over flight files
for i, flight_file in enumerate(flight_files):
	print((i/len(flight_files))*100, '%')	
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	alt = flight_data['altitude']
	speed = flight_data['speed']
	head = flight_data['heading']
	fname = filenames[i]	
	flight_num = fname[9:18]
	date = fname[0:8]

	# Iterate over seismometer data
	for s in range(len(seismo_data)):
		dist_lim = 2.01
		con = False
		# Iterate over flight data
		for l in range(len(flight_data)):
			if l == 0:
				continue
			elif l >= len(flight_data)-1:
				continue
			else:
				dist = distance(flight_latitudes[l], flight_longitudes[l], seismo_latitudes[s], seismo_longitudes[s])
				if dist < 2.01:
					clat, clon, dist, ctime = closest_encounter(flight_latitudes, flight_longitudes, l, time, seismo_latitudes[s], seismo_longitudes[s])
					if clat != None:
						# Check if the distance is less than the current minimum distance
						if dist <= dist_lim:
							timeF = ctime
							altF = alt[l]
							speedF = speed[l]
							seismo_staF = seismo_sta[s]
							headF = head[l]
							con = True
							closest_lat = clat
							closest_lon = clon
							dist_lim = dist

						else:
							continue
					else:
						continue
				else:
					continue
		if con == True:
			# Write data to the output file
			output.write(str(date)+','+ str(flight_num)+','+str(clat)+','+str(clon)+','+str(timeF)+','+str(altF)+','+str(speedF)+','+str(seismo_staF)+','+str(headF)+',\n')
			print('Flight:', flight_num, 'Station:', seismo_staF, 'Distance:', dist_lim)
		else:
			continue

output.close()

