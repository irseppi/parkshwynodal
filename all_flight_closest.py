import pandas as pd
from obspy.geodetics import gps2dist_azimuth
import numpy as np
from prelude import load_flights
import pandas as pd
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from prelude import load_flights

# Function to calculate distance between two coordinates
def distance(lat1, lon1, lat2, lon2):
	dist,_,_ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist/1000
	return dist_km

# Load flight files and filenames
flight_files, filenames = load_flights(2, 4, 11, 27)

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_sta = seismo_data['Station']

# Open output file for writing
output = open('all_station_crossing_db_updated.txt', 'w')

# Iterate over flight files
for i, flight_file in enumerate(flight_files):
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
	timestamp = time

	# Create a dictionary to store flight data for faster lookup
	flight_dict = {}
	for l in range(len(flight_data)):
		if l == 0 or l >= len(flight_data) - 1:
				continue

		flight_dict[l] = {
			'clat': flight_latitudes[l],
			'clon': flight_longitudes[l],
			'c2lat': flight_latitudes[l + 1],
			'c2lon': flight_longitudes[l + 1],
			'sclat': flight_latitudes[l - 1],
			'sclon': flight_longitudes[l - 1]
		}

	# Iterate over seismometer data
	for s in range(len(seismo_data)):
		dist_lim = 2.01
		seismo_lat = seismo_latitudes[s]
		seismo_lon = seismo_longitudes[s]

		# Iterate over flight data
		for l in range(len(flight_data)):
			if l == 0 or l >= len(flight_data) - 1:
				continue
			flight_info = flight_dict[l]
			clat = flight_info['clat']
			clon = flight_info['clon']
			for gg in range(0,2):
				if gg == 0:
					sclat = flight_info['sclat']
					sclon = flight_info['sclon']
					l2 = l + 1
				else:
					sclat = flight_info['c2lat']
					sclon = flight_info['c2lon']
					l2 = l - 1
            
			x = [clon, sclon]
			y = [clat, sclat]

			if clon < sclon:
				m = (y[0] - y[1]) / (x[0] - x[1])
				b = y[0] - m * x[0]
				lon_range = np.arange(clon, sclon, 0.000001)
			else:
				m = (y[1] - y[0]) / (x[1] - x[0])
				b = y[0] - m * x[0]
				lon_range = np.arange(sclon, clon, 0.000001)

			# Calculate distances for each longitude value
			lat_range = m * lon_range + b
			dist_km_range = []
			for i in range(len(lon_range)):
				dist_km_range.append(distance(lat_range[i], lon_range[i], seismo_lat, seismo_lon))

			# Find the minimum distance within the range
			min_dist_idx = np.argmin(dist_km_range)
			min_dist = dist_km_range[min_dist_idx]
			closest_lat = lat_range[min_dist_idx]
			closest_lon = lon_range[min_dist_idx]

			if min_dist < dist_lim:
				dist_lim = min_dist
				c2lat = sclat
				c2lon = sclon
				l2 = l2

		if dist_lim < 2.01:
			# Iterate over longitude values within a small range
			lon_range = np.arange(closest_lon - 0.000001, closest_lon + 0.000001, 0.0000000001)
			lat_range = m * lon_range + b
			for i in range(len(lon_range)):
				dist_km_range.append(distance(lat_range[i], lon_range[i], seismo_lat, seismo_lon))

			# Find the minimum distance within the range
			min_dist_idx = np.argmin(dist_km_range)
			min_dist = dist_km_range[min_dist_idx]
			closest_lat = lat_range[min_dist_idx]
			closest_lon = lon_range[min_dist_idx]

			if min_dist < 1:
				# Iterate over longitude values within a smaller range
				lon_range = np.arange(closest_lon - 0.0000000001, closest_lon + 0.0000000001, 0.0000000000001)
				lat_range = m * lon_range + b
				for i in range(len(lon_range)):
					dist_km_range.append(distance(lat_range[i], lon_range[i], seismo_lat, seismo_lon))

				# Find the minimum distance within the range
				min_dist_idx = np.argmin(dist_km_range)
				min_dist = dist_km_range[min_dist_idx]
				closest_lat = lat_range[min_dist_idx]
				closest_lon = lon_range[min_dist_idx]

			dist_old_new = distance(closest_lat, closest_lon, clat, clon)
			dist_old_old = distance(c2lat, c2lon, clat, clon)
			ratio = dist_old_new / dist_old_old
			timestamp = timestamp[l] + (timestamp[l] - timestamp[l2]) * ratio

			timeF = timestamp
			altF = (alt[l2] + alt[l]) / 2
			speedF = (speed[l2] + speed[l]) / 2
			seismo_staF = seismo_sta[s]

			with open('/scratch/irseppi/nodal_data/flightradar24/' + date + '_flights.csv') as text:
				# Iterate over lines in the text file
				for line in text:
					val = line.split(',')
					if val[0] == flight_num:
						# Write data to the output file
						output.write(
							f"{date},{flight_num},{val[1]},{timeF},{altF},{speedF},{seismo_staF},{val[3]},\n")
						break

	print((i / len(flight_files)) * 100, '%')

output.close()
