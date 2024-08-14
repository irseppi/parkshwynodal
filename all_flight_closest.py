import pandas as pd
from obspy.geodetics import gps2dist_azimuth
import numpy as np
from prelude import load_flights

# Function to calculate distance between two coordinates
def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km

# Load flight files and filenames
flight_files = load_flights(2, 4, 11, 27)[0]
filenames = load_flights(2, 4, 11, 27)[1]

# Open output file for writing
output = open('all_station_crossing_db_updated.txt','w')

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_sta = seismo_data['Station']

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

	# Iterate over seismometer data
	for s in range(len(seismo_data)):
		dist_lim = 2.01

		# Iterate over flight data
		for l in range(len(flight_data)):
			if l == 0:
				continue
			elif l >= len(flight_data)-1:
				continue
			else:
				pass
			clat = flight_latitudes[l]
			clon = flight_longitudes[l]
			
			closest_lat = clat
			closest_lon = clon
			c2lat = flight_latitudes[l+1]
			c2lon = flight_longitudes[l+1]
			l2 = l + 1
			# Iterate over two points for linear interpolation
			for tr in range(0,2):
				if tr == 0:
					sclat = flight_latitudes[l-1]
					sclon = flight_longitudes[l-1]  
				else:
					sclat = flight_latitudes[l+1]
					sclon = flight_longitudes[l+1]

				x = [clon, sclon]
				y = [clat, sclat]
				
				if clon < sclon:
					m = (y[0]-y[1])/(x[0]-x[1])
					b = y[0] - m*x[0]
					# Iterate over longitude values between two points
					for point in np.arange(clon, sclon, 0.000001):
						lon = point
						lat = m*lon + b
						dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])
						if dist_km < dist_lim:
							dist_lim = dist_km
							closest_lat = lat
							closest_lon = lon
							c2lat = flight_latitudes[l-1]
							c2lon = flight_longitudes[l-1]
							l2 = l - 1
						else:
							continue
				else:
					m = (y[1]-y[0])/(x[1]-x[0])
					b = y[0] - m*x[0]
					# Iterate over longitude values between two points
					for point in np.arange(sclon, clon, 0.000001):
						lon = point
						lat = m*lon + b
						dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])

						if dist_km < dist_lim:
							dist_lim = dist_km
							closest_lat = lat
							closest_lon = lon
							c2lat = flight_latitudes[l+1]
							c2lon = flight_longitudes[l+1]
							l2 = l + 1
						else:
							continue
			print(dist_lim)
		if dist_lim < 2.01:
			print('here')
			# Iterate over longitude values within a small range
			for location in np.arange((closest_lon-0.000001),(closest_lon+0.000001),0.0000000001):
				lon = location
				lat = m*lon + b
				dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])
				if dist_km < dist_lim:
					dist_lim = dist_km
					closest_lon = lon
					closest_lat = lat
				else:
					continue
			if dist_lim < 1:
				# Iterate over longitude values within a smaller range
				for location in np.arange((closest_lon-0.0000000001),(closest_lon+0.0000000001),0.0000000000001):
					lon = location
					lat = m*lon + b
					dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])
					if dist_km < dist_lim:
						dist_lim = dist_km
						closest_lon = lon
						closest_lat = lat
					else:
						continue
			else:
				pass
			dist_old_new = distance(closest_lat, closest_lon, clat, clon)
			dist_old_old = distance(c2lat, c2lon, clat, clon)
			ratio = dist_old_new/dist_old_old
			timestamp = timestamp[l] + (timestamp[l] - timestamp[l2])*ratio

			timeF = timestamp
			altF = (alt[l2]+alt[l])/2
			speedF = (speed[l2]+speed[l])/2
			seismo_staF = seismo_sta[s]



			text = open('/scratch/irseppi/nodal_data/flightradar24/'+date+'_flights.csv')
			# Iterate over lines in the text file
			for line in text.readlines():
				val = line.split(',')
				if val[0] == flight_num:
					# Write data to the output file
					output.write(str(date)+','+ str(flight_num)+','+val[1]+','+str(timeF)+','+str(altF)+','+str(speedF)+','+str(seismo_staF)+','+val[3]+',\n')
					break
				else:
					continue
		else:
			continue
	print((i/len(flight_files))*100, '%')
output.close()

