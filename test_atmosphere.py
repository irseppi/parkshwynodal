from datetime import datetime, timezone
import os
from pyproj import Proj
from prelude import *
import pandas as pd

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

def find_closest_point(flight_utm, seismo_utm):
	"""
	Find the closest point on a flight path to a seismic station.

	Args:
		flight_utm (list): List of UTM coordinates of the flight path.
		seismo_utm (tuple): UTM coordinates of the seismic station.

	Returns:
		tuple: A tuple containing the closest point on the flight path, the distance between the flight path and the station, and the index of the closest point.
	"""
	min_distance = np.inf
	closest_point = None
	index = None
	min_station = None
	for i in range(len(flight_utm)-1):
		for g in range(len(seismo_utm)):
			flight_utm_x1, flight_utm_y1 = flight_utm[i]
			flight_utm_x2, flight_utm_y2 = flight_utm[i + 1]
			seismo_utm_x = seismo_utm[g][0]
			seismo_utm_y = seismo_utm[g][1]
			point, d = closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, seismo_utm_x, seismo_utm_y)
		
			if point == None:
				continue
			elif d < min_distance:
				min_distance = d
				closest_point = point
				index = i
				min_station = (seismo_utm_x, seismo_utm_y)
			else:
				continue

	return closest_point, min_distance, index, min_station	
# Load the seismometer location data
seismo_data = pd.read_csv('input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_stations = seismo_data['Latitude']
sta = seismo_data['Station']
elevation = seismo_data['Elevation']

seismo_utm = [utm_proj(lon, lat) for lat, lon in zip(seismo_latitudes, seismo_longitudes)]
seismo_utm_x, seismo_utm_y = zip(*seismo_utm)

# Convert UTM coordinates to kilometers
seismo_utm_x_km = [x / 1000 for x in seismo_utm_x]
seismo_utm_y_km = [y / 1000 for y in seismo_utm_y]

flight_files, filenames = load_flights(2,3,12,13)
for i, flight_file in enumerate(flight_files):
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	timestamps = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']
	utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')
	for t in timestamps:
		time_UTC = datetime.fromtimestamp(t, timezone.utc)
	# Convert flight latitude and longitude to UTM coordinates
	flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
	flight_utm_x, flight_utm_y = zip(*flight_utm)

	# Convert UTM coordinates to kilometers
	flight_utm_x_km = [x / 1000 for x in flight_utm_x]
	flight_utm_y_km = [y / 1000 for y in flight_utm_y]
	flight_path = [(x,y) for x, y in zip(flight_utm_x_km, flight_utm_y_km)]

	seismometer = (seismo_utm_x_km, seismo_utm_y_km)  
	closest_p, dist_km, index, min_sta = find_closest_point(flight_path, seismometer)
	
	if dist_km <= 2:
		closest_x, closest_y = closest_p
		station_x, station_y = min_sta
		#Calculate the time of the closest point
		flight_utm_x1, flight_utm_y1 = flight_path[index]
		flight_utm_x2, flight_utm_y2 = flight_path[index + 1]

		x_timestamp_dif_vec = flight_utm_x2 - flight_utm_x1
		y_timestamp_dif_vec = flight_utm_y2 - flight_utm_y1

		cx_timestamp_dif_vec =  closest_x - flight_utm_x1
		cy_timestamp_dif_vec = closest_y - flight_utm_y1

		line_vector = (x_timestamp_dif_vec, y_timestamp_dif_vec)
		cline_vector = (cx_timestamp_dif_vec, cy_timestamp_dif_vec)

		line_magnitude = np.sqrt(line_vector[0] ** 2 + line_vector[1] ** 2)
		cline_magnitude = np.sqrt(cline_vector[0] ** 2 + cline_vector[1] ** 2)

		length_ratio = cline_magnitude / line_magnitude
		closest_time = timestamps[index] + length_ratio*(timestamps[index+1] - timestamps[index])
		time_UTC = datetime.fromtimestamp(closest_time, timezone.utc)
		alt = (alt[index]+alt[index+1])/2
		sp = (speed[index]+speed[index+1])/2

		alt_m = alt * 0.3048
		speed_mps = sp * 0.514444
		dist_m = dist_km * 1000
	else:
		continue
	
	closest_x_m = closest_x / 1000
	closest_y_m = closest_y / 1000
	station_x_m = station_x / 1000
	station_y_m = station_y / 1000

	# Convert UTM coordinates to latitude and longitude
	lon_p, lat_p = utm_proj(closest_x_m, closest_y_m, inverse=True)
	lon_s, lat_s = utm_proj(closest_x_m, closest_y_m, inverse=True)

	h = time_UTC.hour
	if h != 10:
		continue
	year = time_UTC.year
	month = time_UTC.month
	day = time_UTC.day
	time_stamp = time_UTC.timestamp()

	output = str(time_stamp) + '_' + str(lat_p) + '_' + str(lon_p) + '.dat' 

	comand = f'ncpag2s.py line --date {year}-{month}-{day} --hour {h} --startlat {lat_p} --startlon {lon_p} --endlat {lat_s} --endlon {lon_s} --points 21 --output {output}'

	os.system(comand)