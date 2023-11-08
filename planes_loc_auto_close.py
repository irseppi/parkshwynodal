import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.core import UTCDateTime
import datetime
from pathlib import Path

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()
		
def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km

sta_f = open('input/flight_name_sta.txt','r')


# Load the seismometer location data
seismo_data = pd.read_csv('input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
sta = seismo_data['Station']

#Loop through each station in text file that we already know comes within 2km of the nodes			
for line in sta_f.readlines():
	val = line.split(',')
	date = val[0]
	flight = val[1]
	station = val[2]

	flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + date + '_positions/' + date + '_' + flight + '.csv'
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']

	for line  in range(len(sta)):
		if int(sta[line]) == int(station):
			
			min_lon = seismo_longitudes[line]-0.05
			max_lon = seismo_longitudes[line]+0.05
			min_lat = seismo_latitudes[line]-0.05
			max_lat = seismo_latitudes[line]+0.05
			
			for l in range(len(flight_latitudes)):
				dist = distance(seismo_latitudes[line], seismo_longitudes[line], flight_latitudes[l], flight_longitudes[l])
				if dist <= 2:
					fig = plt.figure() #figsize=(24,30))

					# Create a scatter plot for the seismometer locations
					for sd in range(len(seismo_data)):
						plt.scatter(seismo_longitudes[sd], seismo_latitudes[sd], c='red')
					# Create a scatter plot for the timestamp locations		
					for fd in range(len(flight_data)):
						plt.scatter(flight_longitudes[fd], flight_latitudes[fd], c='b')

					#Label station
					plt.text(seismo_longitudes[line], seismo_latitudes[line], sta[line], fontsize=10)
					plt.scatter(seismo_longitudes[line], seismo_latitudes[line], c='pink')
					
					ht = datetime.datetime.utcfromtimestamp(time[l])

					#Label timestamp 
					plt.text(flight_longitudes[l], flight_latitudes[l], ht, fontsize=10)
					plt.scatter(flight_longitudes[l], flight_latitudes[l], c='orange')
					
					y =[flight_latitudes[l],  seismo_latitudes[line]]
					x = [flight_longitudes[l], seismo_longitudes[line]]
				
					plt.plot(x,y, '--', c='orange')
					yy = sum(y)/len(y)
					xx = sum(x)/len(x)
					plt.text(xx,yy, str(round(dist, 2))+'km', fontsize=10)

					# Set labels and title
					plt.xlim(min_lon, max_lon)
					plt.ylim(min_lat, max_lat)
					plt.xlabel('Longitude')
					plt.ylabel('Latitude')

					#Save
					plt.title('Date: ' + date + ' | Flight: ' + flight + ' | Station: ' + station + ' | Speed: '+str(speed[l])+'knts | Altitude: '+str(alt[l])+'ft')
					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/map_zoom/' + date + '/'+flight+'/'
					make_base_dir(BASE_DIR)
					plt.savefig('/scratch/irseppi/nodal_data/plane_info/map_zoom/'+ date + '/'+flight+'/'+flight+'_'+station+'_' + str(time[l]) + '.png')
					plt.close()
					
				else:
					continue
			
		else:
			continue

	

