import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.projections import geo
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

 
sta_f = open('all_station_crossing_db.txt','r')

min_lon = -150.5
max_lon = -147.5
min_lat = 62
max_lat = 65.5

# Load the seismometer location data
seismo_data = pd.read_csv('perm_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
sta = seismo_data['Station']

#Loop through each station in text file that we already know comes within 2km of the nodes			
for line in sta_f.readlines():
	val = line.split(',')
	date = val[0]
	flight = val[1]
	station = val[5]
	
	if val[5].isdigit() == False:
		print(station)
		flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + date + '_positions/' + date + '_' + flight + '.csv'
		flight_data = pd.read_csv(flight_file, sep=",")
		flight_latitudes = flight_data['latitude']
		flight_longitudes = flight_data['longitude']
		time = flight_data['snapshot_id']
		speed = flight_data['speed']
		alt = flight_data['altitude']

		for line  in range(len(sta)):
			if sta[line] == station:
				for l in range(len(flight_latitudes)):
					dist = distance(seismo_latitudes[line], seismo_longitudes[line], flight_latitudes[l], flight_longitudes[l])
					if dist <= 2:
						y =[flight_latitudes[l],  seismo_latitudes[line]]
						x = [flight_longitudes[l], seismo_longitudes[line]]
						yy = sum(y)/len(y)
						xx = sum(x)/len(x)

						#f = (np.arccos((62/111.32)*np.pi/180)/110.32)
						f = 1.0/np.cos(62*np.pi/180)
	 
						fig = plt.figure(figsize=(24,30))
						plt.gca().set_aspect(f)
						
						plt.scatter(seismo_longitudes, seismo_latitudes, c='red')
						
						plt.scatter(flight_longitudes, flight_latitudes, c='b')

						#Label station
						plt.text(seismo_longitudes[line], seismo_latitudes[line], sta[line], fontsize=9, fontweight='bold')
						plt.scatter(seismo_longitudes[line], seismo_latitudes[line], c='pink')
						
						ht = datetime.datetime.utcfromtimestamp(time[l])

						#Label timestamp 
						plt.text(flight_longitudes[l], flight_latitudes[l], ht, fontsize=9, fontweight='bold')
						plt.scatter(flight_longitudes[l], flight_latitudes[l], c='orange')
						
						
						#plt.plot(x,y, '--', c='orange')
						
						#plt.text(xx,yy, str(round(dist, 2))+'km', fontsize=8, fontweight='bold')
						
						# Set labels and title
						plt.xlim(min_lon, max_lon)
						plt.ylim(min_lat, max_lat)
						#ax.tick_params(axis='both', which='major', labelsize=9)
						#Save
						#plt.title('Date: ' + date + ' | Flight: ' + flight + ' | Station: ' + station + '\n | Speed: '+str(round(speed[l]*0.514444,2))+'m/s | Altitude: '+str(round(alt[l]*0.3048,2))+'m')
						BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/pmap/' + date + '/'+flight+'/'+station+'/'
						make_base_dir(BASE_DIR)
						plt.savefig('/scratch/irseppi/nodal_data/plane_info/pmap/'+ date + '/'+flight+'/'+station+'/zmap_'+flight+'_' + str(time[l]) + '.png')
						plt.close()
						
					else:
						continue
				
			else:
				continue

