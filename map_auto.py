import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.projections import geo
import os
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.core import UTCDateTime
import datetime
from pathlib import Path
import pytz
import obspy

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

sta_f = open('input/all_station_crossing_db.txt','r')

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_stations = seismo_data['Latitude']
sta = seismo_data['Station']

flight_num = 0
#Loop through each station in text file that we already know comes within 2km of the nodes			
for line in sta_f.readlines():
	val = line.split(',')
	date = val[0]
	flight = val[1]
	station = val[5]
	if flight != flight_num:
		flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + date + '_positions/' + date + '_' + flight + '.csv'
		flight_data = pd.read_csv(flight_file, sep=",")
		flight_latitudes = flight_data['latitude']
		flight_longitudes = flight_data['longitude']
		
		f = 1.0/np.cos(62*np.pi/180)

		fig = plt.figure()#figsize=(24,30))
		plt.gca().set_aspect(f)
		
		plt.scatter(seismo_longitudes, seismo_latitudes, c='red', s = 30, label='seismometers')
		plt.plot(flight_longitudes, flight_latitudes, '-o', c='c', lw=3, ms = 6, label='flight path')
		
		# Set labels and title
		plt.xlim(min_lon, max_lon)
		plt.ylim(min_lat, max_lat)
		plt.yticks(fontsize = 14)
		plt.xticks(fontsize = 14)
		
		#Save
		BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/pmap2/' + date + '/'
		make_base_dir(BASE_DIR)
		plt.savefig('/scratch/irseppi/nodal_data/plane_info/pmap2/'+ date + '/map_'+flight+'.png')
		flight_num = flight
		plt.close()	
	else:
		continue
	

