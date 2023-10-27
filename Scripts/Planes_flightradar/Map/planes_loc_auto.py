import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.core import UTCDateTime
import datetime
import pytz
import obspy
import math 
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

def dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes):
	f = False
	for s in range(len(flight_latitudes)):
		for l in range(len(seismo_latitudes)):
			dist = distance(seismo_latitudes[l], seismo_longitudes[l], flight_latitudes[s], flight_longitudes[s])
			if dist <= 2:
				f = True
				break
			else:
				continue
	return f

write_out = open('stats_w.txt','w')

min_lon = -150.5
max_lon = -148.5
min_lat = 62.227
max_lat = 64.6


flight_files=[]
filenames = []

# Load the seismometer location data
seismo_data = pd.read_csv('nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_stations = seismo_data['Latitude']
sta = seismo_data['Station']

for month in range(2,4):
	if month == 2:
		month = '02'
		for day in range(11,29):
			day = str(day)
			# assign directory
			directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_positions'

			# iterate over files in directory
			for filename in os.listdir(directory):
				filenames.append(filename)
				f = os.path.join(directory, filename)
				
				# checking if it is a file
				if os.path.isfile(f):
					flight_files.append(f)
	elif month == 3:
		month = '03'
		for day in range(1, 27):
			if day < 10:
				day = '0' + str(day)
				# assign directory
				directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_positions'
			
				# iterate over files in directory
				for filename in os.listdir(directory):
					filenames.append(filename)
					f = os.path.join(directory, filename)
					
					# checking if it is a file
					if os.path.isfile(f):
						flight_files.append(f)
			else:
				day = str(day)
				# assign directory
				directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_positions'
				
				# iterate over files in directory
				for filename in os.listdir(directory):
					filenames.append(filename)
					f = os.path.join(directory, filename)
					
					# checking if it is a file
					if os.path.isfile(f):
						flight_files.append(f)
				

for i, flight_file in enumerate(flight_files):
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']
	fname = filenames[i]	
	flight_num = fname[:19]

	con = dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes)
	if con == True:	
		write_out.writelines("In file "+filenames[i]+":\n")

		dist_on = 2

		fig = plt.figure(figsize=(24,30))
		# Create a scatter plot for the seismometer locations
		for sd in range(len(seismo_data)):
			plt.scatter(seismo_longitudes[sd], seismo_latitudes[sd], c='red')		
		for fd in range(len(flight_data)):
			plt.scatter(flight_longitudes[fd], flight_latitudes[fd], c='b')
		# Create a scatter plot for the seismometer locations
		for sd in range(len(seismo_data)):		
			for fd in range(len(flight_data)):
				dist = distance(seismo_latitudes[sd], seismo_longitudes[sd], flight_latitudes[fd], flight_longitudes[fd])
				if dist <= 2:
				
					station = str(sta[sd])

					ht = datetime.datetime.utcfromtimestamp(time[fd])
					
					write_out.writelines("Station "+str(sta[sd])+" is "+str(dist)+" km away from the nearest time stamp at time "+str(ht)+'\n')
					
					#Label stations
					plt.text(seismo_longitudes[sd], seismo_latitudes[sd], sta[sd], fontsize=6)
					plt.scatter(seismo_longitudes[sd], seismo_latitudes[sd], c='pink')

					#Label time stamps with time
					plt.scatter(flight_longitudes[fd], flight_latitudes[fd], c='orange')
						
					if dist < dist_on:
						dist_on = dist
						index_p = fd
						index_s = sd
						time_use = ht
			
					else:
						continue
					
				else:
					continue
	
		plt.text(flight_longitudes[index_p], flight_latitudes[index_p], time_use, fontsize=8)
		t = 'At station '+str(sta[index_s])+': Speed '+str(speed[index_p])+'kts at '+str(alt[index_p])+'ft'
		plt.text(-150.45, 64.5, t, fontsize = 20)

		# Set labels and title
		plt.xlim(min_lon, max_lon)
		plt.ylim(min_lat, max_lat)
		plt.xlabel('Longitude')
		plt.ylabel('Latitude')
		plt.title(filenames[i])
		mon = str(flight_num[4:6])
		da = str(flight_num[6:8])
		BASE_DIR = "/scratch/irseppi/nodal_data/Plane_info/Plane_map/2019-"+mon+"-"+da
		make_base_dir(BASE_DIR)
						
		plt.savefig('/scratch/irseppi/nodal_data/Plane_info/Plane_map/2019-'+mon+'-'+da+'/map_'+flight_num+'png')

	print((i/len(flight_files))*100, '% Done')	
