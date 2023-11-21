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
import matplotlib.dates as mdates
from pathlib import Path

def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km

flight_files=[]
filenames = []

for day in range(11,29):
	day = str(day)
	# assign directory
	directory = '/scratch/irseppi/nodal_data/flightradar24/201902'+day+'_positions'

	# iterate over files in directory
	for filename in os.listdir(directory):
		filenames.append(filename)
		f = os.path.join(directory, filename)

		# checking if it is a file
		if os.path.isfile(f):
			flight_files.append(f)

for day in range(1, 10):
	day = '0' + str(day)
	# assign directory
	directory = '/scratch/irseppi/nodal_data/flightradar24/201903'+day+'_positions'

	# iterate over files in directory
	for filename in os.listdir(directory):
		filenames.append(filename)
		f = os.path.join(directory, filename)

		# checking if it is a file
		if os.path.isfile(f):
			flight_files.append(f)
for day in range(10, 27):
	day = str(day)
	# assign directory
	directory = '/scratch/irseppi/nodal_data/flightradar24/201903'+day+'_positions'

	# iterate over files in directory
	for filename in os.listdir(directory):
		filenames.append(filename)
		f = os.path.join(directory, filename)

		# checking if it is a file
		if os.path.isfile(f):
			flight_files.append(f)

output = open('all_station_crossing_db.txt','w')

# Load the seismometer location data
seismo_data = pd.read_csv('all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_sta = seismo_data['Station']
		
for i, flight_file in enumerate(flight_files):
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	alt = flight_data['altitude']
	speed = flight_data['speed']
	fname = filenames[i]	
	flight_num = fname[9:18]
	date = fname[0:8]
	
	for l in range(len(flight_data)):
		for s in range(len(seismo_data)):
			dist = distance(flight_latitudes[l], flight_longitudes[l], seismo_latitudes[s], seismo_longitudes[s])
			if dist <= 2:
				text = open('/scratch/irseppi/nodal_data/flightradar24/'+date+'_flights.csv')
				for line in text.readlines():
					val = line.split(',')
					if val[0] == flight_num:		
						output.write(str(date)+','+ str(flight_num)+','+str(time[l])+','+str(alt[l])+','+str(speed[l])+','+str(seismo_sta[s])+','+val[3]+'\n')
	print((i/len(flight_files))*100, '% Done')	
output.close()	
