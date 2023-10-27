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

def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km

def convert_UTC_to_epoch(timestamp):
	tz_UTC = pytz.timezone('UTC')
	time_format = "%Y-%m-%d %H:%M:%S"
	naive_timestamp = datetime.datetime.strptime(timestamp, time_format)
	aware_timestamp = tz_UTC.localize(naive_timestamp)
	epoch = aware_timestamp.strftime("%s")
	return (int) (epoch)
flight_files=[]
filenames = []
#Pick date to view
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
					flight_files=[]
					filenames = []
					# iterate over files in directory
					for filename in os.listdir(directory):
						filenames.append(filename)
						f = os.path.join(directory, filename)
						
						# checking if it is a file
						if os.path.isfile(f):
							flight_files.append(f)

sstations = 1001
estations = 1306

min_lon = -162.5
max_lon = -142.0
min_lat = 50
max_lat = 68


s_stamp = "00:00:00"
e_stamp = "23:59:59"
s_epoch = convert_UTC_to_epoch("2019-"+month+"-"+day+" "+s_stamp)
e_epoch = convert_UTC_to_epoch("2019-"+month+"-"+day+" "+e_stamp)

color=[]
#Read in color text file to get different flights to be diffrent colors on map
with open('colors.txt','r') as c_in:
		for line in c_in:
			c=str(line[0:7])
			color.append(c)


# Load the seismometer location data
seismo_data = pd.read_csv('nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_stations = seismo_data['Latitude']
sta = seismo_data['Station']

for i, flight_file in enumerate(flight_files):
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']
	print("Wait for map")
	for s in range(len(flight_data)-1):
		if s_epoch <= int(time[s]) <= e_epoch and min_lat <= flight_latitudes[s] <= max_lat and min_lon <= flight_longitudes[s] <= max_lon:
			for l in range(len(seismo_data)):
				if sstations <= sta[l] <= estations:
				
					station = str(sta[l])
					dist = distance(seismo_latitudes[l], seismo_longitudes[l], flight_latitudes[s], flight_longitudes[s])
					print('...')
					if dist <= 5:
					
					
						t = UTCDateTime(float(time[s]))
						ht = datetime.datetime.fromtimestamp(time[s])

						if month == '03' and int(day) > 10:
							h = ht.hour+8
						else: 
							h = ht.hour+9
						h_u = str(h+1)
						
						print("In file "+filenames[i]+":")
						fig = plt.figure()
						# Create a scatter plot for the seismometer locations
						for sd in range(len(seismo_data)):
							plt.scatter(seismo_longitudes[sd], seismo_latitudes[sd], c='red')
						#plt.scatter(flight_longitudes, flight_latitudes, 'o-',c='blue')
						for fd in range(len(flight_data)-1):
							plt.scatter(flight_longitudes[fd], flight_latitudes[fd], c=color[i % len(color)])
				
						for stations_up in range(len(seismo_data)):
							for flights_up in range(len(flight_data)-1):
								dist = distance(seismo_latitudes[stations_up], seismo_longitudes[stations_up], flight_latitudes[flights_up], flight_longitudes[flights_up])
								if dist <= 5:
									htn = str(UTCDateTime(float(time[flights_up])))
									htnn = datetime.datetime.fromtimestamp(time[flights_up])
									
									print("Station", sta[stations_up], "is", dist,"km away from the nearest time stamp at time "+htn)
									
									
									#Label stations
									plt.text(seismo_longitudes[stations_up], seismo_latitudes[stations_up], sta[stations_up], fontsize=5)
									
									plt.scatter(seismo_longitudes[stations_up], seismo_latitudes[stations_up], c='pink')
									#Label time stamps with epoch time
									plt.text(flight_longitudes[flights_up], flight_latitudes[flights_up], htn, fontsize=5)
									plt.scatter(flight_longitudes[flights_up], flight_latitudes[flights_up], c='orange')
									xx = np.vstack([seismo_longitudes[stations_up], seismo_latitudes[stations_up]])
									yy = np.vstack([flight_longitudes[flights_up],flight_latitudes[flights_up]])
									
								else:
									continue
										
						plt.plot(xx,yy, '-.', c='orange')
						# Set labels and title
						plt.xlim(-153, -142)
						plt.ylim(60, 65)
						plt.xlabel('Longitude')
						plt.ylabel('Latitude')
						plt.title(filenames[i])	
						plt.show(block=False)
						
						#add arival over location here
						#here you see spectrograms add flight data (ie. return type of plane, speed, altitude approximate arrival over station(plot as line on spectrogram
						
						sta_spec = 
						stations_up = int(int(sta_spec)-1000)
						index_on = 0
						dist_on = 100
						for flights_up in range(len(flight_data)-1):
							dist = distance(seismo_latitudes[stations_up], seismo_longitudes[stations_up], flight_latitudes[flights_up], flight_longitudes[flights_up])
							if dist < dist_on:
								dist_on = dist
								index_on = flights_up
							else:
								continue
						htn = str(UTCDateTime(float(time[index_on])))
						htnn = datetime.datetime.fromtimestamp(time[index_on])
						if month == '03' and int(day) > 10:
							h = htnn.hour+8
						else: 
							h = htnn.hour+9
						h_u = str(h+1)
						n = "/scratch/naalexeev/NODAL/2019-"+month+"-"+day+"T"+str(h)+":00:00.000000Z.2019-"+month+"-"+day+"T"+h_u+":00:00.000000Z."+str(sta_spec)+".mseed"
						if os.path.isfile(n):
							
							tr = obspy.read(n)
							tr.trim(tr[2].stats.starttime + int(str(htnn.minute))*60 - 60+int(str(htnn.second)), tr[2].stats.starttime + int(str(htnn.minute))*60 + 60+int(str(htnn.second)))
							fig1, ax1 = plt.subplots()
							tr[2].spectrogram(axes = ax1,log=False,dbscale=True,cmap='hsv')
							fig1.show()

							# Show the plot
							tr[2].plot()
							print("The plane is traveling at an altitude of", alt[index_on], "at ", speed[index_on] ,"km per hour.")
						else:
							print('No data for this station')
					
						break
					else:
						continue
					break
					
				else:
					continue
		else:
			continue
			
