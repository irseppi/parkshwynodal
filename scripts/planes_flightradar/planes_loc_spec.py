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

def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km

def dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes):
	f = False
	for s in range(len(flight_latitudes)):
		for l in range(len(seismo_latitudes)):
			dist = distance(seismo_latitudes[l], seismo_longitudes[l], flight_latitudes[s], flight_longitudes[s])
			if dist <= 5:
				f = True
				break
			else:
				continue
	return f

def calc_time(t, l, vo):
	#t is epoch time at time wave is generated by aircraft
	#l is the shortest distance to between the station and aircraft
	#vo is the aircraft velocity
	#c is the speed of aucostic wave propogation
	c = 0.343
	to=t+(math.sqrt(l**2+(vo*t)**2))/c
	return to
#Read in color text file to get different flights to be diffrent colors on map
with open('colors.txt','r') as c_in:
	for line in c_in:
		c=str(line[0:7])
		color.append(c)
write_out = open('stats_w.txt','w')

min_lon = -150.5
max_lon = -148.5
min_lat = 62.227
max_lat = 64.6


flight_files=[]
filenames = []

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

for month in (2,4):
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
	print(flight_file)
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']

	con = dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes)
	if con == True:	
		write_out.writelines("In file "+filenames[i]+":\n")

		spec_stations =[]
		spec_hours1 = []
		spec_hours2 =[]
		spec_month1 = []
		spec_month2 =[]
		spec_day1 =[]
		spec_day2 =[]
		spec_mins = []
		spec_sec = []
		dist_ps = []
		speed_p = []

		dist_on = 5

		fig = plt.figure(figsize=(24,30))
		# Create a scatter plot for the seismometer locations
		for sd in range(len(seismo_data)):
			plt.scatter(seismo_longitudes[sd], seismo_latitudes[sd], c='red')		
		for fd in range(len(flight_data)):
			plt.scatter(flight_longitudes[fd], flight_latitudes[fd], c='b')
		# Create a scatter plot for the seismometer locations
		for sd in range(len(seismo_data)):
			print(sta[sd])		
			for fd in range(len(flight_data)):
				dist = distance(seismo_latitudes[sd], seismo_longitudes[sd], flight_latitudes[fd], flight_longitudes[fd])
				if dist <= 5:
					print(dist)
					dist_ps.append(dist)
					speed_p.append(speed[fd])
					station = str(sta[sd])
					spec_stations.append(station)

					ht = datetime.datetime.utcfromtimestamp(time[fd])
					
					write_out.writelines("Station "+str(sta[sd])+" is "+str(dist)+" km away from the nearest time stamp at time "+str(ht)+'\n')
					
					h = ht.hour
					month = ht.month
					day = ht.day
					mins = ht.minute
					secs = ht.second
					spec_mins.append(mins)
					spec_sec.append(secs)
					
					month2 = str(month)
					if month == 3 and day < 10:
						day1 = '0'+str(day)
						print(month, day1)
					else:
						day1 = str(day)
						print(month, day1)
					if h < 23:
						h_u = str(h+1)			
						day2 = day1
						print(day2, h_u)
					else:
						h_u = '00'
						if month == '02' and day == '28':
							month2 = '03'
							day2 = '01'
							print(month2, day2, h_u) 
						else:
							day2 = str(day + 1)
							print(month2, day2, h_u) 
					
					spec_hours1.append(str(h))
					spec_hours2.append(h_u)
					spec_month1.append(str(month))
					spec_month2.append(month2)
					spec_day1.append(day1)
					spec_day2.append(day2)
					
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
						print(dist, fd, sd)
					else:
						print('where')
						continue
					
				else:
					continue
	
		plt.text(flight_longitudes[index_p], flight_latitudes[index_p], time_use, fontsize=8)
		t = 'At station '+str(sta[index_s])+': Speed '+str(speed[index_p])+'km/h at '+str(alt[index_p])+'ft'
		plt.text(-150.45, 64.5, t, fontsize = 20)

		# Set labels and title
		plt.xlim(min_lon, max_lon)
		plt.ylim(min_lat, max_lat)
		plt.xlabel('Longitude')
		plt.ylabel('Latitude')
		plt.title(filenames[i])
		plt.savefig('/scratch/irseppi/nodal_data/Plane_map_spec/map_'+filenames[i]+'.png')

		for spec in range(len(spec_stations)):
			#here you see spectrograms add flight data (ie. return type of plane, speed, altitude approximate arrival over station(plot as line on spectrogram)
			n = "/scratch/naalexeev/NODAL/2019-0"+spec_month1[spec]+"-"+spec_day1[spec]+"T"+spec_hours1[spec]+":00:00.000000Z.2019-0"+spec_month2[spec]+"-"+spec_day2[spec]+"T"+spec_hours2[spec]+":00:00.000000Z."+spec_stations[spec]+".mseed"
			print(n)
			if os.path.isfile(n):
				print('made it')
				tr = obspy.read(n)
				t = spec_mins[spec] * 60 + 60 + spec_sec[spec]
				l = dist_ps[spec]
				vo = speed_p[spec]
				xloc = calc_time(t, l, vo)

				tr[2].trim(tr[2].stats.starttime + spec_mins[spec] * 60 - 60 + spec_sec[spec], tr[2].stats.starttime + spec_mins[spec] * 60 + 60 + spec_sec[spec])
				fig, ax = plt.subplots()
				ax.set_xlabel('Time')
				ax.axvline(x=xloc, ls = '--')
				ax.set_ylabel('Frequency (Hz)')
				ax.xaxis_date('UTC')
				#ax.set_title('Time stamp: '+tr[2].stats.starttime + spec_mins[spec] * 60 - 60 + spec_sec[spec])
				# Save the plot
				tr[2].spectrogram(axes = ax,log=False, outfile = 'Plane_map_spec/spec_'+spec_stations[spec]+'_'+filenames[i]+'.png', fmt='.png', dbscale=True,cmap='hsv', show=False)
				#, title = tr[2].stats.starttime + spec_mins[spec] * 60 - 60 + spec_sec[spec])
				fig.saveflig('/scratch/irseppi/nodal_data/Plane_map_spec/spec_'+spec_stations[spec]+'_'+filenames[i]+'.png')

				tr[2].plot(outfile = '/scratch/irseppi/nodal_data/Plane_map_spec/trace_'+spec_stations[spec]+'_'+filenames[i]+'.png', show = False)
	print(i/len(flight_files), '% Done')	
write_out.close()	