import pandas as pd
import os
from obspy.geodetics import gps2dist_azimuth
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from matplotlib.patches import Rectangle

def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km
def func(y):
    return -1.0/np.cos((y)*np.pi/180)
def rfunc(y):
    return -180/np.pi*np.arccos(1/y)
flight_files=[]
filenames = []

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3
cenlat = (min_lat + max_lat)/2
cenlon = (min_lon + max_lon)/2

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

#output = open('all_station_crossing_db_updated.txt','w')

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_sta = seismo_data['Station']

# Iterate over flight files
for i, flight_file in enumerate(flight_files):
	print((i/len(flight_files))*100, '%')	
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
		dist_lim = 2.1
		con = False
		# Iterate over flight data
		for l in range(len(flight_data)):
			if l == 0:
				continue
			elif l >= len(flight_data)-1:
				continue

			clat = flight_latitudes[l]
			clon = flight_longitudes[l]
			
			closest_lat = clat
			closest_lon = clon
			c2lat = flight_latitudes[l+1]
			c2lon = flight_longitudes[l+1]
			l2 = l + 1
			# Iterate over two points for linear interpolation
			for tr in range(0,2):
				if  flight_longitudes[l+1] < flight_longitudes[l-1]:
					if tr == 0:
						sclat = flight_latitudes[l-1]
						sclon = flight_longitudes[l-1]  

						x = [clon, sclon]
						y = [clat, sclat]
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
					elif tr == 1:
						sclat = flight_latitudes[l+1]
						sclon = flight_longitudes[l+1]

						x = [clon, sclon]
						y = [clat, sclat]
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
				else:
					if tr == 0:
						sclat = flight_latitudes[l-1]
						sclon = flight_longitudes[l-1]  

						x = [clon, sclon]
						y = [clat, sclat]
						m = (y[0]-y[1])/(x[0]-x[1])
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
								c2lat = flight_latitudes[l-1]
								c2lon = flight_longitudes[l-1]
								l2 = l - 1
							else:
								continue
					elif tr == 1:
						sclat = flight_latitudes[l+1]
						sclon = flight_longitudes[l+1]

						x = [clon, sclon]
						y = [clat, sclat]
						m = (y[1]-y[0])/(x[1]-x[0])
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
								c2lat = flight_latitudes[l+1]
								c2lon = flight_longitudes[l+1]
								l2 = l - 1
							else:
								continue
		for location in np.arange((closest_lon-0.000001),(closest_lon+0.000001),0.0000000001):
			lon = location
			lat = m*lon + b
			dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])
			if dist_km < dist_lim:
				dist_lim = dist_km
				closest_lon = lon
				closest_lat = lat
		if dist_lim < 1:
			for location in np.arange((closest_lon-0.0000000001),(closest_lon+0.0000000001),0.0000000000001):
				lon = location
				lat = m*lon + b
				dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])
				if dist_km < dist_lim:
					dist_lim = dist_km
					closest_lon = lon
					closest_lat = lat
		dist_old_new = distance(closest_lat, closest_lon, clat, clon)
		dist_old_old = distance(c2lat, c2lon, clat, clon)
		ratio = dist_old_new/dist_old_old
		timestamp = timestamp[l] + (timestamp[l] - timestamp[l2])*ratio

		timeF = timestamp
		altF = (alt[l2]+alt[l])/2
		speedF = (speed[l2]+speed[l])/2
		seismo_staF = seismo_sta[s]
		con = True
		# Create a figure with two subplots side by side
		fig, axs = plt.subplots(1, 2) 
		fig.subplots_adjust(wspace=0.5)  # Adjust the spacing between subplots


		y =[closest_lat,  seismo_latitudes[s]]
		x = [closest_lon, seismo_longitudes[s]]
		yy = sum(y)/len(y)
		xx = sum(x)/len(x)
		if dist_km < 0.1:
			minl = np.round((xx - 0.001),4)
			maxl = np.round((xx + 0.001),4)
			minla = np.round((yy - 0.001),4)
			maxla =  np.round((yy + 0.001),4)
			axs[1].set_xticks(np.arange(minl, maxl, 0.0005))
			axs[1].set_yticks(np.arange(minla, maxla, 0.0005))
		else:
			minl = np.round((xx - 0.025), 2)
			maxl = np.round((xx + 0.025), 2)
			minla = np.round((yy - 0.025), 2)
			maxla =  np.round((yy + 0.025), 2)
			axs[1].set_xticks(np.arange(minl, maxl, 0.01))
			axs[1].set_yticks(np.arange(minla, maxla, 0.01))
		axs[0].set_xticks(np.arange(min_lon, max_lon, 1))
		axs[0].set_yticks(np.arange(min_lat, max_lat, 1))
		axs[0].set_xscale('function', functions=(func,rfunc))
		axs[1].set_xscale('function', functions=(func,rfunc))
		axs[0].grid(True,linestyle='dotted',color='gray')
		axs[1].grid(True,linestyle='dotted',color='gray')
		axs[0].scatter(seismo_longitudes, seismo_latitudes, c='#e41a1c', s = 3, label='seismometers')
		axs[0].plot(flight_longitudes, flight_latitudes, '-', c='#377eb8', lw=1, ms = 1, label='flight path')
		for i in range(1, len(flight_latitudes)-1, int(len(flight_latitudes)/5)):
			direction = np.arctan2(flight_latitudes[i+1] - flight_latitudes[i], flight_longitudes[i+1] - flight_longitudes[i])
			m = (flight_latitudes[i+1] - flight_latitudes[i])/(flight_longitudes[i+1] - flight_longitudes[i])
			b = flight_latitudes[i] - m*flight_longitudes[i]
			axs[0].quiver((flight_latitudes[i]-b)/m, flight_latitudes[i], np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', headwidth = 5)
		# Set labels and title
		axs[0].set_xlim(min_lon, max_lon)
		axs[0].set_ylim(min_lat, max_lat)
		axs[0].tick_params(axis='both', which='major', labelsize=12)
		

		heading = np.deg2rad(head[l])
		
		rect = Rectangle((minl, minla), 0.1, 0.06, ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
		axs[0].add_patch(rect)
		axs[1].plot(x,y, '--', c='#ff7f00')

		# Draw the zoomed in map on the second subplot
		axs[1].plot(flight_longitudes, flight_latitudes, c='#377eb8',linestyle ='dotted')
		axs[1].set_xlim(minl, maxl)
		axs[1].set_ylim(minla, maxla)
		axs[1].tick_params(axis='both', which='major', labelsize=9)
		axs[1].ticklabel_format(useOffset=False, style='plain')
		if len(axs[1].get_xticklabels()) > 4:
			axs[1].set_xticklabels([round(x, 4) for x in axs[1].get_xticks()], rotation=20, fontsize=9)
			c2lat, c2lon, clat, clon
		m = (c2lat - clat)/(c2lon - clon)
		b = clat - m*clon
		direction = np.arctan2(c2lat - clat, c2lon - clon)
			
		axs[1].quiver(closest_lon, closest_lat, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=7)

		axs[1].quiver(closest_lon, closest_lat, np.cos(heading), np.sin(heading), angles='xy', scale = 7, color='#999999')
		axs[1].scatter(closest_lon, closest_lat, c='#377eb8', s=50, zorder=3)
		axs[1].scatter(seismo_longitudes[s], seismo_latitudes[s], c='#e41a1c', s=50, zorder=3)

		axs[1].text(xx,yy, str(round(dist_km, 3))+' km', fontsize=15, fontweight='bold')

		# Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
		con = mpatch.ConnectionPatch(xyA=(minl, minla), xyB=(maxl, minla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
		fig.add_artist(con)
		con = mpatch.ConnectionPatch(xyA=(minl, maxla), xyB=(maxl, maxla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
		fig.add_artist(con)
		plt.show()
		
				
		#if con == True:
		#	text = open('/scratch/irseppi/nodal_data/flightradar24/'+date+'_flights.csv')
		#	for line in text.readlines():
		#		val = line.split(',')
		#		if val[0] == flight_num:
		#			print('Found')
		#			output.write(str(date)+','+ str(flight_num)+','+val[1]+','+str(timeF)+','+str(altF)+','+str(speedF)+','+str(seismo_staF)+','+val[3]+',\n')
		#else:
		#	continue


#output.close()	
