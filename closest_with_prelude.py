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

def closest_encounter(flight_latitudes, flight_longitudes, index, timestamp, seismo_latitude, seismo_longitude):

	clat = flight_latitudes[index]
	clon = flight_longitudes[index]
		
	closest_lat = clat
	closest_lon = clon
	dist_lim = 2.01
	c2lat = flight_latitudes[index+1]
	c2lon = flight_longitudes[index+1]
	index2 = index + 1
	for tr in range(0,2):
		if index == 0 or index == -1:
			continue
		if  flight_longitudes[index+1] < flight_longitudes[index-1]:
			if tr == 0:
				sclat = flight_latitudes[index-1]
				sclon = flight_longitudes[index-1]  

				x = [clon, sclon]
				y = [clat, sclat]
				m = (y[0]-y[1])/(x[0]-x[1])
				b = y[0] - m*x[0]
				for point in np.arange(clon, sclon, 0.000001):
					lon = point
					lat = m*lon + b
					dist_km = distance(lat, lon, seismo_latitude, seismo_longitude)

					if dist_km < dist_lim:
						dist_lim = dist_km
						closest_lat = lat
						closest_lon = lon
						c2lat = flight_latitudes[index-1]
						c2lon = flight_longitudes[index-1]
						index2 = index - 1
			elif tr == 1:
				sclat = flight_latitudes[index+1]
				sclon = flight_longitudes[index+1]

				x = [clon, sclon]
				y = [clat, sclat]
				m = (y[1]-y[0])/(x[1]-x[0])
				b = y[0] - m*x[0]


				for point in np.arange(sclon, clon, 0.000001):
					lon = point
					lat = m*lon + b
					dist_km = distance(lat, lon, seismo_latitude, seismo_longitude)

					if dist_km < dist_lim:
						dist_lim = dist_km
						closest_lat = lat
						closest_lon = lon
						c2lat = flight_latitudes[index+1]
						c2lon = flight_longitudes[index+1]
						index2 = index + 1
		else:
			if tr == 0:
				sclat = flight_latitudes[index-1]
				sclon = flight_longitudes[index-1]  

				x = [clon, sclon]
				y = [clat, sclat]
				m = (y[0]-y[1])/(x[0]-x[1])
				b = y[0] - m*x[0]
				for point in np.arange(sclon, clon, 0.000001):
					lon = point
					lat = m*lon + b
					dist_km = distance(lat, lon, seismo_latitude, seismo_longitude)

					if dist_km < dist_lim:
						dist_lim = dist_km
						closest_lat = lat
						closest_lon = lon
						c2lat = flight_latitudes[index-1]
						c2lon = flight_longitudes[index-1]
						index2 = index - 1
			elif tr == 1:
				sclat = flight_latitudes[index+1]
				sclon = flight_longitudes[index+1]

				x = [clon, sclon]
				y = [clat, sclat]
				m = (y[1]-y[0])/(x[1]-x[0])
				b = y[0] - m*x[0]


				for point in np.arange(clon, sclon, 0.000001):
					lon = point
					lat = m*lon + b
					dist_km = distance(lat, lon, seismo_latitude, seismo_longitude)

					if dist_km < dist_lim:
						dist_lim = dist_km
						closest_lat = lat
						closest_lon = lon
						c2lat = flight_latitudes[index+1]
						c2lon = flight_longitudes[index+1]
						index2 = index + 1
	if dist_lim < 2.01:
		for location in np.arange((closest_lon-0.000001),(closest_lon+0.000001),0.0000000001):
			lon = location
			lat = m*lon + b
			dist = distance(lat, lon, seismo_latitude, seismo_longitude)
			if dist < dist_lim:
				dist_lim = dist
				closest_lon = lon
				closest_lat = lat
		if dist_lim < 1:
			for location in np.arange((closest_lon-0.0000000001),(closest_lon+0.0000000001),0.0000000000001):
				lon = location
				lat = m*lon + b
				dist = distance(lat, lon, seismo_latitude, seismo_longitude)
				if dist < dist_lim:
					dist_lim = dist
					closest_lon = lon
					closest_lat = lat
		dist_old_new = distance(closest_lat, closest_lon, clat, clon)
		dist_old_old = distance(c2lat, c2lon, clat, clon)
		ratio = dist_old_new/dist_old_old
		timestamp = timestamp[index] + (timestamp[index] - timestamp[index2])*ratio

	return  closest_lat, closest_lon, dist_lim, timestamp

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3
cenlat = (min_lat + max_lat)/2
cenlon = (min_lon + max_lon)/2

output = open('all_station_crossing_db_updated1.txt','w')
input = open('input/all_station_crossing_db.txt','r')

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_sta = seismo_data['Station']
		
for line in input.readlines():
	lin = line.split(',')
	flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/'+lin[0] + '_positions/'+lin[0]+'_'+lin[1]+'.csv', sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	alt = flight_data['altitude']
	heading = flight_data['heading']
	speed = flight_data['speed']
	flight_num = lin[1]
	date = lin[0]

	seismometer = lin[6]
	for s in range(len(seismo_sta)):
		if seismometer == seismo_sta[s]:
			dist_lim = 2.01
			con = False
			for l in range(len(flight_data)):
				if l == 0 or l == -1:
					continue
				index = l
				timestamp = time[l]
				closest_info = closest_encounter(flight_latitudes, flight_longitudes, index, time, seismo_latitudes[s], seismo_longitudes[s])
				dist = closest_info[2]
				if dist <= dist_lim:
					timeF = closest_info[3]
					altF = alt[l]
					speedF = speed[l]
					seismo_staF = seismo_sta[s]
					headF = heading[l]
					con = True
					closest_lat = closest_info[0]
					closest_lon = closest_info[1]
					dist_lim = dist
					c2lat = flight_latitudes[l]
					c2lon = flight_longitudes[l]
				else:
					continue
	if con == True:
		text = open('/scratch/irseppi/nodal_data/flightradar24/'+date+'_flights.csv')
		for ln in text.readlines():
			val = ln.split(',')
			if val[0] == flight_num:
				print('Found')
				output.write(str(date)+','+ str(flight_num)+','+val[1]+','+str(timeF)+','+str(altF)+','+str(speedF)+','+str(seismo_staF)+','+closest_lat+','+closest_lon+','+val[3]+',\n')
		# Create a figure with two subplots side by side
		fig, axs = plt.subplots(1, 2) 
		fig.subplots_adjust(wspace=0.5)  # Adjust the spacing between subplots


		y =[closest_lat,  seismo_latitudes[s]]
		x = [closest_lon, seismo_longitudes[s]]
		yy = sum(y)/len(y)
		xx = sum(x)/len(x)
		if dist_lim < 0.1:
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


		heading = np.deg2rad(headF)

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

		m = (closest_lat - c2lat)/(closest_lon - c2lon)
		b = closest_lat - m*closest_lon
		direction = np.arctan2(c2lat - c2lat, c2lon - c2lon)

		axs[1].quiver(closest_lon, closest_lat, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=7)

		axs[1].quiver(closest_lon, closest_lat, np.cos(heading), np.sin(heading), angles='xy', scale = 7, color='#999999')
		axs[1].scatter(closest_lon, closest_lat, c='#377eb8', s=50, zorder=3)
		axs[1].scatter(seismo_longitudes[s], seismo_latitudes[s], c='#e41a1c', s=50, zorder=3)

		axs[1].text(xx,yy, str(round(dist_lim, 3))+' km', fontsize=15, fontweight='bold')

		# Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
		con = mpatch.ConnectionPatch(xyA=(minl, minla), xyB=(maxl, minla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
		fig.add_artist(con)
		con = mpatch.ConnectionPatch(xyA=(minl, maxla), xyB=(maxl, maxla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
		fig.add_artist(con)
		plt.show()
	else:
		continue
	print(line)	
input.close()
output.close()	
