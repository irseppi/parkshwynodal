import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from matplotlib.patches import Rectangle
from prelude import distance, load_flights

# Load flight files and filenames
flight_files,filenames = load_flights(2, 4, 11, 27)

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_sta = seismo_data['Station']

#Assign the boundries for the full nodal array map
min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3
cenlat = (min_lat + max_lat)/2
cenlon = (min_lon + max_lon)/2

#Define functions for lat and lon scaling at higher latitudes
def f(y):
    return -1.0/np.cos((y)*np.pi/180)
def rf(y):
    return (-180*np.arccos(-1/y))/np.pi

# Iterate over flight files
for i, flight_file in enumerate(flight_files):
	print((i/len(flight_files))*100, '%')	# Print progress

	# Load flight data and extract relevant columns
	flight_data = pd.read_csv(flight_file, sep=",") 
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	timestamps = flight_data['snapshot_id']
	alt = flight_data['altitude']
	speed = flight_data['speed']
	head = flight_data['heading']
	fname = filenames[i]	
	flight_num = fname[9:18]
	date = fname[0:8]

	# Iterate over seismometer data
	for s in range(len(seismo_data)):
		dist_lim = 2.01
		dist_final = float('inf')
		con = False
		# Iterate over flight data
		for index in range(len(flight_data)):
			if index == 0:
				continue
			elif index >= len(flight_data)-1:
				continue
			else:
				# Check if the distance between the initial flight timestamp and seismometer is less than 2km
				dist = distance(flight_latitudes[index], flight_longitudes[index], seismo_latitudes[s], seismo_longitudes[s])
				if dist < 2.01:
					#
					clat = flight_latitudes[index]
					clon = flight_longitudes[index]
						
					closest_lat = clat
					closest_lon = clon

					for tr in range(0,2):
						if  flight_latitudes[index+1] < flight_latitudes[index-1]:
							if tr == 0:
								sclat = flight_latitudes[index-1]
								sclon = flight_longitudes[index-1]  

								x = [clon, sclon]
								y = [clat, sclat]
								m = (y[1]-y[0])/(x[1]-x[0])
								b = y[0] - m*x[0]
								
								for point in np.arange(clat, sclat, 0.000001):
									lat = point
									lon = (lat - b)/m
									dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])

									if dist_km < dist_lim:
										dist_lim = dist_km
										closest_lat = lat
										closest_lon = lon
										c2lat = flight_latitudes[index-1]
										c2lon = flight_longitudes[index-1]
										index2 = index - 1
									else:
										continue
							elif tr == 1:
								sclat = flight_latitudes[index+1]
								sclon = flight_longitudes[index+1]

								x = [clon, sclon]
								y = [clat, sclat]
								m = (y[0]-y[1])/(x[0]-x[1])
								b = y[0] - m*x[0]
								
								for point in np.arange(sclat, clat, 0.000001):
									lat = point
									lon = (lat - b)/m
									dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])

									if dist_km < dist_lim:
										dist_lim = dist_km
										closest_lat = lat
										closest_lon = lon
										c2lat = flight_latitudes[index+1]
										c2lon = flight_longitudes[index+1]
										index2 = index + 1
									else:
										continue
						elif flight_latitudes[index+1] > flight_latitudes[index-1]:
							if tr == 0:
								sclat = flight_latitudes[index-1]
								sclon = flight_longitudes[index-1]  

								x = [clon, sclon]
								y = [clat, sclat]
								m = (y[0]-y[1])/(x[0]-x[1])
								b = y[0] - m*x[0]

								for point in np.arange(sclat, clat, 0.000001):
									lat = point
									lon = (lat - b)/m
									dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])

									if dist_km < dist_lim:
										dist_lim = dist_km
										closest_lat = lat
										closest_lon = lon
										c2lat = flight_latitudes[index-1]
										c2lon = flight_longitudes[index-1]
										index2 = index - 1
									else:
										continue
							elif tr == 1:
								sclat = flight_latitudes[index+1]
								sclon = flight_longitudes[index+1]

								x = [clon, sclon]
								y = [clat, sclat]
								m = (y[1]-y[0])/(x[1]-x[0])
								b = y[0] - m*x[0]

								for point in np.arange(clat, sclat, 0.000001):
									lat = point
									lon = (lat - b)/m
									dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])

									if dist_km < dist_lim:
										dist_lim = dist_km
										closest_lat = lat
										closest_lon = lon
										c2lat = flight_latitudes[index+1]
										c2lon = flight_longitudes[index+1]
										index2 = index + 1
									else:
										continue
						else:
							continue
					
					if closest_lat != clat:
						for location in np.arange((closest_lon-0.000001),(closest_lon+0.000001),0.0000000001):
							lon = location
							lat = m*lon + b
							dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])
							if dist_km < dist_lim:
								dist_lim = dist_km
								closest_lon = lon
								closest_lat = lat
							else:
								continue
						if dist_lim < 1:
							for location in np.arange((closest_lon-0.0000000001),(closest_lon+0.0000000001),0.0000000000001):
								lon = location
								lat = m*lon + b
								dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])
								if dist_km < dist_lim:
									dist_lim = dist_km
									closest_lon = lon
									closest_lat = lat
								else:
									continue
						dist_old_new = distance(closest_lat, closest_lon, clat, clon)
						dist_old_old = distance(c2lat, c2lon, clat, clon)
						ratio = dist_old_new/dist_old_old
						timestamp = timestamps[index] + (timestamps[index] - timestamps[index2])*ratio
						index_final = index
					else:
						continue
				else:
					continue
			if dist_final < dist_lim:
				con = True
				dist_final = dist_lim
				print(dist_final)
			else:
				continue
		if con == True:             
			clat = closest_lat
			clon = closest_lon
			fig, axs = plt.subplots(1, 2) 
			fig.subplots_adjust(wspace=0.5)  

			y =[clat,  seismo_latitudes[s]]
			x = [clon, seismo_longitudes[s]]
			yy = sum(y)/len(y)
			xx = sum(x)/len(x)
			if dist_final < 0.1:
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
			axs[0].set_xscale('function', functions=(f,rf))
			axs[1].set_xscale('function', functions=(f,rf))
			axs[0].grid(True,linestyle='dotted',color='gray')
			axs[1].grid(True,linestyle='dotted',color='gray')
			axs[0].scatter(seismo_longitudes, seismo_latitudes, c='#e41a1c', s = 3, label='seismometers')
			axs[1].scatter(seismo_longitudes, seismo_latitudes, c='#e41a1c', s = 3, label='seismometers')
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
			
			heading = np.deg2rad(head[index_final])
			
			rect = Rectangle((minl, minla), (maxl-minl), (maxla-minla), ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
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
			m = (flight_latitudes[index_final+1] - flight_latitudes[index_final])/(flight_longitudes[index_final+1] - flight_longitudes[index_final])
			b = flight_latitudes[index_final] - m*flight_longitudes[index_final]

			direction = np.arctan2(flight_latitudes[index_final+1] - flight_latitudes[index_final], flight_longitudes[index_final+1] - flight_longitudes[index_final])
				
			axs[1].quiver(clon, clat, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=7)

			axs[1].quiver(clon, clat, np.cos(heading), np.sin(heading), angles='xy', scale = 7, color='#999999')
			axs[1].scatter(clon, clat, c='#377eb8', s=50, zorder=3)
			axs[1].scatter(seismo_longitudes[s], seismo_latitudes[s], c='#e41a1c', s=50, zorder=3)

			axs[1].text(xx,yy, str(round(dist_final, 3))+' km', fontsize=15, fontweight='bold')

			# Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
			con = mpatch.ConnectionPatch(xyA=(minl, minla), xyB=(maxl, minla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
			fig.add_artist(con)
			con = mpatch.ConnectionPatch(xyA=(minl, maxla), xyB=(maxl, maxla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
			fig.add_artist(con)
			plt.show()

		else:
			continue