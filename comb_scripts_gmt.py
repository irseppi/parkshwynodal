import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from matplotlib.patches import Rectangle
from prelude import distance, load_flights
import pygmt

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
		dist_final = 2.01
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
			if dist_final > dist_lim:
				con = True
				dist_final = dist_lim
				print(dist_final)
				print(seismo_sta[s])
				print(flight_num)
			else:
				continue
		if con == True:     
			clat = closest_lat
			clon = closest_lon
 

			y =[clat,  seismo_latitudes[s]]
			x = [clon, seismo_longitudes[s]]
			yy = sum(y)/len(y)
			xx = sum(x)/len(x)
			if dist_final < 0.01:
				minl = np.round((xx - 0.001),4)
				maxl = np.round((xx + 0.001),4)
				minla = np.round((yy - 0.001),4)
				maxla =  np.round((yy + 0.001),4)
				
			else:
				minl = np.round((xx - 0.025), 2)
				maxl = np.round((xx + 0.025), 2)
				minla = np.round((yy - 0.025), 2)
				maxla =  np.round((yy + 0.025), 2)


			fig = pygmt.Figure()
			fig.basemap(region=[minl, maxl, minla, maxla], projection="Eclon/clat/c", frame=True)
			fig.coast(
				region=[minl, maxl, minla, maxla],
				shorelines=True,
				land="lightgreen",
				water="lightblue",
			)

	
			fig.plot(x=x, y=y) 

			# Draw the zoomed in map on the second subplot
			fig.plot(x=flight_longitudes, y=flight_latitudes)


			fig.plot(x=clon, y=clat, style="c0.2c", fill="blue", pen="black")
			fig.plot(x=seismo_longitudes[s], y=seismo_latitudes[s],  style="c0.2c", fill="red", pen="black")

			

			fig.show()

		else:
			continue