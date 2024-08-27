import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from prelude import distance, load_flights

# Load flight files and filenames
flight_files,filenames = load_flights(2, 4, 11, 27)

# Load the seismometer location data
seismo_data = pd.read_csv('input/nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']


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
def rf(x):
    return (-180*np.arccos(-1/x))/np.pi

# Iterate over flight files
for i, flight_file in enumerate(flight_files):
	print((i/len(flight_files))*100, '%')	# Print progress

	# Load flight data and extract relevant columns
	flight_data = pd.read_csv(flight_file, sep=",") 
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']


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
				dist = distance(flight_latitudes[index], flight_longitudes[index], seismo_latitudes[s], seismo_longitudes[s])
				if dist < 2.01:
					clat = flight_latitudes[index]
					clon = flight_longitudes[index]
					
					closest_lat = clat
					closest_lon = clon

					for tr in range(-1,2,2):
							
						sclat = flight_latitudes[index+tr]
						sclon = flight_longitudes[index+tr]  
						
						x = [clon, sclon]
						y = [clat, sclat]
						m = (y[1]-y[0])/(x[1]-x[0])
						b = y[0] - m*x[0]
						print(clat, clon)
						print(sclat, sclon)
						print(m, b)
						if m < 0:
							ggg = -0.0000001
						elif m > 0:
							ggg = 0.0000001
						for point in np.arange(clat+ggg, sclat+ggg, ggg):
							lat = point
							lon = (lat - b)/m
							dist_km = distance(lat, lon, seismo_latitudes[s], seismo_longitudes[s])
							
							if dist_km < dist_lim:
								dist_lim = dist_km
								closest_lat = lat
								closest_lon = lon
								c2lat = flight_latitudes[index+tr]
								c2lon = flight_longitudes[index+tr]
								index2 = index + tr
								index_final = index
							else:
								continue


		if dist_final > dist_lim:
			con = True
			dist_final = dist_lim

		else:
			continue

		if con == True:             
			clat = closest_lat
			clon = closest_lon

			fig, axs = plt.subplots(1, 1) 

			y =[clat,  seismo_latitudes[s]]
			x = [clon, seismo_longitudes[s]]

			fig.set_xscale('function', functions=(f,rf))

			fig.grid(True,linestyle='dotted',color='gray')

			fig.plot(x,y, '--', c='#ff7f00')
			fig.scatter(x,y, '--', c='#ff7f00')
			# Draw the zoomed in map on the second subplot
			fig.plot(flight_longitudes, flight_latitudes, c='#377eb8',linestyle ='--')
			fig.scatter(flight_longitudes, flight_latitudes, c='green')
			fig.scatter(flight_longitudes[index_final], flight_latitudes[index_final], c='pink')

			fig.scatter(clon, clat, c='#377eb8', s=50, zorder=3)
			fig.scatter(seismo_longitudes[s], seismo_latitudes[s], c='#e41a1c', s=50, zorder=3)

			plt.show()

		else:
			continue