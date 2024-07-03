import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.patches as mpatch
from obspy.geodetics import gps2dist_azimuth
from matplotlib.patches import Rectangle

sta_f = open('input/all_station_crossing_db.txt','r')

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
elevations = seismo_data['Elevation']
sta = seismo_data['Station']
'''

# Python program to find the distance between
# a given point and a given line in 2 D.
 
import math
 
# Function to find distance
def shortest_distance(x1, y1, a, b, c): 
      
    d = abs((a * x1 + b * y1 + c)) / (math.sqrt(a * a + b * b))
    print("Perpendicular distance is"),d
     
 
# Driver Code 
x1 = 5
y1 = 6
a = -2
b = 3
c = 4
shortest_distance(x1, y1, a, b, c)  
'''
def find_closest_distance(point1, point2, point3, fourth_point):
	# Sort the three points based on their timestamps
	sorted_points = sorted([point1, point2, point3], key=lambda p: p[2])

	# Calculate the equations of the two lines formed by connecting the three points
	line1 = calculate_line_equation(sorted_points[0], sorted_points[1])
	line2 = calculate_line_equation(sorted_points[1], sorted_points[2])

	# Find the closest point on line1 to the fourth point
	closest_point1 = find_closest_point(line1, fourth_point)

	# Find the closest point on line2 to the fourth point
	closest_point2 = find_closest_point(line2, fourth_point)

	# Calculate the distance between the closest points and the fourth point
	distance1 = calculate_distance(closest_point1[1], closest_point1[0], fourth_point[1], fourth_point[0]) 
	distance2 = calculate_distance(closest_point2[1], closest_point2[0], fourth_point[1], fourth_point[0]) 
	# Return the minimum distance
	#return min(distance1, distance2), 
	if distance1 < distance2:
		return distance1, closest_point1
	else:
		return distance2, closest_point2

def find_closest_point(line, point):
    # Calculate the coordinates of the closest point on the line to the given point
    x = (line['slope'] * (point[1] - line['intercept']) + point[0]) / (line['slope']**2 + 1)
    y = line['slope'] * x + line['intercept']
    
    return (x, y)

def calculate_line_equation(point1, point2):
    # Calculate the equation of a line given two points
    x1 = point1[0]
    y1 = point1[1]
    x2 = point2[0]
    y2 = point2[1]

    slope = (y2 - y1) / ((x2 - x1)*1.0/np.cos((y1+((y2 - y1)/2)*np.pi/180)))
    intercept = y1 - slope * x1

    return {'slope': slope, 'intercept': intercept}

def calculate_distance(lat1, lon1, lat2, lon2):
	# Calculate the distance between two points in space
    distance, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)

    return distance
# Loop through each station in text file that we already know comes within 2km of the nodes
for line in sta_f.readlines():
    val = line.split(',')
    date = val[0]
    
    flight = val[1]
    station = val[5]
    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + date + '_positions/' + date + '_' + flight + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    time = flight_data['snapshot_id']
    head = flight_data['heading']
    for l in range(len(time)):
        if str(time[l]) == str(val[2]):
            for t  in range(len(sta)):
                if sta[t] == station:

                    point1=[flight_longitudes[l], flight_latitudes[l], time[l]]
                    point2=[flight_longitudes[l+1], flight_latitudes[l+1], time[l+1]]
                    point3=[flight_longitudes[l-1], flight_latitudes[l+1], time[l-1]]
                    point4=[seismo_longitudes[t], seismo_latitudes[t]]

                    dist, closest_point = find_closest_distance(point1, point2, point3, point4)
                    clat = closest_point[1]
                    clon = closest_point[0]
                    f = 1.0/np.cos((clat+((clat-seismo_latitudes[t])/2))*np.pi/180)
                    ht = datetime.datetime.utcfromtimestamp(time[l])

                    # Create a figure with two subplots side by side
                    fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

                    axs[0].set_aspect(f)
                    axs[1].set_aspect(f)

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
			
                    y =[clat,  seismo_latitudes[t]]
                    x = [clon, seismo_longitudes[t]]
                    yy = sum(y)/len(y)
                    xx = sum(x)/len(x)

                    minl = xx - 0.05
                    maxl = xx + 0.05
                    minla = yy - 0.03
                    maxla = yy + 0.03
                    heading = np.deg2rad(head[l])
                    
                    rect = Rectangle((minl, minla), 0.1, 0.06, ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
                    axs[0].add_patch(rect)
                    axs[1].plot(x,y, '--', c='#ff7f00')

                    # Draw the zoomed in map on the second subplot
                    axs[1].plot(flight_longitudes, flight_latitudes, c='#377eb8',linestyle ='dotted')
                    axs[1].set_xlim(minl, maxl)
                    axs[1].set_ylim(minla, maxla)
                    #axs[1].text(seismo_longitudes[t], seismo_latitudes[t], sta[t], fontsize=11, fontweight='bold')
                    axs[1].tick_params(axis='both', which='major', labelsize=10)
                    
                    
                    #axs[1].text(flight_longitudes[l], flight_latitudes[l], ht, fontsize=11, fontweight='bold')
                    m = (flight_latitudes[l+1] - flight_latitudes[l])/(flight_longitudes[l+1] - flight_longitudes[l])
                    b = flight_latitudes[l] - m*flight_longitudes[l]

                    direction = np.arctan2(flight_latitudes[l+1] - flight_latitudes[l], flight_longitudes[l+1] - flight_longitudes[l])
                        
                    axs[1].quiver(clon, clat, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=10)

                    axs[1].quiver(clon, clat, np.cos(heading), np.sin(heading), angles='xy', scale = 10, color='#999999')
                    axs[1].scatter(clon, clat, c='#377eb8', s=50, zorder=3)
                    axs[1].scatter(seismo_longitudes[t], seismo_latitudes[t], c='#e41a1c', s=50, zorder=3)

                    axs[1].text(xx,yy, str(round(dist/1000, 2))+' km', fontsize=15, fontweight='bold')

                    # Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
                    con = mpatch.ConnectionPatch(xyA=(minl, minla), xyB=(maxl, minla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                    con = mpatch.ConnectionPatch(xyA=(minl, maxla), xyB=(maxl, maxla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                    plt.show()
               
                    
                else:
                    continue

        else:
            continue

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import datetime
import matplotlib.patches as mpatch

from matplotlib.patches import Rectangle
from pathlib import Path
from prelude import make_base_dir, distance, closest_encounter
from math import radians, cos, sin, asin, sqrt

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    return c * r

sta_f = open('input/all_station_crossing_db.txt','r')

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3
cenlat = (min_lat + max_lat)/2
cenlon = (min_lon + max_lon)/2
# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
sta = seismo_data['Station']

# Loop through each station in text file that we already know comes within 2km of the nodes
for line in sta_f.readlines():
    val = line.split(',')
    date = val[0]
    
    flight = val[1]
    station = val[5]
    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + date + '_positions/' + date + '_' + flight + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    time = flight_data['snapshot_id']
    head = flight_data['heading']
    for l in range(len(time)):
        if str(time[l]) == str(val[2]):
            for t  in range(len(sta)):
                if sta[t] == station:

                    clat, clon, dist, ctime = closest_encounter(flight_latitudes, flight_longitudes, l, time, seismo_latitudes[t], seismo_longitudes[t])
                    #dist = haversine(clon, clat, seismo_longitudes[t], seismo_latitudes[t])
                    f = 1.0/np.cos((clat-(seismo_latitudes[t]-clat)/2)*np.pi/180)
                    ht = datetime.datetime.utcfromtimestamp(time[l])
                
                    # Create a figure with two subplots side by side
                    fig = plt.figure(figsize=(10, 5))


                    axs0 = fig.add_subplot(1,2,1,projection=ccrs.AlbersEqualArea(central_longitude=cenlon,central_latitude=cenlat)) #projection=ccrs.UTM(zone=5))
                    axs1 = fig.add_subplot(1,2,2,projection=ccrs.AlbersEqualArea(central_longitude=clon,central_latitude=clat)) #projection=ccrs.UTM(zone=5))

                    axs0.set_extent([min_lon, max_lon, min_lat, max_lat], ccrs.AlbersEqualArea(central_longitude=cenlon,central_latitude=cenlat))

                    g0 = axs0.gridlines(draw_labels=True,dms=True)
                    g0.top_labels = False
                    g0.right_labels = False
                    g0.rotate_labels = False
                    print(1.0/np.cos((min_lat+cenlat)*np.pi/180))
                 
                    #axs0.set_aspect(1.0/np.cos((max_lat-cenlat)*np.pi/180))
                    #axs1.set_aspect(f)

                    axs0.scatter(seismo_longitudes, seismo_latitudes, c='#e41a1c', s = 3, label='seismometers')
                    axs0.plot(flight_longitudes, flight_latitudes, '-', c='#377eb8', lw=1, ms = 1, label='flight path')
                    for i in range(1, len(flight_latitudes)-1, int(len(flight_latitudes)/5)):
                        direction = np.arctan2(flight_latitudes[i+1] - flight_latitudes[i], flight_longitudes[i+1] - flight_longitudes[i])
                        m = (flight_latitudes[i+1] - flight_latitudes[i])/(flight_longitudes[i+1] - flight_longitudes[i])
                        b = flight_latitudes[i] - m*flight_longitudes[i]
                        axs0.quiver((flight_latitudes[i]-b)/m, flight_latitudes[i], np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', headwidth = 5)
                    # Set labels and title
                    #axs0.set_xlim(min_lon, max_lon)
                    #axs0.set_ylim(min_lat, max_lat)
                    axs0.tick_params(axis='both', which='major', labelsize=12)
			
                    y =[clat,  seismo_latitudes[t]]
                    x = [clon, seismo_longitudes[t]]
                    yy = sum(y)/len(y)
                    xx = sum(x)/len(x)

                    minl = xx - 0.03
                    maxl = xx + 0.03
                    minla = yy - 0.03
                    maxla = yy + 0.03
                    axs1.set_extent([minl, maxl, minla, maxla], ccrs.AlbersEqualArea(central_longitude=clon,central_latitude=clat))
                    heading = np.deg2rad(head[l])
                    
                    rect = Rectangle((minl, minla), 0.1, 0.06, ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
                    axs0.add_patch(rect)
                    axs1.plot(x,y, '--', c='#ff7f00')

                    # Draw the zoomed in map on the second subplot
                    axs1.plot(flight_longitudes, flight_latitudes, c='#377eb8',linestyle ='dotted')
                    #axs1.set_xlim(minl, maxl)
                    #axs1.set_ylim(minla, maxla)
                    #axs1.tick_params(axis='both', which='major', labelsize=10)
                    
                     
                    m = (flight_latitudes[l+1] - flight_latitudes[l])/(flight_longitudes[l+1] - flight_longitudes[l])
                    b = flight_latitudes[l] - m*flight_longitudes[l]

                    direction = np.arctan2(flight_latitudes[l+1] - flight_latitudes[l], flight_longitudes[l+1] - flight_longitudes[l])
                        
                    #axs1.quiver(clon, clat, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=10)

                    #axs1.quiver(clon, clat, np.cos(heading), np.sin(heading), angles='xy', scale = 10, color='#999999')
                    axs1.scatter(clon, clat, c='#377eb8', s=50, zorder=3)
                    axs1.scatter(seismo_longitudes[t], seismo_latitudes[t], c='#e41a1c', s=50, zorder=3)

                    axs1.text(xx,yy, str(round(dist/1000, 2))+' km', fontsize=15, fontweight='bold')

                    g1 = axs1.gridlines(draw_labels=True,dms=True)
                    g1.rotate_labels = False                 
                    g1.top_labels = False
                    g1.right_labels = False

                    # Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
                    con = mpatch.ConnectionPatch(xyA=(minl, minla), xyB=(maxl, minla), coordsA="data", coordsB="data", axesA=axs1, axesB=axs0, color="black", linestyle="--")
                    fig.add_artist(con)
                    con = mpatch.ConnectionPatch(xyA=(minl, maxla), xyB=(maxl, maxla), coordsA="data", coordsB="data", axesA=axs1, axesB=axs0, color="black", linestyle="--")
                    fig.add_artist(con)
                    plt.show()
               
                    
                else:
                    continue

        else:
            continue

