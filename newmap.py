import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.patches as mpatch
from obspy.geodetics import gps2dist_azimuth
from matplotlib.patches import Rectangle
from pathlib import Path
from prelude import make_base_dir, distance, closest_encounter

sta_f = open('input/all_station_crossing_db.txt','r')
def calculate_distance(lat1, lon1, lat2, lon2):

	distance, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)  # distance in meters
	return distance

def closest_encounter(flight_latitudes, flight_longitudes, index, timestamp, seismo_latitude, seismo_longitude):

    clat = flight_latitudes[index]
    clon = flight_longitudes[index]
        
    closest_lat = clat
    closest_lon = clon
    dist_lim = 2
    c2lat = flight_latitudes[index+1]
    c2lon = flight_longitudes[index+1]
    index2 = index + 1
    for tr in range(0,2):
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
                    dist_km = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000

                    if dist_km < dist_lim:
                        dist_lim = dist_km
                        closest_lat = lat
                        closest_lon = lon
                        c2lat = flight_latitudes[index-1]
                        c2lon = flight_longitudes[index-1]
                        index2 = index - 1
            if tr == 1:
                sclat = flight_latitudes[index+1]
                sclon = flight_longitudes[index+1]

                x = [clon, sclon]
                y = [clat, sclat]
                m = (y[1]-y[0])/(x[1]-x[0])
                b = y[0] - m*x[0]


                for point in np.arange(sclon, clon, 0.000001):
                    lon = point
                    lat = m*lon + b
                    dist_km = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000

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
                    dist_km = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000

                    if dist_km < dist_lim:
                        dist_lim = dist_km
                        closest_lat = lat
                        closest_lon = lon
                        c2lat = flight_latitudes[index-1]
                        c2lon = flight_longitudes[index-1]
                        index2 = index - 1
            if tr == 1:
                sclat = flight_latitudes[index+1]
                sclon = flight_longitudes[index+1]

                x = [clon, sclon]
                y = [clat, sclat]
                m = (y[1]-y[0])/(x[1]-x[0])
                b = y[0] - m*x[0]


                for point in np.arange(clon, sclon, 0.000001):
                    lon = point
                    lat = m*lon + b
                    dist_km = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000

                    if dist_km < dist_lim:
                        dist_lim = dist_km
                        closest_lat = lat
                        closest_lon = lon
                        c2lat = flight_latitudes[index+1]
                        c2lon = flight_longitudes[index+1]
                        index2 = index - 1
    for location in np.arange((closest_lon-0.000001),(closest_lon+0.000001),0.000000000000001):
        lon = location
        lat = m*lon + b
        dist = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000
        if dist < dist_lim:
            dist_lim = dist
            closest_lon = lon
            closest_lat = lat
    dist_old_new = calculate_distance(closest_lat, closest_lon, clat, clon)/1000
    dist_old_old = calculate_distance(c2lat, c2lon, clat, clon)/1000
    ratio = dist_old_new/dist_old_old
    timestamp = timestamp[index] + (timestamp[index] - timestamp[index2])*ratio
    return  closest_lat, closest_lon, dist_lim,timestamp

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
def f(y):
    return -1.0/np.cos((y)*np.pi/180)
def rf(y):
    return -180/np.pi*np.arccos(1/y)
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
                    ht = datetime.datetime.fromtimestamp(time[l], datetime.timezone.utc)
                    print(dist)
                    dist = dist*1000
                    # Create a figure with two subplots side by side
                    fig, axs = plt.subplots(1, 2,constrained_layout=True)

                    y =[clat,  seismo_latitudes[t]]
                    x = [clon, seismo_longitudes[t]]
                    yy = sum(y)/len(y)
                    xx = sum(x)/len(x)
                    if dist < 100:
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
                    axs[1].tick_params(axis='both', which='major', labelsize=10)
                    
                     
                    m = (flight_latitudes[l+1] - flight_latitudes[l])/(flight_longitudes[l+1] - flight_longitudes[l])
                    b = flight_latitudes[l] - m*flight_longitudes[l]

                    direction = np.arctan2(flight_latitudes[l+1] - flight_latitudes[l], flight_longitudes[l+1] - flight_longitudes[l])
                        
                    axs[1].quiver(clon, clat, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=7)

                    axs[1].quiver(clon, clat, np.cos(heading), np.sin(heading), angles='xy', scale = 7, color='#999999')
                    axs[1].scatter(clon, clat, c='#377eb8', s=50, zorder=3)
                    axs[1].scatter(seismo_longitudes[t], seismo_latitudes[t], c='#e41a1c', s=50, zorder=3)

                    axs[1].text(xx,yy, str(round(dist/1000, 3))+' km', fontsize=15, fontweight='bold')

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

'''

# Create a figure with two subplots side by side
fig = plt.figure() 

min_lon = -153
max_lon = -146
min_lat = 62.2
max_lat = 65.3
cenlat = (min_lat + max_lat)/2
cenlon = (min_lon + max_lon)/2
def f(y):
return -1.0/np.cos((y)*np.pi/180)
def rf(y):
return -180/np.pi*np.arccos(1/y)


axs0 = fig.add_subplot(1,2,1,projection=ccrs.UTM(zone=5)) 
axs1 = fig.add_subplot(1,2,2,projection=ccrs.UTM(zone=5))

axs0.set_extent([min_lon, max_lon, min_lat, max_lat]) #, ccrs.UTM(zone=5))
axs0.set_xscale('function', functions=(f,rf))
axs1.set_xscale('function', functions=(f,rf))
g0 = axs0.gridlines(draw_labels=True,dms=True)
g0.top_labels = False
g0.right_labels = False
g0.rotate_labels = False


y =[clat,  seismo_latitudes[t]]
x = [clon, seismo_longitudes[t]]
yy = sum(y)/len(y)
xx = sum(x)/len(x)

minl = xx - 0.02
maxl = xx + 0.02
minla = yy - 0.02
maxla = yy + 0.02
axs1.set_extent([minl, maxl, minla, maxla]) 

g1 = axs1.gridlines(draw_labels=True,dms=True)
g1.rotate_labels = False                 
g1.top_labels = False
g1.right_labels = False


'''