import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.patches as mpatch

from matplotlib.patches import Rectangle
from pathlib import Path
from prelude import make_base_dir, distance

sta_f = open('input/all_station_crossing_db.txt','r')

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3
f = 1.0/np.cos(62*np.pi/180)

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
sta = seismo_data['Station']

# Loop through each station in text file that we already know comes within 2km of the nodes
for line in sta_f.readlines():
    val = line.split(',')
    date = val[0]
    print(date)
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
                    dist = distance(seismo_latitudes[t], seismo_longitudes[t], flight_latitudes[l], flight_longitudes[l])
                    ht = datetime.datetime.utcfromtimestamp(time[l])
                    # Create a figure with two subplots side by side
                    fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

                    axs[0].set_aspect(f)
                    axs[1].set_aspect(f)

                    axs[0].scatter(seismo_longitudes, seismo_latitudes, c='#e41a1c', s = 3, label='seismometers')
                    axs[0].plot(flight_longitudes, flight_latitudes, '-', c='#377eb8', lw=1, ms = 1, label='flight path')
                    for i in range(int(len(flight_latitudes)/5), len(flight_latitudes)-2, int(len(flight_latitudes)/5)):
                        direction = np.arctan2(flight_latitudes[i+1] - flight_latitudes[i], flight_longitudes[i+1] - flight_longitudes[i])
                        axs[0].quiver(flight_longitudes[i], flight_latitudes[i], np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', headwidth = 10, headlength = 7)
                    # Set labels and title
                    axs[0].set_xlim(min_lon, max_lon)
                    axs[0].set_ylim(min_lat, max_lat)
                    axs[0].tick_params(axis='both', which='major', labelsize=12)

			
                    y =[flight_latitudes[l],  seismo_latitudes[t]]
                    x = [flight_longitudes[l], seismo_longitudes[t]]
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
                    axs[1].plot(flight_longitudes, flight_latitudes, c='#377eb8',linestyle ='dotted', label='Flight Path')
                    axs[1].set_xlim(minl, maxl)
                    axs[1].set_ylim(minla, maxla)
                    axs[1].tick_params(axis='both', which='major', labelsize=10)
                    axs[1].quiver(flight_longitudes[l], flight_latitudes[l], np.cos(heading), np.sin(heading), angles='xy', color='#999999', scale = 8, headwidth = 5, headlength = 6,label = 'Heading Direction')
                    if l < len(time)-1:
                        direction = np.arctan2(flight_latitudes[l+1] - flight_latitudes[l], flight_longitudes[l+1] - flight_longitudes[l])
                        axs[1].quiver(flight_longitudes[l], flight_latitudes[l], np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=8, headwidth = 5, headlength = 6, label = 'Direction of Motion')
                    
                    axs[1].scatter(flight_longitudes[l], flight_latitudes[l], c='#377eb8', s=50, zorder=3, label='Timestamp: ' + str(ht))
                    axs[1].scatter(seismo_longitudes[t], seismo_latitudes[t], c='#e41a1c', s=50, zorder=3, label='Station '+str(sta[t]))
                    axs[1].legend(loc='lower right', fontsize=8)
                    axs[1].text(xx, yy, str(round(dist, 2))+' km', fontsize=15, fontweight='bold')

                    # Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
                    con = mpatch.ConnectionPatch(xyA=(minl, minla), xyB=(maxl, minla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                    con = mpatch.ConnectionPatch(xyA=(minl, maxla), xyB=(maxl, maxla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                    
                    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/map_all/' + date + '/'+flight+'/'+station+'/'
                    make_base_dir(BASE_DIR)
                    plt.savefig('/scratch/irseppi/nodal_data/plane_info/map_all/'+ date + '/'+flight+'/'+station+'/map_'+flight+'_' + str(time[l]) + '.png')
                    plt.close()
                    
                else:
                    continue

        else:
            continue


