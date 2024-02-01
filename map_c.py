import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib.patches as mpatch
from matplotlib.patches import Rectangle
import os
from obspy.geodetics import gps2dist_azimuth
from obspy.core import UTCDateTime
import datetime
from pathlib import Path

def make_base_dir(base_dir):
    base_dir = Path(base_dir)
    if not base_dir.exists():
        current_path = Path("/")
        for parent in base_dir.parts:
            current_path = current_path/parent
            if not current_path.exists():
                current_path.mkdir()

def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km


flight_num = [530342801,528485724,528473220,528407493,528293430]
time = [1551066051,1550172833,1550168070,1550165577,1550089044]
station = [1022,1272,1173,1283,1004]
day = [25,14,14,14,13]

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
for line in range(0,5):
    date = '201902'+str(day[line])
    flight = str(flight_num[line])
    print(flight)
    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + date + '_positions/' + date + '_' + flight + '.csv'
    #print(flight_file)
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    tm = flight_data['snapshot_id']

    for l in range(len(tm)):
        if str(tm[l]) == str(time[line]):
            for t  in range(len(sta)):
                if sta[t] == station[line]:
                    dist = distance(seismo_latitudes[t], seismo_longitudes[t], flight_latitudes[l], flight_longitudes[l])
                    ht = datetime.datetime.utcfromtimestamp(tm[l])
                    # Create a figure with two subplots side by side
                    fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

                    axs[0].set_aspect(f)
                    axs[1].set_aspect(f)

                    axs[0].scatter(seismo_longitudes, seismo_latitudes, c='red', s = 3, label='seismometers')
                    axs[0].plot(flight_longitudes, flight_latitudes, '-o', c='c', lw=1, ms = 1, label='flight path')

                    # Set labels and title
                    axs[0].set_xlim(min_lon, max_lon)
                    axs[0].set_ylim(min_lat, max_lat)
                    axs[0].tick_params(axis='both', which='major', labelsize=13)

			
                    y =[flight_latitudes[l],  seismo_latitudes[t]]
                    x = [flight_longitudes[l], seismo_longitudes[t]]
                    yy = sum(y)/len(y)
                    xx = sum(x)/len(x)

                    minl = xx - 0.05
                    maxl = xx + 0.05
                    minla = yy - 0.03
                    maxla = yy + 0.03

                    rect = Rectangle((minl, minla), 0.1, 0.06, ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
                    axs[0].add_patch(rect)
                    axs[1].plot(x,y, '--', c='orange')
                    # Draw the zoomed in map on the second subplot
                    axs[1].scatter(seismo_longitudes, seismo_latitudes, c='red')
                    axs[1].plot(flight_longitudes, flight_latitudes, c='c',linestyle ='dotted')
                    axs[1].set_xlim(minl, maxl)
                    axs[1].set_ylim(minla, maxla)
                    axs[1].text(seismo_longitudes[t], seismo_latitudes[t], sta[t], fontsize=9, fontweight='bold')
                    axs[1].tick_params(axis='both', which='major', labelsize=13)
                    axs[1].text(flight_longitudes[l], flight_latitudes[l], ht, fontsize=9, fontweight='bold')
                    axs[1].scatter(flight_longitudes[l], flight_latitudes[l], c='lawngreen')
                    axs[1].scatter(seismo_longitudes[t], seismo_latitudes[t], c='pink')

                    

                    axs[1].text(xx,yy, str(round(dist, 2))+'km', fontsize=8, fontweight='bold')

                    # Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
                    con = mpatch.ConnectionPatch(xyA=(minl, minla), xyB=(maxl, minla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                    con = mpatch.ConnectionPatch(xyA=(minl, maxla), xyB=(maxl, maxla), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                
                    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5map_all/' + date + '/'+flight+'/'+station+'/'
                    make_base_dir(BASE_DIR)
                    plt.savefig('/scratch/irseppi/nodal_data/plane_info/5map_all/'+ date + '/'+flight+'/'+station+'/map_'+flight+'_' + str(tm[l]) + '.png')
                    plt.close()
                    #print(sta[t], flight)
                else:
                    continue

        else:
            continue

import datetime
from math import radians, sin, cos, sqrt, atan2

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0  # radius of the Earth in kilometers

    lat1_rad = radians(lat1)
    lon1_rad = radians(lon1)
    lat2_rad = radians(lat2)
    lon2_rad = radians(lon2)

    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    a = sin(dlat / 2)**2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    return R * c  # returns distance in kilometers

def closest_encounter(aircraft_data, seismometer_location):
    closest_distance = float('inf')
    closest_time = None

    for timestamp, aircraft_location in aircraft_data:
        distance = haversine(aircraft_location[0], aircraft_location[1], seismometer_location[0], seismometer_location[1])
        if distance < closest_distance:
            closest_distance = distance
            closest_time = timestamp

    return closest_time, closest_distance

def calculate_wave_arrival(closest_time, closest_distance):
    speed_of_sound = 343  # speed of sound in m/s
    time_for_sound_to_travel = (closest_distance * 1000) / speed_of_sound  # convert distance to meters

    wave_arrival_time = closest_time + datetime.timedelta(seconds=time_for_sound_to_travel)

    return wave_arrival_time
from geopy.distance import gps2dist_azimuth



import datetime

def closest_encounter(aircraft_data, seismometer_location):
    closest_distance = float('inf')
    closest_time = None
    closest_altitude = None

    for timestamp, aircraft_location, altitude in aircraft_data:
        _, _, distance = gps2dist_azimuth(aircraft_location[0], aircraft_location[1], seismometer_location[0], seismometer_location[1])
        distance = distance / 1000  # convert distance to kilometers
        if distance < closest_distance:
            closest_distance = distance
            closest_time = timestamp
            closest_altitude = altitude

    return closest_time, closest_distance, closest_altitude

def calculate_wave_arrival(closest_time, closest_distance, closest_altitude):
    speed_of_sound = 343  # speed of sound in m/s at sea level
    time_for_sound_to_travel = sqrt((closest_distance * 1000)**2 + (closest_altitude * 0.3048)**2) / speed_of_sound  # convert distance to meters and altitude to meters

    wave_arrival_time = closest_time + datetime.timedelta(seconds=time_for_sound_to_travel)

    return wave_arrival_time


from math import radians, cos

def calculate_wave_arrival(closest_time, closest_distance, closest_altitude, aircraft_speed, aircraft_heading, seismometer_heading):
    speed_of_sound = 343  # speed of sound in m/s at sea level
    aircraft_speed_mps = aircraft_speed * 0.514  # convert aircraft speed from knots to m/s

    # calculate the component of the aircraft's speed in the direction of the seismometer
    relative_heading = radians(aircraft_heading - seismometer_heading)
    relative_speed = aircraft_speed_mps * cos(relative_heading)

    # calculate the time it takes for the sound to travel from the aircraft to the seismometer
    time_for_sound_to_travel = sqrt((closest_distance * 1000)**2 + (closest_altitude * 0.3048)**2) / (speed_of_sound + relative_speed)

    wave_arrival_time = closest_time + datetime.timedelta(seconds=time_for_sound_to_travel)

    return wave_arrival_time



