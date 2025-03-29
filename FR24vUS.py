import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from pyproj import Proj
from prelude import *

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt','r')

latc = []
lonc = []
altc = []
times_list = []
speeds_list = []
dists_list = []
stat_list = []
for line in file_in.readlines():
    text = line.split(',')
    x =  float(text[2])  
    y = float(text[3])  

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)
    latc.append(lat)
    lonc.append(lon)

    altc.append(float(text[4])*0.0003048) #convert between feet and km
    times_list.append(float(text[5]))
    speeds_list.append(float(text[7]))
    dists_list.append(float(text[4]))
    stat_list.append(text[9])
file_in.close()

file_list = ['C185data_1o.txt','C185data_atm_1o.txt','C185data_full.txt','C185data_atm_full.txt']
title = ['One Harmonic/Fixed Atmosphere Correction','One Harmonic/Time Varying Atmosphere Correction','Full Harmonics/Fixed Atmosphere Correction','Full Harmonics/Time Varying Atmosphere Correction']

fig, axs = plt.subplots(4, 3, figsize=(18, 24), sharey=False)
fig.suptitle("Flight Radar (Y-axis) vs Nodal Data Inversion (X_axis)")

for idx, fil in enumerate(file_list):
    data = open(fil, 'r')
 
    time_new = []
    v0_new = []
    distance_new = []

    times_org = []
    speeds_org = []
    dists_org = []
    
    date = []
    y=0
    temp_c = []

    for line in data.readlines():
        y += 1
        counts = []
        lines = line.split(',')
        time = float(lines[3])
        flight_num = lines[1]
        date_lab = lines[0]
        for i in range(len(latc)):
            if times_list[i] == time:
                index_UTC = i
                lat = latc[i]
                lon = lonc[i]
                alt = altc[i]
                sta = stat_list[i]
                break

        input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat'

        try:
            file =  open(input_files, 'r') 
        except:
            continue
        data = json.load(file)

        # Extract metadata
        metadata = data['metadata']
        sourcefile = metadata['sourcefile']
        datetim = metadata['time']['datetime']
        latitude = metadata['location']['latitude']
        longitude = metadata['location']['longitude']
        parameters = metadata['parameters']

        # Extract data
        data_list = data['data']

        # Convert data to a DataFrame
        data_frame = pd.DataFrame(data_list)

        # Find the "Z" parameter and extract the value at index
        z_index = None
        hold = np.inf
        for item in data_list:
            if item['parameter'] == 'Z':
                for i in range(len(item['values'])):
                    if abs(float(item['values'][i]) - float(alt)) < hold:
                        hold = abs(float(item['values'][i]) - float(alt))
                        z_index = i
        for item in data_list:
            if item['parameter'] == 'T':
                Tc = - 273.15 + float(item['values'][z_index])
        c = speed_of_sound(Tc)

        flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date_lab) + '_positions/' + str(date_lab) + '_' + str(flight_num) + '.csv'
        flight_data = pd.read_csv(flight_file, sep=",")
        flight_latitudes = flight_data['latitude']
        flight_longitudes = flight_data['longitude']
        time = flight_data['snapshot_id']
        timestamps = flight_data['snapshot_id']
        speed = flight_data['speed']
        altitude = flight_data['altitude']

        closest_x, closest_y, dist_km, closest_time, tarrive, alt, sp, elevation, speed_mps, height_m, dist_m, tmid = closest_approach_UTM(seismo_latitudes, seismo_longitudes, flight_latitudes, flight_longitudes, timestamps, altitude, speed, stations, elevations, c, sta)
        if closest_x == None:
            continue
        
        #To set the initial window of arrival correct picks your start end Must use the tarrive time to get the correct data
        ta_old = calc_time(tmid,dist_m,height_m,343)

        temp_c.append(Tc)
        time_relative = True
        if time_relative:
            time_new.append(float(lines[4]))
            times_org.append(tarrive - ta_old)
        else:
            time_new.append(float(lines[3]))
            times_org.append(tarrive)

        v0_new.append(float(lines[5]))
        distance_new.append(float(lines[6]))

        speeds_org.append(speeds_list[index_UTC])
        dists_org.append(dists_list[index_UTC])
        date.append(y)
    scatter1 = axs[idx, 0].scatter(v0_new, speeds_org, c=temp_c, cmap='coolwarm')
    axs[idx, 0].set_title(f"{title[idx]}: Velocity", fontsize=10)
    axs[idx, 0].set_xlim(50, 80)
    axs[idx, 0].axline((0, 0), slope=1, color='black', linestyle='--')
    axs[idx, 0].set_ylim(50, 80)
    axs[idx, 0].set_aspect('equal')
    axs[idx, 0].set_xticks(np.arange(50, 81, 10))
    axs[idx, 0].set_yticks(np.arange(50, 81, 10))

    scatter2 = axs[idx, 1].scatter(distance_new, dists_org, c=temp_c, cmap='coolwarm')
    axs[idx, 1].set_title(f"{title[idx]}: Distance", fontsize=10)
    axs[idx, 1].set_xlim(0, 2000)
    axs[idx, 1].set_ylim(0, 2000)
    axs[idx, 1].axline((0, 0), slope=1, color='black', linestyle='--')
    axs[idx, 1].set_aspect('equal', adjustable='box')
    axs[idx, 1].set_xticks(np.arange(0, 2001, 1000))
    axs[idx, 1].set_yticks(np.arange(0, 2001, 1000))

    if time_relative:
        scatter3 = axs[idx, 2].scatter(np.array(time_new), np.array(times_org), c=temp_c, cmap='coolwarm')
        axs[idx, 2].set_title(f"{title[idx]}: Time", fontsize=10)
    else:
        scatter3 = axs[idx, 2].scatter(np.array(time_new), np.array(times_org), c=temp_c, cmap='coolwarm')
        axs[idx, 2].set_title(f"{title[idx]}: Time", fontsize=10)
        axs[idx, 2].set_xscale('log')
        axs[idx, 2].set_yscale('log')
    axs[idx, 2].set_aspect('equal', adjustable='box')

    # Add a single colorbar for the entire figure
    cbar = fig.colorbar(scatter1, ax=axs[idx, 2], orientation='vertical', pad=0.1)
    cbar.set_label('Temperature (°C)')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.subplots_adjust(hspace=0.3, wspace=0.4)
plt.show()
