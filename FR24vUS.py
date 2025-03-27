import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from pyproj import Proj

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt','r')

latc = []
lonc = []
altc = []
times_list = []
speeds_list = []
dists_list = []

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
        for i in range(len(latc)):
            if times_list[i] == time:
                index_UTC = i
                lat = latc[i]
                lon = lonc[i]
                alt = altc[i]
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
        temp_c.append(Tc)
        time_relative = False
        if time_relative:
            time_new.append(float(lines[4]))
            times_org.append(120)
        elif not time_relative:
            time_new.append((float(lines[3])-120) + float(lines[4]))
            times_org.append(times_list[index_UTC])

        v0_new.append(float(lines[5]))
        distance_new.append(float(lines[6]))

        speeds_org.append(speeds_list[index_UTC])
        dists_org.append(dists_list[index_UTC])
        date.append(y)


    scatter1 = axs[idx, 0].scatter(v0_new, speeds_org, c=temp_c, cmap='coolwarm')
    axs[idx, 0].set_title(f"{title[idx]}: Velocity", fontsize=10)
    axs[idx, 0].set_xlim(40, 80)
    axs[idx, 0].set_ylim(40, 80)
    #axs[idx, 0].set_xlabel("Velocity New (m/s)", fontsize=8)
    #axs[idx, 0].set_ylabel("Speeds Org (m/s)", fontsize=8)

    scatter2 = axs[idx, 1].scatter(distance_new, dists_org, c=temp_c, cmap='coolwarm')
    axs[idx, 1].set_title(f"{title[idx]}: Distance", fontsize=10)
    axs[idx, 1].set_xlim(0, 3000)
    axs[idx, 1].set_ylim(0, 2000)
    #axs[idx, 1].set_xlabel("Distance New (m)", fontsize=8)
    #axs[idx, 1].set_ylabel("Dists Org (m)", fontsize=8)
    tryT = 'no'
    if tryT == 'yes':
        scatter3 = axs[idx, 2].scatter(np.array(time_new) - (1.55*(10**9)), np.array(times_org) - (1.55*(10**9)), c=temp_c, cmap='coolwarm')
        axs[idx, 2].set_title(f"{title[idx]}: Time", fontsize=10)
        axs[idx, 2].set_xscale('log')
        axs[idx, 2].set_yscale('log')
    elif tryT == 'no':
        scatter3 = axs[idx, 2].scatter(np.array(time_new) - np.array(times_org), date, c=temp_c, cmap='coolwarm')
        axs[idx, 2].set_title(f"{title[idx]}: Diff in Time", fontsize=10)
        axs[idx, 2].set_xlim(-10.5, 0)
    else:
        scatter3 = axs[idx, 2].scatter(np.array(time_new), np.array(times_org), c=temp_c, cmap='coolwarm')
        axs[idx, 2].set_title(f"{title[idx]}: Time", fontsize=10)
    #axs[idx, 2].set_xlabel("Time New (s)", fontsize=8)
    #axs[idx, 2].set_ylabel("Times Org  (s)", fontsize=8)

    # Add a single colorbar for the entire figure
    cbar = fig.colorbar(scatter1, ax=axs[idx, 2], orientation='vertical', pad=0.1)
    cbar.set_label('Temperature (°C)')



plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.subplots_adjust(hspace=0.3, wspace=0.4)
plt.show()
