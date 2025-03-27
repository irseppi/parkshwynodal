import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from pyproj import Proj

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt','r')
latc = []
lonc = []
timec = []
altc = []
for line in file_in.readlines():
    text = line.split(',')

    timec.append(float(text[5]))
    x =  float(text[2])  # Replace with your UTM x-coordinate
    y = float(text[3])  # Replace with your UTM y-coordinate

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)
    latc.append(lat)
    lonc.append(lon)
    altc.append(float(text[4])*0.0003048) #convert between feet and km
file_in.close()

file_list = ['1o_atmc_v_2c.txt','2c_1o_v_full.txt','full_atmc_v_2c.txt','atmc_1o_v_full.txt']
title = ['One Harmonic/Fixed Temp v One Harmonic/Varying Temp','One Harmonic/Fixed Temp v Full Harmonics/Fixed Temp','Full Harmonics/Fixed Temp v Full Harmonics/Varying Temp','One Harmonic/Varying Temp v Full Harmonics/Varying Temp']
for gg, file in enumerate(file_list):
    file = open(file, 'r') 
    date = []
    y=0
    time_new = []
    time_old = []

    v0_old = []
    v0_new = []

    distance_old = []
    distance_new = []

    delf0_new = []
    delf0_old= []

    pppp_old = []
    pppp_new = []
    med_old = []
    med_new = []
    date_old = []
    date_new = []
    date_all = []
    temp_c = []
    for line in file.readlines():
        y += 1
        counts = []

        lines = line.split(',')
        time = float(lines[3])
        for i in range(len(latc)):
            if timec[i] == time:
                lat = latc[i]
                lon = lonc[i]
                alt = altc[i]
                break

        input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat'
        try:
            file =  open(input_files, 'r') #as file:
        except:
            print('No file for: ', date, flight_num)
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

        time_old.append(float(lines[4]))
        time_new.append(float(lines[9]))
        v0_old.append(float(lines[5]))
        v0_new.append(float(lines[10]))
        flight_num = float(lines[1])
        distance_old.append(float(lines[6]))
        distance_new.append(float(lines[11]))
        date.append(y)

        peaks_old = np.array(lines[7])

        peaks_old = str(peaks_old)
        peaks_old = np.char.replace(peaks_old, '[', '')
        peaks_old = np.char.replace(peaks_old, ']', '')

        peaks_old = str(peaks_old)
        peaks_old = np.array(peaks_old.split(' '))

        peaks_new = np.array(lines[12])

        peaks_new= str(peaks_new)
        peaks_new = np.char.replace(peaks_new, '[', '')
        peaks_new = np.char.replace(peaks_new, ']', '')

        peaks_new = str(peaks_new)
        peaks_new = np.array(peaks_new.split(' '))

        f1 = []
        for i in range(len(peaks_old)):
            peak_old = peaks_old[i]
            pppp_old.append(float(peak_old))
            date_old.append(y)
            if i == 0:
                continue
            diff = float(peaks_old[i]) - float(peaks_old[i-1])
            if diff > 22 or diff < 18:
                continue
            f1.append(diff)
        med_old.append(np.nanmedian(f1))
        f2 = []
        for i in range(len(peaks_new)):
            peak_new = peaks_new[i]
            pppp_new.append(float(peak_new))
            date_new.append(y)
            if i == 0:
                continue
            diff = float(peaks_new[i]) - float(peaks_new[i-1])
            if diff > 22 or diff < 18:
                continue
            f2.append(diff)
        med_new.append(np.nanmedian(f2))
        date_all.append(y)
    fig, ax1 = plt.subplots(1, 1, sharex=False, figsize=(8, 6))
    fig.suptitle(title[gg], fontsize=16)

    ax1.margins(x=0)
    ax2 = fig.add_axes([0.83, 0.11, 0.07, 0.77], sharey=ax1)

    ax1.scatter(pppp_old, date_old, c='b')
    ax1.scatter(pppp_new, date_new,c='r')

    ax2.scatter(med_old, date_all, c='b')
    ax2.scatter(med_new, date_all, c='r')

    ax1.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
    ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)

    ax1.set_xlabel('Frequency')
    ax2.set_xlabel('\u0394'+'F')
    ax1.legend(loc='upper left',fontsize = 'x-small')
    ax1.set_xlim(0, 300)
    ax1.set_xticks(range(0, 251, 25)) 
    
    plt.show()

    fig, axs = plt.subplots(2, 3, figsize=(8, 10), sharey=True)

    # First subplot
    axs[0, 0].scatter(v0_old, date, c='b', label='v0_old')
    axs[0, 0].scatter(v0_new, date, c='r', label='v0_new')
    axs[0, 0].set_title('v0')
    axs[0, 0].legend()
    axs[0, 0].set_ylabel('Index')

    # Second subplot
    scatter = axs[1,0].scatter((np.array(v0_new) - np.array(v0_old)), date, c=temp_c, cmap='coolwarm', label='Velocity Residuals')
    axs[1, 0].set_title("Velocity Residuals")

    axs[0, 0].set_ylabel('Index')
    # Third subplot
    axs[0, 1].scatter(distance_old, date, c='b', label='distance_old')
    axs[0, 1].scatter(distance_new, date, c='r', label='distance_new')
    axs[0, 1].set_title('Distance')
    axs[0, 1].legend()


    # Fourth subplot
    scatter = axs[1, 1].scatter((np.array(distance_new) - np.array(distance_old)), date, c=temp_c, cmap='coolwarm', label='Distance Residuals')
    axs[1, 1].set_title("Distance Residuals")

    # Fifth subplot
    axs[0, 2].scatter(time_old, date, c='b', label='time_old')
    axs[0, 2].scatter(time_new, date, c='r', label='time_new')
    axs[0, 2].set_title('Time')
    axs[0, 2].legend()

    # Sixth subplot
    scatter = axs[1, 2].scatter((np.array(time_new) - np.array(time_old)), date, c=temp_c, cmap='coolwarm', label='Time Residuals')
    axs[1, 2].set_title("Time Residuals")
    fig.colorbar(scatter, ax=axs[1, 2], orientation='vertical', label='Temperature (°C)')


    plt.tight_layout()
    plt.show()