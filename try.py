import numpy as np
import pandas as pd
import obspy
from datetime import datetime, timezone
from pyproj import Proj
from prelude import *
from scipy.signal import spectrogram
from matplotlib import pyplot as plt
from plot_func import remove_median

equip = 'C185'
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')
sta_f = open('input/all_station_crossing_db_' + equip + '.txt','r')
second_column = []
for line in sta_f.readlines():
    val = line.split(',')
    if len(val) >= 2:
        second_column.append(val[1])
sta_f.close()
second_column_array = np.array(second_column)

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt','r')
for li in file_in.readlines():
    text = li.split(',')
    flight_num = text[1]
    if flight_num not in second_column_array:
        continue
    date = text[0]
    try:
        sta = int(text[9])
    except:
        continue
    time = float(text[5])

    ht = datetime.fromtimestamp(time, tz=timezone.utc)
    h = ht.hour

    alt = float(text[4])*0.0003048 #convert between feet and km
    x =  float(text[2])  # Replace with your UTM x-coordinate
    y = float(text[3])  # Replace with your UTM y-coordinate

    Tc = -2
    c = speed_of_sound(Tc)
    sound_speed = c

    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight_num) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
 
    timestamps = flight_data['snapshot_id']
    speed = flight_data['speed']
    altitude = flight_data['altitude']

    closest_x, closest_y, dist_km, closest_time, tarrive, alt, sp, elevation, speed_mps, height_m, dist_m, tmid = closest_approach_UTM(seismo_latitudes, seismo_longitudes, flight_latitudes, flight_longitudes, timestamps, altitude, speed, stations, elevations, c, sta)
    if closest_x == None:
        continue
    #tarrive is different each time, need a standard base time that in not changed by temperature and speed of sound
    ht = datetime.fromtimestamp(tarrive, tz=timezone.utc)
    mins = ht.minute
    secs = ht.second
    month = ht.month
    day = ht.day

    h_u = str(h+1)
    if h < 23:			
        day2 = str(day)
        if h < 10:
            h_u = '0'+str(h+1)
            h = '0'+str(h)
        else:
            h_u = str(h+1)
            h = str(h)
    else:
        h_u = '00'
        day2 = str(day+1)
    if len(str(day)) == 1:
        day = '0'+str(day)
        day2 = day

    try:
        p = "/scratch/naalexeev/NODAL/2019-0"+str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month)+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(sta)+".mseed"
        tr = obspy.read(p)
    except:
        continue

    tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - 120, tr[2].stats.starttime + (mins * 60) + secs + 120)
    data = tr[2][:]
    fs = int(tr[2].stats.sampling_rate)
    title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						


    # Compute spectrogram
    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 

    spec, MDF = remove_median(Sxx)
    
    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 
 

    file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'
                
    if Path(file_name).exists():
        coords = []
        with open(file_name, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                coords.append((float(pick_data[0]), float(pick_data[1])))
        file.close()  

    else:
        continue

    coords_array = np.array(coords) 

    output2 = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/overtonepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight_num)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight_num)+'.csv'
    if Path(output2).exists():

        peaks = []
        freqpeak = []
        with open(output2, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                peaks.append(float(pick_data[1]))
                freqpeak.append(float(pick_data[0]))
        file.close()  

    else:
        continue

    plt.figure()
    plt.title(title)
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=0, vmax=vmax)
    plt.plot(coords_array[:,0], coords_array[:,1], 'x', color='red')
    plt.plot(freqpeak, peaks, 'x', color='red')
    plt.show()

