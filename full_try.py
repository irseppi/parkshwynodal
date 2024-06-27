import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import datetime
from prelude import *
from scipy.signal import find_peaks, spectrogram

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station = seismo_data['Station']
flight_num = [530342801,528485724,528473220,528407493,528293430] 
time = [1551066051,1550172833,1550168070,1550165577,1550089044] 
sta = [1022,1272,1173,1283,1004]
day = [25,14,14,14,13]
month = [2,2,2,2,2]

for n in range(0,5):
    ht = datetime.utcfromtimestamp(time[n])
    mins = ht.minute
    secs = ht.second
    h = ht.hour
    tim = 120	
    h_u = str(h+1)
    if h < 23:			
        day2 = str(day[n])
        if h < 10:
            h_u = '0'+str(h+1)
            h = '0'+str(h)
        else:
            h_u = str(h+1)
            h = str(h)
    else:
        h_u = '00'
        day2 = str(day[n]+1)
    if len(str(day[n])) == 1:
        day[n] = '0'+str(day[n])
        day2 = day[n]
    flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/20190'+str(month[n])+str(day[n])+'_positions/20190'+str(month[n])+str(day[n])+'_'+str(flight_num[n])+'.csv', sep=",")

    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    tm = flight_data['snapshot_id']
    speed = flight_data['speed']
    alt = flight_data['altitude']
    head = flight_data['heading']

    for line in range(len(tm)):
        if str(tm[line]) == str(time[n]):
            speed = flight_data['speed'][line]
            speed_mps = speed * 0.514444
            alt = flight_data['altitude'][line]
            alt_m = alt * 0.3048

            for y in range(len(station)):
                if str(station[y]) == str(sta[n]):
                    dist = distance(seismo_latitudes[y], seismo_longitudes[y], flight_latitudes[line], flight_longitudes[line])	

                    p = "/scratch/naalexeev/NODAL/2019-0"+str(month[n])+"-"+str(day[n])+"T"+str(h)+":00:00.000000Z.2019-0"+str(month[n])+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(station[y])+".mseed"
                    tr = obspy.read(p)
                    tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
                    data = tr[2][:]
                    fs = int(tr[2].stats.sampling_rate)
                    title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
                    torg = tr[2].times()
                        
                    clat, clon, dist_m, tmid = closest_encounter(flight_latitudes, flight_longitudes,line, tm, seismo_latitudes[y], seismo_longitudes[y])
                    tarrive = tim + (time[n] - calc_time(tmid,dist_m,alt_m))

                    # Compute spectrogram
                    frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9, detrend = 'constant') 

                    a, b = Sxx.shape

                    MDF = np.zeros((a,b))
                    for row in range(len(Sxx)):
                        m = len(Sxx[row])
                        p = sorted(Sxx[row])
                        median = p[int(m/2)]
                        for col in range(m):
                            MDF[row][col] = median
                    spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))
                    
                    if isinstance(sta[n], int):
                        spec = np.zeros((a,b))
                        for col in range(0,b):
                            p = sorted(Sxx[:, col])
                            median = p[int(len(p)/2)]

                            for row in range(len(Sxx)):
                                spec[row][col] = 10 * np.log10(Sxx[row][col]) - ((10 * np.log10(MDF[row][col])) + ((10*np.log10(median))))
                    
                    middle_index = len(times) // 2
                    middle_column = spec[:, middle_index]
                    vmin = 0  
                    vmax = np.max(middle_column) 
                    p, _ = find_peaks(middle_column, distance=10)

                    if n == 0:
                        tprime0 = 112
                        fnot = [93, 115, 153, 172, 228]
                        tpr = np.arange(0, 241, 1)
                        c = 343
                        v0 = 68
                        l = 2135

                    if n == 1:
                        fnot = [37, 56, 73, 110, 146, 165, 182, 218, 238, 256, 275]
                        tprime0 = 107
                        tpr = np.arange(0, 241, 1)
                        c = 343
                        v0 = 100
                        l = 2700

                    if n == 2:
                        fnot = [79,131,261]
                        tprime0 = 93
                        tpr = np.arange(0, 241, 1)
                        c = 343
                        v0 = 139
                        l = 4650

                    if n == 3:
                        fnot = [36,73,121,136,144]
                        tprime0 = 116
                        tpr = np.arange(80, 170, 1)
                        tpr = np.arange(0, 241, 1)
                        c = 343
                        v0 = 142
                        l = 2450

                    if n == 4:
                        fnot = [13,27,40,47,54,60,67,74,80,87,90,94,101,108,114,121,127,134,148,160,177,189,202,223,239,247,270]
                        tprime0 = 140
                        tpr = np.arange(40, 230, 1)
                        tpr = np.arange(0, 241, 1)
                        c = 343
                        v0 = 67
                        l = 580
            points = []
            plt.figure()
            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
            for t in range(0, len(times)):
                col = spec[50:250, t]
                maxx = np.max(spec)
                p = sorted(col)
                median = p[int(len(p)/2)]
                peaks, _ = find_peaks(col, prominence=15, distance=15, height=median + maxx/5)

                for p in peaks:
                    plt.scatter(times[t], frequencies[p]+50, color='black', marker='x')
                    points.append((times[t], frequencies[p]+50))
            plt.show()