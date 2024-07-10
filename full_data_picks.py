import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import datetime
from prelude import *
from scipy.signal import find_peaks, spectrogram
import scipy.linalg as la


seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station = seismo_data['Station']
flight_num = [530342801,528485724,528473220,528407493,528293430,531605202,531715679,529805251] 
time = [1551066051,1550172833,1550168070,1550165577,1550089044,1551662362,1551736354,1550803701] 
sta = [1022,1272,1173,1283,1004,1010,1021,1006]
day = [25,14,14,14,13,4,4,22]
month = [2,2,2,2,2,3,3,2]


for n in range(0,8):
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

                    c = 343    
                    if n == 0:
                        f0_array = [38, 57, 76, 96, 116, 135, 154, 173, 231] 
                        tprime0 = 112
                        v0 = 63
                        l = 1645
            
                    if n == 1:
                        f0_array = [36, 55, 73, 109, 146, 164, 183, 218, 236, 254, 273]
                        tprime0 = 106
                        v0 = 106
                        l = 3176

                    if n == 2:
                        f0_array = [78,119,130, 258]
                        tprime0 = 93
                        v0 = 142
                        l = 4992

                    if n == 3:
                        f0_array = [35,70,103,119,133,139]
                        tprime0 = 115
                        v0 = 160
                        l = 3832

                    if n == 4:
                        f0_array = [14,27,40,53,67,80,94,108,122,135,147,161,174,187,202,225,240,248,270]
                        tprime0 = 140
                        v0 = 62
                        l = 504

                    if n == 5:
                        f0_array = [38, 57, 76, 96, 116, 135, 154, 173, 192, 211, 231]
                        tprime0 = 112
                        v0 = 53 
                        l = 831

                    if n == 6:
                        f0_array = [19,40,59,79,100,120,140,160,180,200,221,241,261]
                        tprime0 = 118
                        v0 = 59
                        l = 479

                    if n == 7:
                        f0_array = [14,32,43,48,64,80,86,96,112,129,145,158,161,180,194,202,210,227,243,260,277]
                        tprime0 = 110
                        v0 = 89
                        l = 1307
                        
                    c = 343
                    
                    corridor_width = 6 
                    if n == 3: #if it is a Boeing Jet
                        corridor_width = 3
                    elif n == 4: # if it is a helicopter
                        corridor_width = 4
                    elif n == 7: # if it is a C46: CURTISS COMMANDO
                        corridor_width = 3 
                    middle_index =  len(times) // 2
                    middle_column = spec[:, middle_index]
                    vmin = 0  
                    vmax = np.max(middle_column) 


                    peaks_assos = []
                    fobs = []
                    tobs = []

                    for i in range(len(f0_array)):
                        f0 = f0_array[i]

                        ft = calc_ft(times,  tprime0, f0, v0, l, c)
                       
                        maxfreq = []
                        coord_inv = []
                        ttt = []

                        f01 = f0 + corridor_width
                        f02 = f0  - corridor_width
                        upper = calc_ft(times,  tprime0, f01, v0, l, c)
                        lower = calc_ft(times,  tprime0, f02, v0, l, c)

                        for t_f in range(len(times)):

                            try:      
                                tt = spec[int(np.round(lower[t_f],0)):int(np.round(upper[t_f],0)), t_f]
                                try:
                                    if n != 2 or n != 3 or n != 4 or n != 7:
                                        max_amplitude_index,_ = find_peaks(tt, prominence = 15, wlen=10, height=vmax*0.1)
                                    elif n == 7:
                                        max_amplitude_index,_ = find_peaks(tt, prominence = 1, wlen=25, height=vmax*0.2)
                                    else:
                                        max_amplitude_index,_ = find_peaks(tt, prominence = 25, wlen=5, height=vmax*0.5)
                                    maxa = np.argmax(tt[max_amplitude_index])
                                    max_amplitude_frequency = frequencies[int(max_amplitude_index[maxa])+int(np.round(lower[t_f],0))]
                                except:
                                    if n >= 3:
                                        if np.max(tt) > vmax*0.4: #This is used for the boeing jet
                                            max_amplitude_index = np.argmax(tt)
                                            max_amplitude_frequency = max_amplitude_index+int(np.round(lower[t_f],0))
                                        else:
                                            continue
                                    
                                    else:
                                        continue
                                maxfreq.append(max_amplitude_frequency)
                                coord_inv.append((times[t_f], max_amplitude_frequency))
                                ttt.append(times[t_f])

                            except:
                                continue
                        
                        if f0 < 200:
                            coord_inv_array = np.array(coord_inv)
                            mtest = [f0,v0, l, tprime0]
                            mtest,_ = invert_f(mtest, coord_inv_array, num_iterations=4)
                            ft = calc_ft(ttt,  mtest[3], mtest[0], mtest[1], mtest[2], c)
                        else:
                            ft = calc_ft(ttt,  tprime0, f0, v0, l, c)

                        delf = np.array(ft) - np.array(maxfreq)

                        count = 0
                        for i in range(len(delf)):
                            if np.abs(delf[i]) <= (3):
                                fobs.append(maxfreq[i])
                                tobs.append(ttt[i])
                                count += 1

                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                    plt.scatter(tobs,fobs, color='black', marker='x')
                    plt.show()
                    