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
elevations = seismo_data['Elevation']

sta_f = open('input/all_station_crossing_db.txt','r')

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3
f = 1.0/np.cos(62*np.pi/180)

tim = 120

# Loop through each station in text file that we already know comes within 2km of the nodes
for line in sta_f.readlines():
    val = line.split(',')
    date = val[0]
    flight = val[1]
    sta = val[5]
    equipment = val[6][0:4]

    tm = float(val[2])

    ht = datetime.utcfromtimestamp(tm)
    mins = ht.minute
    secs = ht.second
    month = ht.month
    day = ht.day

    h = ht.hour
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

    flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/20190'+str(month)+str(day)+'_flights.csv', sep=",")
    flight_id = flight_data['flight_id']
    callsign = flight_data['callsign'] 
    fly = flight_data['flight']

    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    time = flight_data['snapshot_id']
    head = flight_data['heading']
    speed = flight_data['speed']
    altitude = flight_data['altitude']
    try:
        p = "/scratch/naalexeev/NODAL/2019-0"+str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month)+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(sta)+".mseed"
        tr = obspy.read(p)
    except:
        continue

    tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
    data = tr[2][:]
    fs = int(tr[2].stats.sampling_rate)
    title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
    torg = tr[2].times()

    for n in range(len(time)):         
        if str(tm) == str(time[n])+'.0':
            spd = speed[n]
            speed_mps = spd * 0.514444
            alt = altitude[n]
            alt_m = alt * 0.3048
            index = n

            for e in range(len(station)):
                if station[e] == sta:
                    elevation = elevations[e]
                    clat, clon, dist_m, tmid = closest_encounter(flight_latitudes, flight_longitudes, index, time, seismo_latitudes[e], seismo_longitudes[e])
                else:
                    continue
        else:
            continue
    height_m = alt_m - elevation 
    tarrive = 120+ (time[n] - calc_time(tmid,dist_m,height_m))

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

    middle_index =  len(times) // 2
    middle_column = spec[:, middle_index]
    vmin = 0  
    vmax = np.max(middle_column) 

    c = 343    
    
    tprime0 = tarrive
    v0 = speed_mps
    l = np.sqrt(dist_m**2 + (alt_m-elevation)**2)

    if equipment == 'C185':
        f0_array = [38, 57, 76, 96, 116, 135, 154, 173, 192, 211, 231]


    elif equipment == 'PA31':
        f0_array = [36, 55, 73, 109, 146, 164, 183, 218, 236, 254, 273]


    elif equipment == 'SW4':
        f0_array = [78,119,130, 258]


    elif equipment == 'B736':
        f0_array = [35,70,103,119,133,139]


    if equipment == 'R44':
        f0_array = [14,27,40,53,67,80,94,108,122,135,147,161,174,187,202,225,240,248,270]

    elif equipment == 'GA8':
        f0_array = [19,40,59,79,100,120,140,160,180,200,221,241,261]


    elif equipment == 'C46':
        f0_array = [14,32,43,48,64,80,86,96,112,129,145,158,161,180,194,202,210,227,243,260,277]

    else:
        continue

    plt.figure()
    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax) 
    plt.show()

    con = input("Do you want to use this? (y or n)")
    if con == 'n':
        continue
    else:
        corridor_width = 6 
        if equipment == 'B736': #if it is a Boeing Jet
            corridor_width = 3
        elif equipment == 'R44': # if it is a helicopter
            corridor_width = 4
        elif equipment == 'C46': # if it is a C46: CURTISS COMMANDO
            corridor_width = 3 

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
                        if equipment != 'SW4' or equipment != 'B736' or equipment != 'R44' or equipment != 'C46':
                            max_amplitude_index,_ = find_peaks(tt, prominence = 15, wlen=10, height=vmax*0.1)
                        elif equipment == 'C46':
                            max_amplitude_index,_ = find_peaks(tt, prominence = 1, wlen=25, height=vmax*0.2)
                        else:
                            max_amplitude_index,_ = find_peaks(tt, prominence = 25, wlen=5, height=vmax*0.5)
                        maxa = np.argmax(tt[max_amplitude_index])
                        max_amplitude_frequency = frequencies[int(max_amplitude_index[maxa])+int(np.round(lower[t_f],0))]
                    except:
                        if equipment == 'B736' or len(f0_array) > 11: #This is used for the boeing jet and any other flight with more than 11 fundamental frequencies
                            if np.max(tt) > vmax*0.4: 
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
            if len(coord_inv) > 1:
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
            else:
                continue
        if len(fobs) == 0:
            print('No picks found')
            continue
        time_pick = False
        if time_pick == True:
            set_time = []
            plt.figure()
            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
            plt.scatter(tobs,fobs, color='black', marker='x')
            def onclick(event):
                global coords
                set_time.append(event.xdata) 
                plt.scatter(event.xdata, event.ydata, color='red', marker='x')  # Add this line
                plt.draw() 
                print('Clicked:', event.xdata, event.ydata)  

            cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
            plt.show(block=True)

            start_time = set_time[0]
            end_time = set_time[1]

            ftobs = []
            ffobs = []
            for i in range(len(tobs)):
                if tobs[i] >= start_time and tobs[i] <= end_time:
                    ftobs.append(tobs[i])
                    ffobs.append(fobs[i])
            tobs = ftobs
            fobs = ffobs

        plt.figure()
        plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
        plt.scatter(tobs,fobs, color='black', marker='x')
        plt.show()
