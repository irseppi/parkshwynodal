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
elevations = seismo_data['Elevation']

sta_f = open('input/all_station_crossing_db_C185.txt','r')

# Loop through each station in text file that we already know comes within 2km of the nodes
for line in sta_f.readlines():
    val = line.split(',')
    date = val[0]
    flight = val[1]
    sta = val[5]
    equipment = val[6][0:4]

    tm = float(val[2])
    if datetime.utcfromtimestamp(tm).month == 3:

        flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight) + '.csv'
        flight_data = pd.read_csv(flight_file, sep=",")
        flight_latitudes = flight_data['latitude']
        flight_longitudes = flight_data['longitude']
        time = flight_data['snapshot_id']

        speed = flight_data['speed']
        altitude = flight_data['altitude']

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
        tarrive = calc_time(tmid,dist_m,height_m)

        ht = datetime.utcfromtimestamp(tarrive)
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

        try:
            p = "/scratch/naalexeev/NODAL/2019-0"+str(month)+"-"+str(day)+"T"+str(h)+":00:00.000000Z.2019-0"+str(month)+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(sta)+".mseed"
            tr = obspy.read(p)
        except:
            continue

        tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - 120, tr[2].stats.starttime + (mins * 60) + secs + 120)
        data = tr[2][:]
        fs = int(tr[2].stats.sampling_rate)
        title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
        torg = tr[2].times()

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
        try_median = False
        if try_median == True:
            if isinstance(sta, int):
                spec = np.zeros((a,b))
                for col in range(0,b):
                    p = sorted(Sxx[:, col])
                    median = p[int(len(p)/2)]

                    for row in range(len(Sxx)):
                        spec[row][col] = 10 * np.log10(Sxx[row][col]) - ((10 * np.log10(MDF[row][col])) + ((10*np.log10(median))))
        middle_index =  len(times) // 2
        middle_column = spec[:, middle_index]
        vmin = 0  
        vmax = np.max(middle_column) 

        c = 343
        tprime0 = 120
        v0 = speed_mps
        l = np.sqrt(dist_m**2 + (height_m)**2)
        
        f0_array = [38, 57, 76, 96, 116, 135, 154, 173, 192, 211, 231]
        #f0_array = [62,82,105,124,145,167,186,209,250]

        tf = np.arange(0, 240, 1)

        fig, ax1 = plt.subplots(1, 1)   

        # Plot spectrogram
        cax = ax1.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)		
		
        ax1.set_xlabel('Time (s)')

        for pp in range(len(f0_array)):
            f0 = f0_array[pp]
 
            ft = calc_ft(times, tprime0, f0, v0, l, c)

            ax1.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) 

        ax1.axvline(x=120, color='black', linestyle='--', linewidth=0.7)


        ax1.margins(x=0)
        ax2 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

        plt.colorbar(mappable=cax, cax=ax2)


        ax1.margins(x=0)
        ax1.set_xlim(0, 240)
        ax1.set_ylim(0, int(fs/2))

        plt.show()

        con = input("Do you want to use this? (y or n)")
        if con == 'n':
            continue
        else:

            corridor_width = 6 

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
                            max_amplitude_index,_ = find_peaks(tt, prominence = 15, wlen=10, height=vmax*0.1)
                            maxa = np.argmax(tt[max_amplitude_index])
                            max_amplitude_frequency = frequencies[int(max_amplitude_index[maxa])+int(np.round(lower[t_f],0))]
                        except:
                            continue
                        maxfreq.append(max_amplitude_frequency)
                        coord_inv.append((times[t_f], max_amplitude_frequency))
                        ttt.append(times[t_f])

                    except:
                        continue
                if len(coord_inv) > 0:
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

    else:
        continue