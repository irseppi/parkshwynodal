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
flight_num = [530342801,528485724,528473220,528407493,528293430] 
time = [1551066051,1550172833,1550168070,1550165577,1550089044] 
sta = [1022,1272,1173,1283,1004]
day = [25,14,14,14,13]
month = [2,2,2,2,2]
def plot_matrix(M,gridlines=False,colormap='gray'):
    
    plt.imshow(M,cmap=colormap)
    plt.xticks(ticks=range(np.shape(M)[1]),labels=[str(val) for val in range(1,np.shape(M)[1]+1)])
    plt.yticks(ticks=range(np.shape(M)[0]),labels=[str(val) for val in range(1,np.shape(M)[0]+1)])
    if gridlines:
        xgrid = np.array(range(np.shape(M)[1] + 1)) - 0.5
        ygrid = np.array(range(np.shape(M)[0] + 1)) - 0.5
        for gridline in xgrid:
            plt.axvline(x=gridline,color='k',linewidth=1)
        for gridline in ygrid:
            plt.axhline(y=gridline,color='k',linewidth=1)
    plt.show()
def Sd(m,dobs,dpred,icobs):
    sd = 0.5 * (dpred-dobs).T @ icobs @ (dpred-dobs)
    return sd
# model misfit (related to regularization)
def Sm(m,mprior,icprior):
    sm = 0.5 * (m-mprior).T @ icprior @ (m-mprior)
    return sm
# total misfit
def S(m,dobs,dpred,mprior,icobs,icprior):
    s = Sd(m,dobs,dpred,icobs) + Sm(m,mprior,icprior)
    return s

for n in range(4,5):
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
                    '''
                    if isinstance(sta[n], int):
                            spec = np.zeros((a,b))
                            for col in range(0,b):
                                p = sorted(Sxx[:, col])
                                median = p[int(len(p)/2)]

                                for row in range(len(Sxx)):
                                    spec[row][col] = 10 * np.log10(Sxx[row][col]) - ((10 * np.log10(MDF[row][col])) + ((10*np.log10(median))))
                    '''
                    middle_index = len(times) // 2
                    middle_column = spec[:, middle_index]
                    vmin = 0  
                    vmax = np.max(middle_column) 
                    p, _ = find_peaks(middle_column, distance=10)

                    output1 = '/scratch/irseppi/nodal_data/plane_info/inversepicks/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(time[n])+'_'+str(flight_num[n])+'.csv'
                    coords = []
                    with open(output1, 'r') as file:
                        for line in file:
                            # Split the line using commas
                            pick_data = line.split(',')
                            coords.append((float(pick_data[0]), float(pick_data[1])))
                    file.close()  # Close the file after reading

                    coords_array = np.array(coords)

                    if n == 0:
                        tprime0 = 112
                        fnot = [38, 57, 76, 96, 135, 154, 173, 231]
                        v0 = 63
                        l = 1645
                        f0 = 116

                    if n == 1:
                        fnot = [36, 55, 73, 146, 164, 183, 218, 236, 254, 273]
                        tprime0 = 106
                        f0 = 109
                        v0 = 106
                        l = 3176

                    if n == 2:
                        fnot = [78, 120, 258]
                        tprime0 = 93
                        f0 = 130
                        v0 = 138
                        l = 4992

                    if n == 3:
                        fnot = [34,69,104,134,140]
                        tprime0 = 115
                        f0 = 118
                        v0 = 120
                        l = 2000

                    if n == 4:
                        fnot = [13,26,40,53,67,80,93,110,119,148,160,174,187,202,226,241,249,272]
                        tprime0 = 135
                        f0 = 135
                        v0 = 79
                        l = 580
                   
                    w  = len(fnot)
                    
                    mprior = []
                    mprior.append(v0)
                    mprior.append(l)
                    mprior.append(tprime0)
                    for i in range(w+1):
                        if i == 0:
                            mprior.append(f0)
                        else:
                            mprior.append(fnot[i-1])
                    mprior = np.array(mprior)
                    vnot = v0
                    lnot = l
                    tprimenot = tprime0
                    c = 343
                    m0 = [f0, v0, l, tprime0]
                    m,covm = invert_f(m0, coords_array, num_iterations=4)

                    p, _ = find_peaks(middle_column, distance = 7)
                    corridor_width = 3 #(fs/2) / len(p)
                    #if n == 3 or n == 4:
                    #    corridor_width = 2
                    peaks_assos = []
                    fobs = []
                    tobs = []
                    f0_array = []

                    f0_array.append(f0)
                   
                    ft = calc_ft(times, m[3], m[0], m[1], m[2], c)

                    maxfreq = []
                    coord_inv = []
                    ttt = []
                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                    
                    start_time = 0
                    end_time = 239
                    set_time = []
                    def onclick(event):
                        global coords
                        set_time.append(event.xdata) 
                        plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                        plt.draw() 
                        print('Clicked:', event.xdata, event.ydata)  

                    cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
                    plt.show(block=True)
                    
                    start_time = set_time[0]
                    end_time = set_time[1]

                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                    new_times = times[np.where((times >= start_time) & (times <= end_time))]
                    ft = calc_ft(new_times, m[3], m[0], m[1], m[2], c)
                    for t_f in range(len(new_times)):
                        if not np.isnan(ft[t_f]) and ft[t_f] != np.inf:
                            upper = int(ft[t_f] + corridor_width)
                            lower = int(ft[t_f] - corridor_width)

                            if lower < 0:
                                lower = 0
                            if upper > len(frequencies):
                                upper = len(frequencies)
                            try:
                                tt = spec[lower:upper, t_f]

                                max_amplitude_index = np.argmax(tt)
                                max_amplitude_frequency = frequencies[max_amplitude_index+lower]
                                maxfreq.append(max_amplitude_frequency)
                                coord_inv.append((new_times[t_f], max_amplitude_frequency))
                                ttt.append(new_times[t_f])
                            except:
                                continue
                            
                    coord_inv_array = np.array(coord_inv)
                    if len(coord_inv_array) == 0:
                        continue
                    m,_ = invert_f(m0, coord_inv_array, num_iterations=4)
                    f0 = m[0]
                    v00 = m[1]
                    l0 = m[2]
                    tprime00 = m[3]
                    ft = calc_ft(ttt, tprime00, f0, v00, l0, c)

                    delf = np.array(ft) - np.array(maxfreq)

                    count = 0
                    for i in range(len(delf)):
                        if np.abs(delf[i]) <= (4):
                            fobs.append(maxfreq[i])
                            tobs.append(ttt[i])
                            count += 1
                    
                    peaks_assos.append(count)
                    corridor_width = 3
                    for i in range(len(fnot)):
                        f0 = fnot[i]
                        f0_array.append(f0)
                        m0 = [f0, vnot, lnot, tprimenot]
                        ft = calc_ft(new_times,  tprime00, f0, v00, l0, c)
                        
                        maxfreq = []
                        coord_inv = []
                        ttt = []
                        for t_f in range(len(new_times)):
                            
                            if not np.isnan(ft[t_f]) and ft[t_f] != np.inf:
                                upper = int(ft[t_f] + corridor_width)
                                lower = int(ft[t_f] - corridor_width)

                                if lower < 0:
                                    lower = 0
                                if upper > len(frequencies):
                                    upper = len(frequencies)
                                try:
                                    tt = spec[lower:upper, t_f]

                                    max_amplitude_index = np.argmax(tt)
                                    max_amplitude_frequency = frequencies[max_amplitude_index+lower]
                                    maxfreq.append(max_amplitude_frequency)
                                    coord_inv.append((new_times[t_f], max_amplitude_frequency))
                                    ttt.append(new_times[t_f])
                                    #fobs.append(max_amplitude_frequency)
                                    #tobs.append(new_times[t_f])
                                except:
                                    continue
                        ft = calc_ft(ttt,  tprime00, f0, v00, l0, c)
                        delf = np.array(ft) - np.array(maxfreq)

                        count = 0
                        for i in range(len(delf)):
                            if np.abs(delf[i]) <= (2):
                                fobs.append(maxfreq[i])
                                tobs.append(ttt[i])
                                count += 1
                    plt.scatter(tobs, fobs, color='black', marker='x')

                    plt.show()
                    plt.close()
                   