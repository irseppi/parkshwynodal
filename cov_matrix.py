import numpy.linalg as la
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import obspy
import datetime
from prelude import invert_f, distance, closest_encounter, calc_time, calc_ft,df
from scipy.signal import find_peaks, spectrogram
from matplotlib.colors import LogNorm
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station = seismo_data['Station']
flight_num = [530342801,528485724,528473220,528407493,528293430,527937367,529741194,529776675,529179112,530165646,531605202,531715679,529805251,529948401,530122923]
time = [1551066051,1550172833,1550168070,1550165577,1550089044,1549912188,1550773710,1550787637,1550511447,1550974151,1551662362,1551736354,1550803701,1550867033,1550950429]
sta = [1022,1272,1173,1283,1004,"CCB","F6TP","F4TN","F3TN","F7TV",1010,1021,1006,1109,1298]
day = [25,14,14,14,13,11,21,21,18,24,4,4,22,22,23]
month = [2,2,2,2,2,2,2,2,2,2,3,3,2,2,2]

for n in range(0,4):
    ht = datetime.datetime.utcfromtimestamp(time[n])
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
                    if isinstance(sta[n], str):
                        day_of_year = str((ht - datetime.datetime(2019, 1, 1)).days + 1)
	
                        p = "/aec/wf/2019/0"+day_of_year+"/"+str(sta[n])+".*Z.20190"+day_of_year+"000000+"
                        tr = obspy.read(p)
						
                        tr[0].trim(tr[0].stats.starttime +(int(h) *60 *60) + (mins * 60) + secs - tim, tr[0].stats.starttime +(int(h) *60 *60) + (mins * 60) + secs + tim)
                        data = tr[0][:]
                        fs = int(tr[0].stats.sampling_rate)
                        title    = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'						
                        torg                  = tr[0].times()
                    else:
                        p = "/scratch/naalexeev/NODAL/2019-0"+str(month[n])+"-"+str(day[n])+"T"+str(h)+":00:00.000000Z.2019-0"+str(month[n])+"-"+str(day2)+"T"+str(h_u)+":00:00.000000Z."+str(station[y])+".mseed"
                        tr = obspy.read(p)
                        tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
                        data = tr[2][:]
                        fs = int(tr[2].stats.sampling_rate)
                        title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
                        torg = tr[2].times()
                      
                    dist_m, tmid = closest_encounter(flight_latitudes, flight_longitudes,line, tm, seismo_latitudes[y], seismo_longitudes[y])
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
                    ty = False
                    if ty == True:
                        if isinstance(sta[n], str):
                            spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))
                        else:
                            spec = np.zeros((a,b))
                            for col in range(0,b):
                                p = sorted(Sxx[:, col])
                                median = p[int(len(p)/2)]

                                for row in range(len(Sxx)):

                                    spec[row][col] = 10 * np.log10(Sxx[row][col]) - ((10 * np.log10(MDF[row][col])) + ((10*np.log10(median))))
                    else:
                        spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))		
                    middle_index = len(times) // 2
                    middle_column = spec[:, middle_index]
                    vmin = 0  
                    vmax = np.max(middle_column) 
                    p, _ = find_peaks(middle_column, distance=10)


                    coords = []
                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)

                    def onclick(event):
                        global coords
                        coords.append((event.xdata, event.ydata))
                        plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                        plt.draw() 
                        print('Clicked:', event.xdata, event.ydata)  

                    cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

                    plt.show(block=True)
                    # Convert the list of coordinates to a numpy array
                    coords_array = np.array(coords)
                
                    if n == 0:
                        tprime0 = 112
                        f0 = 115
                        v0 = 68
                        l = 2135

                    if n == 1:
                        f0 = 110
                        tprime0 = 107
                        v0 = 100
                        l = 2700

                    if n == 2:
                        f0 = 131
                        tprime0 = 93
                        v0 = 139
                        l = 4650

                    if n == 3:
                        f0 = 121
                        tprime0 = 116
                        v0 = 142
                        l = 2450

                    if n == 4:
                        f0 = 120
                        tprime0 = 140
                        v0 = 64
                        l = 580

                    if n == 5:
                        f0 = 17.5
                        tprime0 = 123
                        v0 = 112
                        l = 1150

                    if n == 6:
                        f0 = 36
                        tprime0 = 133
                        v0 = 92
                        l = 2400

                    if n == 7:
                        f0 = 26		
                        tprime0 = 122
                        v0 = 126
                        l = 3000

                    if n == 8:
                        f0 = 87.7
                        tprime0 = 100
                        v0 = 67
                        l = 2300

                    if n == 9:
                        f0 = 26
                        tprime0 = 114
                        v0 = 144
                        l = 1900
                    if n > 9:
                        f0 = fs/2
                        tprime0 = tarrive
                        v0 = speed_mps
                        l = np.sqrt(dist_m**2 + alt_m**2)

                    c = 343
                    m0 = [f0, v0, l, tprime0]

                    m = invert_f(m0, coords_array, num_iterations=8)
                    f0 = m[0]
                    v0 = m[1]
                    l = m[2]
                    tprime0 = m[3]
                    
                    ft = calc_ft(times, tprime0, f0, v0, l, c)
                    if isinstance(sta[n], int):
                        peaks = []
                        p, _ = find_peaks(middle_column, distance=7)
                        corridor_width = fs / len(p)                 
                        if len(p) == 0:
                            corridor_width = fs


                        coord_inv = []
                        for t_f in range(len(times)):
                            upper = int(ft[t_f] + corridor_width)
                            lower = int(ft[t_f] - corridor_width)
                            if lower < 0:
                                lower = 0
                            if upper > len(frequencies):
                                upper = len(frequencies)
                            tt = spec[lower:upper, t_f]

                            max_amplitude_index = np.argmax(tt)
                            
                            max_amplitude_frequency = frequencies[max_amplitude_index+lower]
                            peaks.append(max_amplitude_frequency)
                            coord_inv.append((times[t_f], max_amplitude_frequency))



                        coord_inv_array = np.array(coord_inv)

                        m = invert_f(m0, coord_inv_array, num_iterations=12)
                        f0 = m[0]
                        v0 = m[1]
                        l = m[2]
                        tprime0 = m[3]

                        ft = calc_ft(times, tprime0, f0, v0, l, c)
                        
                        delf = np.array(ft) - np.array(peaks)
                        
                        new_coord_inv_array = []
                        for i in range(len(delf)):
                            if np.abs(delf[i]) <= 3:
                                new_coord_inv_array.append(coord_inv_array[i])
                        coord_inv_array = np.array(new_coord_inv_array)


                        m = invert_f(m0, coord_inv_array, num_iterations=12)
                        f0 = m[0]
                        v0 = m[1]
                        l = m[2]
                        tprime0 = m[3]

                        ft = calc_ft(times, tprime0, f0, v0, l, c)

                        
                        
                        w,_ = coord_inv_array.shape
                        fobs = coord_inv_array[:,1]
                        tobs = coord_inv_array[:,0]
                        c = 343
                        G = np.zeros((w,4)) #partial derivative matrix of f with respect to m
                        #partial derivative matrix of f with respect to m 
                        for i in range(0,w):
                            f0 = m[0]
                            v0 = m[1]
                            l = m[2]
                            tprime0 = m[3]
                            tprime = tobs[i]
                            t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
                            ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
                            f_derivef0, f_derivev0, f_derivel, f_derivetprime0 = df(m[0], m[1], m[2], m[3], tobs[i])
                            
                            G[i,0:4] = [f_derivef0, f_derivev0, f_derivel, f_derivetprime0]

                        def plot_matrix(M,gridlines=False,colormap=None):
                            
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

                        #covm = np.cov(G.T@G)
                        sigma = 20
                        covm = (sigma**2)*la.inv(G.T@G)
                        plt.figure()
                        plot_matrix(np.log(np.abs(covm)),True,'binary')
                        plt.xlabel('i')
                        plt.ylabel('i')
                        plt.colorbar(norm=LogNorm())
                        #plt.clim(0, 2000)  # Set the colorbar limits
                        plt.show()
