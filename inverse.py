import numpy as np
from numpy.linalg import inv
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy import signal
import obspy
import math
from obspy.core import UTCDateTime
import datetime
from prelude import make_base_dir, distance
import numpy as np
from scipy.signal import find_peaks

def df(f0,v0,l,tp0,tp):   
    c = 343 # m/sec speed of sound

    #derivative with respect to f0
    f_derivef0 = (1 / (1 - (c * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))) /((c**2 - v0**2) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2))))


    #derivative of f with respect to v0
    f_derivev0 = (-f0 * v0 * (-2 * l**4 * v0**4 + l**2 * (tp - tp0)**2 * v0**6 + c**6 * (tp - tp0) * (2 * l**2 + (tp - tp0)**2 * v0**2) * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) + 
    c**2 * (4 * l**4 * v0**2 - (tp - tp0)**4 * v0**6 + l**2 * (tp - tp0) * v0**4 * (5 * tp - 5 * tp0 - 3 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))) - c**4 * 
    (2 * l**4 - 3 * (tp - tp0)**3 * v0**4 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4)) - l**2 * (tp - tp0) * v0**2 * (-6 * tp + 6 * tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
    (l**2 + (tp - tp0)**2 * v0**2))/c**4)))) / (c * (c - v0) * (c + v0) * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
    (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) - c**2 * np.sqrt(l**2 + (c**4 * v0**2 * 
    (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2))**2))
    


    #derivative of f with respect to l
    f_derivel = ((f0 * l * (tp - tp0) * (c - v0) * v0**2 * (c + v0) * ((-tp + tp0) * v0**2 + c**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))) / 
    (c * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + 
    (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4) - 
    c**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + 
    (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2))**2))


    #derivative of f with respect to tprime0
    f_derivetprime0 = ((f0 * l**2 * (c - v0) * v0**2 * (c + v0) * ((-tp + tp0) * v0**2 + c**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))) / 
    (c * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
    (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) - 
    c**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + 
    (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2))**2))


    return f_derivef0, f_derivev0, f_derivel, f_derivetprime0


def invert_f(m0, coords_array, num_iterations):
    w,_ = coords_array.shape
    fobs = coords_array[:,1]
    tobs = coords_array[:,0]
    m = m0
    n = 0
    c = 343
    while n < num_iterations:
        fnew = []
        G = np.zeros((w,4)) #partial derivative matrix of f with respect to m
        #partial derivative matrix of f with respect to m 
        for i in range(0,w):
            f0 = m[0]
            v0 = m[1]
            l = m[2]
            tprime0 = m[3]
            tprime = tobs[i]
            ft0p = f0*1/(1+(v0/c)*(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
            f_derivef0, f_derivev0, f_derivel, f_derivetprime0 = df(m[0], m[1], m[2], m[3], tobs[i])
            
            G[i,0:4] = [f_derivef0, f_derivev0, f_derivel, f_derivetprime0]

            fnew.append(ft0p) 
    
        m = np.reshape(np.reshape(m0,(4,1))+ np.reshape(inv(G.T@G)@G.T@(np.reshape(fobs, (len(coords_array), 1)) - np.reshape(np.array(fnew), (len(coords_array), 1))), (4,1)), (4,))
        print(m)
        m0 = m
        n += 1
    return m

num_iterations = 8
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station = seismo_data['Station']
flight_num = [530342801,528485724,528473220,528407493,528293430]
time = [1551066051,1550172833,1550168070,1550165577,1550089044]
sta = [1022,1272,1173,1283,1004]
day = [25,14,14,14,13]

for n in range(0,5):
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
    flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/201902'+str(day[n])+'_positions/201902'+str(day[n])+'_'+str(flight_num[n])+'.csv', sep=",")

    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    tm = flight_data['snapshot_id']
    speed = flight_data['speed']
    alt = flight_data['altitude']
    head = flight_data['heading']
    for line in range(len(tm)):
        if str(tm[line]) == str(time[n]):
            speed = flight_data['speed'][line]
            alt = flight_data['altitude'][line]
            for y in range(len(station)):
                if str(station[y]) == str(sta[n]):
                    dist = distance(seismo_latitudes[y], seismo_longitudes[y], flight_latitudes[line], flight_longitudes[line])	

                    p = "/scratch/naalexeev/NODAL/2019-02-"+str(day[n])+"T"+str(h)+":00:00.000000Z.2019-02-"+day2+"T"+h_u+":00:00.000000Z."+station[y]+".mseed"
                    tr = obspy.read(p)
                    tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
                    data = tr[2][0:-1]
                    fs = int(tr[2].stats.sampling_rate)
                    title = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
                    t = tr[2].times()
                    # Time array
                    t = np.arange(len(data)) / fs
                    g = fs*240
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
                    middle_index = len(times) // 2
                    middle_column = spec[:, middle_index]
                    vmin = 0  
                    vmax = np.max(middle_column) 
                    p, _ = signal.find_peaks(middle_column, distance=10)
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
                   
                    c = 343
                    m0 = [f0, v0, l, tprime0]

                    m = invert_f(m0, coords_array, num_iterations=4)
                    ft = []
                    t = np.arange(0, 241, 1)
                    for tprime in times:
                        f0 = m[0]
                        v0 = m[1]
                        l = m[2]
                        tprime0 = m[3]

                        ft0p = f0*1/(1+(v0/c)*(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
                            
                        ft.append(ft0p)
                    
                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                    plt.plot(times, ft, 'g', linewidth=0.5)
                    plt.show()

                    peaks = []
                    p, _ = signal.find_peaks(middle_column, distance=7)
                    corridor_width = 250 / len(p)

                    coord_inv = []
                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)

                    for t_f in range(len(times)):
                        upper = int(ft[t_f] + corridor_width)
                        lower = int(ft[t_f] - corridor_width)

                        plt.scatter(times[t_f], upper, color='pink', marker='x')
                        plt.scatter(times[t_f], lower, color='pink', marker='x')

                        tt = spec[lower:upper, t_f]

                        max_amplitude_index = np.argmax(tt)
                        
                        max_amplitude_frequency = frequencies[max_amplitude_index+lower]
                        peaks.append(max_amplitude_frequency)
                        coord_inv.append((times[t_f], max_amplitude_frequency))
                        plt.scatter(times[t_f], max_amplitude_frequency, color='black', marker='x')
                       
                    plt.show()
                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                    plt.scatter(times, peaks, color='black', marker='x')
                    plt.show()
                    coord_inv_array = np.array(coord_inv)
                    m = invert_f(m0, coord_inv_array, num_iterations=8)
                    ft = []
                    t = np.arange(0, 241, 1)
                    for tprime in times:
                        f0 = m[0]
                        v0 = m[1]
                        l = m[2]
                        tprime0 = m[3]

                        ft0p = f0*1/(1+(v0/c)*(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
                            
                        ft.append(ft0p)
                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                    plt.plot(times, ft, 'g', linewidth=0.5)
                    plt.show()
                    #plt.figure()
                    #plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                    #plt.plot(t, ft, 'g', linewidth=0.5)
                    #plt.show()

                    #make_base_dir('/scratch/irseppi/nodal_data/plane_info/5inv_spec/')
                    #fig.save('/scratch/irseppi/nodal_data/plane_info/5inv_spec/2019-02-'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'.png')
                    #plt.close()
           