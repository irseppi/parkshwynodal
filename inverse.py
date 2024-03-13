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
from prelude import make_base_dir, distance, closest_encounter

def df(f0,v0,l,tprime0,tprime):
    print(f0,v0,l,tprime0,tprime)
    c = 343 # m/sec speed of sound
    #t = np.sqrt(l**2+(v0*tprime0)**2/c)
   
    #t = (tprime - np.sqrt(tprime**2-(1-v0**2/c**2)*(tprime**2-l**2/c**2)))/(1-v0**2/c**2)
    #f = f0*1/(1+(v0/c)*(v0*t/(np.sqrt(l**2+(v0*t)**2))))
   
    ft0p = f0*1/(1+(v0/c)*(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
    #derivative with respect to f0
    #f_derivef0 = 1/(1+(v0/c)*(v0*t/(np.sqrt(l**2+(v0*t)**2))))
    f_derivef0 = np.gradient(ft0p, f0)

    #derivative of f with respect to v0
    f_derivev0 = np.gradient(ft0p, v0) 
    

    #derivative of f with respect to l
    f_derivel = np.gradient(ft0p, l)

    #derivative of f with respect to tprime
    #f_derivetprime0 = -(c**3 * f0 * l**2 * v0**2 * (c * np.sqrt((-(v0**2 * tprime**2) * c**2) + tprime**2 - (l**2 / c**2))) + (v0**2 - c**2) * tprime) / (np.sqrt((-(v0**2 * tprime**2) * c**2) + tprime**2 - (l**2 / c**2)) * np.sqrt(l**2 - v0 * (-np.sqrt(((-(v0**2 * tprime**2) / c**2) + tprime**2 - (l**2 / c**2))) + tprime - (v0**2 / c**2)))**2 * (c**3 * np.sqrt(l**2 - v0 * (-np.sqrt(((-(v0**2 * tprime**2) / c**2) + tprime**2 - (l**2 / c**2))) + tprime - (v0**2 / c**2))) - c**2 * v0**2 * np.sqrt((-(v0**2 * tprime**2) / c**2) + tprime**2 - (l**2 / c**2)) + c**2 * v0**2 * tprime - v0**4)**2)
    f_derivetprime0 = np.gradient(ft0p, tprime0)

    return f_derivef0, f_derivev0, f_derivel, f_derivetprime0





def invert_f(m0, coords_array, num_iterations):
    coords_array = np.array(coords_array)
    w, z = coords_array.shape
    fobs = coords_array[:,1]
    tobs = coords_array[:,0]
    n = num_iterations
    G = np.zeros((w,4)) #partial derivative matrix of f with respect to m
    m = m0

    for t in range(n):
        fnew = []
        f_derivef0, f_derivev0, f_derivel, f_derivetprime = df(m[0], m[1], m[2], m[3], tobs)
        #partial derivative matrix of f with respect to m when m=m0
        for i in range(len(coords_array)):
            
            G[i,0:4] = [f_derivef0[i], f_derivev0[i], f_derivel[i], f_derivetprime[i]]
            fnew.append(m[0]*1/(1+(m[1]/c)*(m[1]*int(tobs[i])/(np.sqrt(m[2]**2+(m[1]*tobs[i]))**2)))) # Convert m[3][i] to integer

        print((G))
        print(G.T.shape)
        m = np.reshape(np.array(m0), (4, 1)) + np.reshape(inv(G.T@G)@G.T@[np.reshape(fobs, ((i+1), 1)) - np.reshape(np.array(fnew), ((i+1), 1))], (4,1)) 
        m0 = m
    return m





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
                        fnot = [93, 115, 153, 172, 228]
                        tpr = np.arange(0, 241, 1)
                        c = 343
                        v0 = 68
                    if n == 1:
                        fnot = [71, 110, 147, 164, 182, 217, 240]
                        tprime0 = 107
                        tpr = np.arange(0, 241, 1)
                        c = 343
                        v0 = 100
                        l = 2700

                    if n == 2:
                        fnot = [131]
                        tprime0 = 93
                        tpr = np.arange(0, 241, 1)
                        c = 343
                        v0 = 139
                    if n == 3:
                        fnot = [36,73,121,136,144]
                        tprime0 = 116
                        tpr = np.arange(80, 170, 1)
                        c = 343
                        v0 = 142
                    if n == 4:
                        fnot = [13,27,40,54,67,79,93,108,120,136,147,159,175,189,202,223,239,247,270]
                        tprime0 = 140
                        tpr = np.arange(40, 230, 1)
                        c = 343
                        v0 = 64
                        l = 580

                    for f0 in fnot:
                        ft = []
                        for tprime in tpr:
                            ft0p = f0*1/(1+(v0/c)*(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
                            #pla distan 
                            ft.append(ft0p)
                        #ax2.plot(tpr, ft, 'g', linewidth=0.5)
                    #t0 = 120
                    #f0 = 120
                    #t = coords_array[:,0]
                    #c = 343
                    #v0 = speed*0.514444
                    #l = closest_encounter(flight_latitudes, flight_longitudes,line, tm, seismo_latitudes[y], seismo_longitudes[y])
                    #print(l)
                    m0 = [f0, v0, l, tprime0]

                    m = invert_f(m0, coords_array, 8)

                    f = m[0]*1/(1+(m[1]/c)*(m[1]*m[3]/(np.sqrt(m[2]**2+(m[1]*m[3])**2))))

                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                    plt.plot(f, m[3])
                    plt.show()

                    make_base_dir('/scratch/irseppi/nodal_data/plane_info/5inv_spec/')
                    plt.save('/scratch/irseppi/nodal_data/plane_info/5inv_spec/2019-02-'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'.png')
                    plt.close()