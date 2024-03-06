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

def df(f0,v0,l,t):
    c = 343 # m/sec speed of sound
    tprime = np.sqrt(l**2-(v0*t)**2/c)
    t = (tprime - np.sqrt(tprime**2-(1-v0**2/c**2)(tprime**2-l**2/c**2)))/(1-v0**2/c**2)

    f = f0*1/(1+(v0/c)*(v0*t/(np.sqrt(l**2-(v0*t)**2))))


    #derivative with respect to f0
    f_derivef0 = 1/(1+(v0/c)*(v0*t/(np.sqrt(l**2-(v0*t)**2))))
    #derivative of f with respect to v0
    f_derivev0 = c*f0*t*v0*(t**2*v0**2-l**2)*1/(np.sqrt(l**2-(v0*t)**2)*(c*np.sqrt(l**2-(v0*t)**2)+t*v0**2)**2)
    #derivative of f with respect to l
    f_derivel = (c*f0*t*v0**2*l)*1/(np.sqrt(l**2-(v0*t)**2)*(c*np.sqrt(l**2-(v0*t)**2)+t*v0**2)**2)
    #derivative of f with respect to tprime
    f_derivetprime = (c**3*f0*l**2*v0**2*(c*np.sqrt((-(v0**2*tprime**2)*c**2)+tprime**2-(l**2/c**2)))+(v0**2-c**2)*tprime)/(np.sqrt((-(v0**2*tprime**2)*c**2)+tprime**2-(l**2/c**2))*np.sqrt(l**2-v0*(-np.sqrt(((-(v0**2*tprime**2)/c**2)+tprime**2-(l**2/c**2)))+tprime-(v0**2/c**2))**2)*(c**3*np.sqrt(l**2-v0*(-np.sqrt(((-(v0**2*tprime**2)/c**2)+tprime**2-(l**2/c**2)))+tprime-(v0**2/c**2)))-c**2*v0**2*np.sqrt((-(v0**2*tprime**2)/c**2)+tprime**2-(l**2/c**2))+c**2*v0**2*tprime-v0**4)**2)
    
    return f_derivef0, f_derivev0, f_derivel, f_derivetprime, tprime


def invert_f(m0, coords_array, num_iterations):
    coords_array = np.array(coords_array)
    w, z = coords_array.shape
    fobs = coords_array[:,1]
    n = num_iterations
    G = np.zeros((w,4)) #partial derivative matrix of f with respect to m
    m = m0
    for t in range(n):
        tpa = []
        for i in range(len(coords_array)):
            f_derivef0, f_derivev0, f_derivel, f_derivetprime, tprime = df(m[0], m[1], m[2], m[3])
            tpa.append(tprime)
            #partial derivative matrix of f with respect to m when m=m0
            G[i,:] = [[f_derivef0, f_derivev0, f_derivel, f_derivetprime]] 

            #update m
            f = m[0]*1/(1+(m[1]/c)*(m[1]*m[3]/(np.sqrt(m[2]**2-(m[1]*m[3])**2))))
            m = m0 +inv(G.T@G)@G.T[fobs - f]
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


                f0 = 120
                t = coords_array[:,0]
                c = 343
                v0 = speed*0.514444
                l = closest_encounter(flight_latitudes, flight_longitudes,line, tm, seismo_latitudes[y], seismo_longitudes[y])
                m0 = [f0, v0, l, t]
                m = invert_f(m0, coords_array, 8)

                f = m[0]*1/(1+(m[1]/c)*(m[1]*m[3]/(np.sqrt(m[2]**2-(m[1]*m[3])**2))))

                plt.figure()
                plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                plt.plot(f, m[3])
                plt.show()

                make_base_dir('/scratch/irseppi/nodal_data/plane_info/5inv_spec/')
                plt.save('/scratch/irseppi/nodal_data/plane_info/5inv_spec/2019-02-'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'.png')
                plt.close()