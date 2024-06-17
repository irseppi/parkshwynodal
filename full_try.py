import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import datetime
from prelude import *
from scipy.signal import find_peaks, spectrogram
from matplotlib.patches import Rectangle
from matplotlib.patches import Rectangle
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
                peaks, _ = find_peaks(col, prominence=15, distance=10, height=median + maxx/5)

                for p in peaks:
                    plt.scatter(times[t], frequencies[p]+50, color='black', marker='x')
                    points.append((times[t], frequencies[p]+50))
            x = []
            y = []
            def onclick(event):
                global coords
                x.append(event.xdata)
                y.append(event.ydata)
                plt.scatter(event.xdata, event.ydata, color='red', marker='x')  # Add this line
                plt.draw() 

            cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

            plt.show(block=True)
            points = np.array(points)

            points_filt = []
            throw_out = []
            plt.figure()
            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
            for ll in range(0,len(x),2):
                for p in range(len(points)):
                    if x[ll] < points[p][0] < x[ll+1] and y[ll+1] <= points[p][1] <= y[ll]:
                        throw_out.append(points[p])
                    else:
                        points_filt.append(points[p])
                        plt.scatter(points[p][0], points[p][1], color='black', marker='x')
            plt.show()

            '''
            qv = 0
            num_iterations = 8
            while qv < num_iterations:
                w  = len(peaks)
                G = np.zeros((0,w+3))
                fnew = []
                if qv == 0:
                    m0 = []
                    m0.append(v0)
                    m0.append(l)
                    m0.append(tprime0)

                    for i in range(w):
                        m0.append(f0_array[i])

                for p in range(w):
                    new_row = np.zeros(w+3)
                    f0 = f0_array[p]
                    if p == 0:
                        for j in range(peaks_assos[p]):
                            tprime = tobs[j]
                            t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
                            ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))

                            f_derivef0, f_derivev0, f_derivel, f_derivetprime0 = df(f0,v0,l,tprime0, tobs[j])
                        
                            new_row[0] = f_derivev0
                            new_row[1] = f_derivel
                            new_row[2] = f_derivetprime0
                            new_row[3+p] = f_derivef0

                            G = np.vstack((G, new_row))

                            fnew.append(ft0p)
                    else:
                        for j in range(peaks_assos[p-1],peaks_assos[p]+peaks_assos[p-1]):
                            tprime = tobs[j]
                            t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
                            ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))

                            f_derivef0, f_derivev0, f_derivel, f_derivetprime0 = df(f0,v0,l,tprime0, tobs[j])
                        
                            new_row[0] = f_derivev0
                            new_row[1] = f_derivel
                            new_row[2] = f_derivetprime0
                            new_row[3+p] = f_derivef0

                            G = np.vstack((G, new_row))

                            fnew.append(ft0p)
                    print(new_row)
                sigma = 5
                
                # Exclude rows with NaN or inf values
                valid_rows = np.isfinite(G).all(axis=1)
                Gv = G[valid_rows]
                
                fobs = np.array(fobs)
                fnew = np.array(fnew)
                fnewv = fnew[valid_rows]
                fobsv = fobs[valid_rows]

                pin = la.pinv(Gv.T@Gv)

                m = np.reshape(np.reshape(m0,(3+w,1))+ np.reshape(pin@Gv.T@(np.reshape(fobsv, (len(fobsv), 1)) - np.reshape(np.array(fnewv), (len(fobsv), 1))), (3+w,1)), (3+w,))
                covmlsq = (sigma**2)*pin
                v0 = m[0]
                l = m[1]
                tprime0 = m[2]
                f0_array = m[3:]

                m0 = m
                print(m)
                qv += 1

            closest_index = np.argmin(np.abs(tprime0 - times))
            arrive_time = spec[:,closest_index]
            for i in range(len(arrive_time)):
                if arrive_time[i] < 0:
                    arrive_time[i] = 0
            vmin = np.min(arrive_time) 
            vmax = np.max(arrive_time) 

            #Plot image
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     
            ax1.plot(torg, data, 'k', linewidth=0.5)
            ax1.set_title(title)
            ax1.margins(x=0)
            ax1.set_position([0.125, 0.6, 0.775, 0.3])  

            cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
            ax2.set_xlabel('Time (s)')
            f0lab = f0_array
            ax2.axvline(x=tprime0, c = '#377eb8', ls = '--', linewidth=0.7,label='Estimated arrival: '+str(np.round(tprime0,2))+' s')

            for pp in range(len(peaks)):
                tprime = freqpeak[pp]

                ft = calc_ft(times, tprime0, f0, v0, l, c)
                ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) #(0,(5,10)),
                
                if np.abs(tprime -tprime0) < 1.5:
                    ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 

            f0lab_sorted = sorted(f0lab)
            covm = np.sqrt(np.diag(covm))

            if len(f0lab_sorted) <= 17:
                fss = 'medium'
            else:
                fss = 'small'
            ax2.set_title("Final Model:\nt0'= "+str(np.round(tprime0,2))+' +/- ' + str(np.round(covm[3],2)) + ' sec, v0 = '+str(np.round(v0,2))+' +/- ' + str(np.round(covm[1],2)) +' m/s, l = '+str(np.round(l,2))+' +/- ' + str(np.round(covm[2],2)) +' m, \n' + 'f0 = '+str(f0lab_sorted)+' +/- ' + str(np.round(covm[0],2)) +' Hz', fontsize=fss)
            ax2.axvline(x=tarrive, c = '#e41a1c', ls = '--',linewidth=0.5,label='Wave arrvial: '+str(np.round(tarrive,2))+' s')

            ax2.legend(loc='upper right',fontsize = 'x-small')
            ax2.set_ylabel('Frequency (Hz)')

            ax2.margins(x=0)
            ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

            plt.colorbar(mappable=cax, cax=ax3)
            ax3.set_ylabel('Relative Amplitude (dB)')

            ax2.margins(x=0)
            ax2.set_xlim(0, 240)
            ax2.set_ylim(0, int(fs/2))

            # Plot overlay
            spec2 = 10 * np.log10(MDF)
            middle_column2 = spec2[:, middle_index]
            vmin2 = np.min(middle_column2)
            vmax2 = np.max(middle_column2)

            # Create ax4 and plot on the same y-axis as ax2
            ax4 = fig.add_axes([0.125, 0.11, 0.07, 0.35], sharey=ax2) 
            ax4.plot(middle_column2, frequencies, c='#ff7f00')  
            ax4.set_ylim(0, int(fs/2))
            ax4.set_xlim(vmax2*1.1, vmin2) 
            ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
            ax4.grid(axis='y')
            plt.show()
            BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5plane_spectotal/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'
            make_base_dir(BASE_DIR)
            fig.savefig('/scratch/irseppi/nodal_data/plane_info/5plane_spectotal/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(time[n])+'_'+str(flight_num[n])+'.png')
            plt.close()
            '''
