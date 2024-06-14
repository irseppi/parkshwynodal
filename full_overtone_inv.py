import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import datetime
from prelude import *
from scipy.signal import find_peaks, spectrogram

def df_mult(f0,v0,l,tp0,tp):   
    """
	Calculate the derivatives of f with respect to f0, v0, l, and tp0.

	Parameters:
	f0 (float): Fundamental frequency produced by the aircraft.
	v0 (float): Velocity of the aircraft.
	l (float): Distance of closest approach between the station and the aircraft.
	tp0 (float): Time of that the central frequency of the overtones occur, when the aircraft is at the closest approach to the station.
	tp (float): Array of times.
	Returns:
	tuple: A tuple containing the derivatives of f with respect to f0, v0, l, and tp0.
	"""

    c = 343 # m/sec speed of sound
    for f0 in f0:
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


        #derivative of f with respect to tp0
        f_derivetprime0 = ((f0 * l**2 * (c - v0) * v0**2 * (c + v0) * ((-tp + tp0) * v0**2 + c**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))) / 
        (c * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
        (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) - 
        c**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + 
        (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2))**2))


    return f_derivef0, f_derivev0, f_derivel, f_derivetprime0

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

                    c = 343
                    f0 = 120
                    m0 = [120, v0, l, tprime0]

                    m,covm = invert_f(m0, coords_array, num_iterations=8)

                    f0 = m[0]
                    v0 = m[1]
                    l = m[2]
                    tprime0 = m[3]
                    
                    ft = calc_ft(times, tprime0, f0, v0, l, c)

                    p, _ = find_peaks(middle_column, distance = 7)
                    corridor_width = (fs/2) / len(p) 
                    output2 = '/scratch/irseppi/nodal_data/plane_info/overtonepicks/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(time[n])+'_'+str(flight_num[n])+'.csv'
                    peaks = []
                    freqpeak = []
                    with open(output2, 'r') as file:
                        for line in file:
                            # Split the line using commas
                            pick_data = line.split(',')
                            peaks.append(float(pick_data[1]))
                            freqpeak.append(float(pick_data[0]))
                    file.close()  # Close the file after reading

                    num_iterations = 12
                    while n < num_iterations:
                        fnew = []

                        w  = len(peaks)
                        G = np.zeros((0,w+3))
                        for i in range(len(peaks)):
                            maxfreq = []
                            f00 = calc_f0(freqpeak[i],tprime0, peaks[i], v0, l, c)
                            v00 = m[1]
                            l0 = m[2]
                            tprime00 = m[3]
                            m0 = [f00, v00, l0, tprime00]
                            ft = calc_ft(times, tprime00, f00, v00, l0, c)
                            coord_inv = []
                            ttt = []
                            for t_f in range(len(times)):
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
                                    coord_inv.append((times[t_f], max_amplitude_frequency))
                                    ttt.append(times[t_f])
                                except:
                                    continue
                            coord_inv_array = np.array(coord_inv)
                            m,_ = invert_f(m0, coord_inv_array, num_iterations=12)
                            ft = calc_ft(ttt, tprime00, f00, v00, l0, c)

                            delf = np.array(ft) - np.array(maxfreq)
                            
                            new_coord_inv_array = []
                            for i in range(len(delf)):
                                if np.abs(delf[i]) <= 3:
                                    new_coord_inv_array.append(coord_inv_array[i])
                            coord_inv_array = np.array(new_coord_inv_array)
                            tobs = coord_inv_array[:,0]
                            fobs = coord_inv_array[:,1]
                            for j in range(0,len(coord_inv_array)):
                                tprime = tobs[j]
                                new_row = np.zeros(w+3)
                                new_row[0] = l
                                new_row[1] = v0
                                new_row[2] = tprime0
                                t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
                                ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
                                for p in range(3,3+w):
                                    f_derivef0, f_derivev0, f_derivel, f_derivetprime0 = df(m[0], m[1], m[2], m[3], tobs[j])
                                    new_row[0] = f_derivel
                                    new_row[1] = f_derivev0
                                    new_row[2] = f_derivetprime0
                                    new_row[w] = f_derivef0
                                G = np.vstack((G, new_row))

                                fnew.append(ft0p) 
                        print(G)
                        sigma = 5
                        m0 = []
                        for i in range(w):
                            m0.append(peaks[i])
                            
                        m0.append(v0)
                        m0.append(l)
                        m0.append(tprime0)
                        m = np.reshape(np.reshape(m0,(3+w,1))+ np.reshape(inv(G.T@G)@G.T@(np.reshape(fobs, (len(coords_array), 1)) - np.reshape(np.array(fnew), (len(coords_array), 1))), (3+w,1)), (3+w,))
                        covmlsq = (sigma**2)*la.inv(G.T@G)
                        
                        m0 = m
                        n += 1

                    closest_index = np.argmin(np.abs(tprime0 - times))
                    arrive_time = spec[:,closest_index]
                    for i in range(len(arrive_time)):
                        if arrive_time[i] < 0:
                            arrive_time[i] = 0
                    vmin = np.min(arrive_time) 
                    vmax = np.max(arrive_time) 

                    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=False, figsize=(8,6))     
                    ax1.plot(torg, data, 'k', linewidth=0.5)
                    ax1.set_title(title)
                    ax1.margins(x=0)
                    ax1.set_position([0.125, 0.6, 0.775, 0.3])  

                    cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
                    ax2.set_xlabel('Time (s)')
                    f0lab = []
                    ax2.axvline(x=tprime0, c = '#377eb8', ls = '--', linewidth=0.7,label='Estimated arrival: '+str(np.round(tprime0,2))+' s')
                    
                    for pp in range(len(peaks)):
                        tprime = freqpeak[pp]
                        ft0p = peaks[pp]
                        f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
                        ft = calc_ft(times, tprime0, f0, v0, l, c)
                        ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) #(0,(5,10)),
                        
                        if np.abs(tprime -tprime0) < 1.5:
                            ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 
                        f0lab.append(int(f0)) 

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
                    
                    fig = plt.figure(figsize=(10,6))
                    plt.grid()

                    plt.plot(frequencies, arrive_time, c='#377eb8')
                   
                    for pp in range(len(peaks)):
                        if np.abs(freqpeak[pp] -tprime0) < 1.5:
                            upper = int(peaks[pp] + 3)
                            lower = int(peaks[pp] - 3)
                            tt = spec[lower:upper, closest_index]
                            ampp = np.max(tt)
                            freqp = np.argmax(tt)+lower
                            plt.scatter(freqp, ampp, color='black', marker='x', s=100)

                            plt.text(freqp - 5, ampp + 0.8, freqp, fontsize=17, fontweight='bold')
 
                    plt.xlim(0, int(fs/2))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    plt.ylim(0,vmax*1.1)
                    plt.xlabel('Frequency (Hz)', fontsize=17)
                    plt.ylabel('Relative Amplitude at t = {:.2f} s (dB)'.format(tprime0), fontsize=17)


                    make_base_dir('/scratch/irseppi/nodal_data/plane_info/5spectotal/20190'+str(month[n])+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/')

                    fig.savefig('/scratch/irseppi/nodal_data/plane_info/5spectotal/20190'+str(month[n])+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(sta[n])+'_' + str(time[n]) + '.png')
                    plt.close() 