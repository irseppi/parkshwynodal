import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import datetime
from prelude import *
from scipy.signal import spectrogram
import scipy.linalg as la

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station = seismo_data['Station']
flight_num = [530342801,528485724,528473220,528407493,528293430,531605202,531715679,529805251,529948401] 
time = [1551066051,1550172833,1550168070,1550165577,1550089044,1551662362,1551736354,1550803701,1550867033] 
sta = [1022,1272,1173,1283,1004,1010,1021,1006,1109]
day = [25,14,14,14,13,4,4,22,22]
month = [2,2,2,2,2,3,3,2,2]

for n in range(4,5): #9):
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
                        f0_array = [78,120,130, 258]
                        tprime0 = 93
                        v0 = 142
                        l = 4992

                    if n == 3:
                        f0_array = [34,69,104,119,134,139]
                        tprime0 = 115
                        v0 = 159
                        l = 3802

                    if n == 4:
                        f0_array = [14,28,41,54,68,81,95,109,123,136,148,161,177,189,203,226,241,249,271]
                        tprime0 = 140
                        v0 = 62
                        l = 500
                    else:
                        #f0_array =
                        tprime0 = tarrive
                        v0 = speed_mps
                        l = np.sqrt(dist_m**2 + alt_m**2)
                    #Choose start and end time
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
                    new_times = times[np.where((times >= start_time) & (times <= end_time))]

                    w  = len(f0_array)
                    
                    mprior = []
                    mprior.append(v0)
                    mprior.append(l)
                    mprior.append(tprime0)
                    for i in range(w):
                        mprior.append(f0_array[i])
                    mprior = np.array(mprior)
                    c = 343
                    
                    peaks_assos = []
                    fobs = []
                    tobs = []

                    plt.figure()
                    plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
                    corridor_width = (fs/2)/len(frequencies)
                    for i in range(w):
                        f0 = f0_array[i]
                        
                        ft = calc_ft(new_times, tprime0, f0, v0, l, c)

                        option = 2
                        if option == 1:
                            count = 0
                        if option == 2:
                            maxfreq = []
                            ttt = []
                            coord_inv = []
                        for t_f in range(len(new_times)):
                            if not np.isnan(ft[t_f]) and ft[t_f] != np.inf:
                                upper = int(ft[t_f] + corridor_width)
                                lower = int(ft[t_f] - corridor_width)

                                if lower < 0:
                                    lower = 0
                                if upper > 250:
                                    upper = 250
                                try:
                                    tt = spec[lower:upper, t_f]
                                    if option == 1:
                                        max_amplitude_index = np.argmax(tt)
                                        max_amplitude_frequency = frequencies[max_amplitude_index+lower]
                                        fobs.append(max_amplitude_frequency)
                                        tobs.append(new_times[t_f])
                                        count += 1
                                    if option == 2:
                                        max_amplitude_index = np.argmax(tt)
                                        max_amplitude_frequency = frequencies[max_amplitude_index+lower]
                                        maxfreq.append(max_amplitude_frequency)
                                        coord_inv.append((new_times[t_f], max_amplitude_frequency))
                                        ttt.append(new_times[t_f])
                                except:
                                    continue
                        if option == 2:
                            coord_inv_array = np.array(coord_inv)
                            mtest = [f0,v0, l, tprime0]
                            mtest,_ = invert_f(mtest, coord_inv_array, num_iterations=12)
                            ft = calc_ft(ttt,  mtest[3], mtest[0], mtest[1], mtest[2], c)
                            delf = np.array(ft) - np.array(maxfreq)

                            count = 0
                            for i in range(len(delf)):
                                if np.abs(delf[i]) <= (2):
                                    fobs.append(maxfreq[i])
                                    tobs.append(ttt[i])
                                    count += 1
                            peaks_assos.append(count)
                        
                    plt.scatter(tobs, fobs, color='black', marker='x')

                    plt.show()
                    
                    qv = 0
                    num_iterations = 4

                    cprior = np.zeros((w+3,w+3))

                    for row in range(len(cprior)):
                        if row == 0:
                            cprior[row][row] = 20**2
                        elif row == 1:
                            cprior[row][row] = 500**2
                        elif row == 2:
                            cprior[row][row] = 20**2
                        else:
                            cprior[row][row] = 1**2
                    
                    Cd = np.zeros((len(fobs), len(fobs)), int)
                    np.fill_diagonal(Cd, 3**2)
                    mnew = np.array(mprior)

                    
                    while qv < num_iterations:
                        G = np.zeros((0,w+3))
                        fnew = []
                        cum = 0
                        for p in range(w):
                            new_row = np.zeros(w+3)
                            f0 = f0_array[p]
                          
                            if p == 0:
                                for j in range(cum,peaks_assos[p]):
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
                                for j in range(cum,cum+peaks_assos[p]):
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
                        
                            cum = cum + peaks_assos[p]
 
                        m = np.array(mnew) + cprior@G.T@la.inv(G@cprior@G.T+Cd)@(np.array(fobs)- np.array(fnew))
                        mnew = m
                        v0 = mnew[0]
                        l = mnew[1]
                        tprime0 = mnew[2]
                        f0_array = mnew[3:]

                        print(m)
                        qv += 1
                    covm = la.inv(G.T@la.inv(Cd)@G + la.inv(cprior))
                    
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
                    ax1.set_position([0.125, 0.6, 0.775, 0.3])  # Move ax1 plot upwards


                    # Plot spectrogram
                    cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)				
                    ax2.set_xlabel('Time (s)')
                    f0lab = []
                    ax2.axvline(x=tprime0, c = '#377eb8', ls = '--', linewidth=0.7,label='Estimated arrival: '+str(np.round(tprime0,2))+' s')
                    f0_array = sorted(f0_array)
                    covm = np.sqrt(np.diag(covm))
                    for pp in range(len(f0_array)):
                        f0 = f0_array[pp]
                        
                        ft = calc_ft(times, tprime0, f0, v0, l, c)

                        ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) 
                        tprime = tprime0
                        t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
                        ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
                        
                        ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 
                        f0lab.append(str(np.round(f0,2)) +'+/-' + str(np.round(covm[3+pp],2))+',') 

                        if len(f0_array) > 16:
                            if pp == int(len(f0_array)/3):
                                f0lab.append('\n')
                            elif pp == int(len((f0_array)/3)+len((f0_array)/3)):
                                f0lab.append('\n')
                        elif len(f0_array) > 8:
                                if pp == int(len(f0_array)/2):
                                    f0lab.append('\n')

                    fss = 'x-small'
                    
                    ax2.set_title("Final Model:\nt0'= "+str(np.round(tprime0,2))+' +/- ' + str(np.round(covm[2],2)) + ' sec, v0 = '+str(np.round(v0,2))+' +/- ' + str(np.round(covm[0],2)) +' m/s, l = '+str(np.round(l,2))+' +/- ' + str(np.round(covm[1],2)) +' m, \n' + 'f0 = '+' '.join(f0lab) +' Hz', fontsize=fss)
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
                    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'
                    make_base_dir(BASE_DIR)
                    fig.savefig('/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-0'+str(month[n])+'-'+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(time[n])+'_'+str(flight_num[n])+'.png')
                    plt.close()
                    
                    fig = plt.figure(figsize=(10,6))
                    plt.grid()

                    plt.plot(frequencies, arrive_time, c='#377eb8')
                       
                    for pp in range(len(f0_array)):
                        f0 = f0_array[pp]
                        if fs/2 < f0:
                            continue
                        tprime = tprime0
                        t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
                        ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))

                        upper = int(ft0p + 5)
                        lower = int(ft0p - 5)
                        tt = spec[lower:upper, closest_index]
                        if upper > 250:
                            freqp = ft0p
                            ampp = np.interp(ft0p, frequencies, arrive_time)
                        elif lower < 0:
                            freqp = ft0p
                            ampp = np.interp(ft0p, frequencies, arrive_time)
                        else:
                            ampp = np.max(tt)
                            freqp = np.argmax(tt)+lower
                        plt.scatter(freqp, ampp, color='black', marker='x', s=100)
                        if isinstance(sta[n], int):
                            plt.text(freqp - 5, ampp + 0.8, freqp, fontsize=17, fontweight='bold')
                        else:
                            plt.text(freqp - 1, ampp + 0.8, freqp, fontsize=17, fontweight='bold')  
                    plt.xlim(0, int(fs/2))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    plt.ylim(0,vmax*1.1)
                    plt.xlabel('Frequency (Hz)', fontsize=17)
                    plt.ylabel('Relative Amplitude at t = {:.2f} s (dB)'.format(tprime0), fontsize=17)


                    make_base_dir('/scratch/irseppi/nodal_data/plane_info/5spec/20190'+str(month[n])+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/')

                    fig.savefig('/scratch/irseppi/nodal_data/plane_info/5spec/20190'+str(month[n])+str(day[n])+'/'+str(flight_num[n])+'/'+str(sta[n])+'/'+str(sta[n])+'_' + str(time[n]) + '.png')
                    plt.close() 