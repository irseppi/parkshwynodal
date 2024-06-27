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
                        fnot = [153, 172, 228, 38, 57, 76, 93, 135]
                        v0 = 68
                        l = 2135
                        f0 = 115

                    if n == 1:
                        fnot = [37, 56, 73,  127, 146, 165, 182, 218, 238, 256, 275]
                        tprime0 = 107
                        f0 = 110
                        v0 = 100
                        l = 2700

                    if n == 2:
                        fnot = [79, 120, 261]
                        tprime0 = 93
                        f0 = 131
                        v0 = 139
                        l = 4650

                    if n == 3:
                        fnot = [36,69,104,136,144]
                        tprime0 = 116
                        f0 = 121
                        v0 = 142
                        l = 2450

                    if n == 4:
                        fnot = [13,27,40,52,67,80,92,108,134,148,160,173,186,202,223,239,247,270]
                        tprime0 = 140
                        f0 = 120
                        v0 = 67
                        l = 5800
                   
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
                    corridor_width = (fs/2) / len(p)

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
                    for t_f in range(len(times)):
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
                                coord_inv.append((times[t_f], max_amplitude_frequency))
                                ttt.append(times[t_f])
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
                        if np.abs(delf[i]) <= (3):
                            fobs.append(maxfreq[i])
                            tobs.append(ttt[i])
                            count += 1
                    
                    peaks_assos.append(count)
                    
                    for i in range(len(fnot)):
                        f0 = fnot[i]
                        f0_array.append(f0)
                        m0 = [f0, vnot, lnot, tprimenot]
                        ft = calc_ft(times,  tprime00, f0, v00, l0, c)
                        
                        maxfreq = []
                        coord_inv = []
                        ttt = []
                        for t_f in range(len(times)):
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
                                    coord_inv.append((times[t_f], max_amplitude_frequency))
                                    ttt.append(times[t_f])
                                except:
                                    continue
                                
                        coord_inv_array = np.array(coord_inv)
                        if len(coord_inv_array) == 0:
                            continue
                        m,_ = invert_f(m0, coord_inv_array, num_iterations=4)
                        
                        ft = calc_ft(ttt, m[3], m[0], m[1], m[2], c)

                        delf = np.array(ft) - np.array(maxfreq)

                        count = 0
                        for i in range(len(delf)):
                            if np.abs(delf[i]) <= (3):
                                fobs.append(maxfreq[i])
                                tobs.append(ttt[i])
                                count += 1
                        
                        peaks_assos.append(count)
                    plt.scatter(tobs, fobs, color='black', marker='x')

                    plt.show()
                    qv = 0
                    num_iterations = 20
                    
                    w  = len(f0_array)
                    cprior = np.zeros((w+3,w+3))
                    m0 = []
                    m0.append(vnot)
                    m0.append(lnot)
                    m0.append(tprimenot)
                    for i in range(w):
                        m0.append(f0_array[i])
                    for row in range(len(cprior)):
                        if row == 0:
                            cprior[row][row] = 20**2
                        elif row == 1:
                            cprior[row][row] = 300**2
                        elif row == 2:
                            cprior[row][row] = 20**2
                        else:
                            cprior[row][row] = 1**2
                    
                    Cd = np.zeros((len(fobs), len(fobs)), int)
                    np.fill_diagonal(Cd, 3**2)
                    mnew = np.array(m0)
                    Sd_vec = np.zeros(num_iterations)
                    Sm_vec = np.zeros(num_iterations)
                    S_vec = np.zeros(num_iterations)
                    iter_vec = np.transpose(np.arange(0,num_iterations))
                    plt.figure()
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
                        
                        #m = np.array(m0) + cprior@G.T@la.inv(G@cprior@G.T+Cd)@(np.array(fobs)- np.array(fnew))
                        print('iteration %i out of %i' % (qv,num_iterations))
                        
                        dpred  = np.array(fnew)
                        Gm     = G
                        icobs  = la.inv(Cd)
                        dobs   = np.array(fobs)
                        
                        icprior = la.inv(cprior)
                        # steepest ascent vector (Eq. 6.307 or 6.312)
                        gamma = cprior @ Gm.T @ icobs @ (dpred - dobs) + (mnew - mprior)  # steepest ascent vector

                        #===================================================
                        # QUASI-NEWTON ALGORITHM (Eq. 6.319, nu=1)
                        
                        # approximate curvature
                        H = np.identity((w+3)) + cprior @ Gm.T @ icobs @ Gm
                        
                        dm   = -inv(H) @ gamma
                        
                        m = m + dm
                        #===================================================
                        
                        # misfit function for new model
                        # note: bookkeeping only -- not used within the algorithm above
                        Sd_vec[qv] = Sd(mnew,dobs,dpred,icobs)
                        Sm_vec[qv] = Sm(mnew,mprior,icprior)
                        S_vec[qv]  = S(mnew,dobs,dpred,mprior,icobs,icprior)
                        print(Sd_vec[qv],Sm_vec[qv],S_vec[qv])
                                
                        #covmlsq = la.inv(G.T@la.inv(Cd)@G + la.inv(cprior))
                        mnew = m
                        v0 = mnew[0]
                        l = mnew[1]
                        tprime0 = mnew[2]
                        f0_array = mnew[3:]
                        plt.scatter(v0,l)
                        
                        print(mnew)
                        qv += 1
                    plt.show()
                    plt.figure()
                    plt.plot(iter_vec,np.log10(S_vec),'k.-',iter_vec,np.log10(Sm_vec),'b.-',iter_vec,np.log10(Sd_vec),'r.-',
                        linewidth=2,markersize=20)
                    plt.legend(('S(mⁿ) = Sd + Sm','Sm(mⁿ)','Sd(mⁿ)'),fontsize=14)
                    plt.xlim([-0.5, num_iterations+0.5])

                    plt.locator_params(axis="x", integer=True, tight=True)
                    plt.xlabel('n, iteration')
                    plt.ylabel(' log10[ S(mⁿ) ], misfit function')
                    plt.title(str(num_iterations) + ' iterations')
                    plt.show()
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
                    
                    for pp in range(len(f0_array)):
                        f0 = f0_array[pp]
                        
                        ft = calc_ft(times, tprime0, f0, v0, l, c)

                        ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) #(0,(5,10)),
                        
                        if np.abs(tprime -tprime0) < 1.5:
                            ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 
                        f0lab.append(int(f0)) 
                        what_if = calc_ft(times, tarrive, fs/4, speed_mps, np.sqrt(dist_m**2 + alt_m**2), c)
                        #ax2.plot(times, what_if, 'red', ls = '--', linewidth=0.4)
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

