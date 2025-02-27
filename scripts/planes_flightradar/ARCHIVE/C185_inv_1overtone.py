import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import datetime
import pyproj
from prelude import *
from scipy.signal import find_peaks, spectrogram
from pathlib import Path

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

seismo_utm = [utm_proj(lon, lat) for lat, lon in zip(seismo_latitudes, seismo_longitudes)]
seismo_utm_x, seismo_utm_y = zip(*seismo_utm)

# Convert UTM coordinates to kilometers
seismo_utm_x_km = [x / 1000 for x in seismo_utm_x]
seismo_utm_y_km = [y / 1000 for y in seismo_utm_y]

seismo_utm_km = [(x, y) for x, y in zip(seismo_utm_x_km, seismo_utm_y_km)]

def closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, seismo_utm_x, seismo_utm_y):
    closest_point = None
    dist_lim = np.inf

    x = [flight_utm_x1, flight_utm_x2]
    y = [flight_utm_y1, flight_utm_y2]

    if (x[1] - x[0]) == 0:
        if (y[1]-y[0]) <= 0:
            ggg = -0.001
        else:
            ggg = 0.001
        for point in np.arange(y[0], y[1], ggg):
            xx = x[0]
            yy = point
            dist_km = np.sqrt((seismo_utm_y-yy)**2 +(seismo_utm_x-xx)**2)
            
            if dist_km < dist_lim:
                dist_lim = dist_km
                closest_point = (xx,yy)
            else:
                continue

    else: 
        m = (y[1]-y[0])/(x[1]-x[0])
        b = y[0] - m*x[0]

        if (x[1] - x[0]) <= 0:
            ggg = -0.001
        else:
            ggg = 0.001
        for point in np.arange(x[0], x[1], ggg):
            xx = point
        
            yy = m*xx + b
            dist_km = np.sqrt((seismo_utm_y-yy)**2 +(seismo_utm_x-xx)**2)
            
            if dist_km < dist_lim:
                dist_lim = dist_km
                closest_point = (xx,yy)
            else:
                continue

    return closest_point, dist_lim


def find_closest_point(flight_utm, seismo_utm):
    min_distance = np.inf
    closest_point = None

    for i in range(len(flight_utm) - 1):
        flight_utm_x1, flight_utm_y1 = flight_utm[i]
        flight_utm_x2, flight_utm_y2 = flight_utm[i + 1]
        seismo_utm_x, seismo_utm_y = seismo_utm
        point, d = closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, seismo_utm_x, seismo_utm_y)
        
        if point == None:
            continue
        elif d < min_distance:
            min_distance = d
            closest_point = point
            index = i
        else:
            continue
   
    return closest_point, min_distance, index

sta_f = open('input/all_station_crossing_db_C185.txt','r')
C185_output = open('output/C185data_1overtone2.csv', 'w')


for line in sta_f.readlines():
    val = line.split(',')
    date = val[0]
    flight = val[1]
    sta = val[6]

    spec_dir = '/scratch/irseppi/nodal_data/plane_info/C185_spec/2019-0'+str(date[5])+'-'+str(date[6:8])+'/'+str(flight)+'/'+str(sta)+'/'

    flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date) + '_positions/' + str(date) + '_' + str(flight) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    time = flight_data['snapshot_id']
    timestamp = flight_data['snapshot_id']
    speed = flight_data['speed']
    altitude = flight_data['altitude']

    # Convert flight latitude and longitude to UTM coordinates
    flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
    flight_utm_x, flight_utm_y = zip(*flight_utm)

    # Convert UTM coordinates to kilometers
    flight_utm_x_km = [x / 1000 for x in flight_utm_x]
    flight_utm_y_km = [y / 1000 for y in flight_utm_y]
    flight_path = [(x,y) for x, y in zip(flight_utm_x_km, flight_utm_y_km)]
    
    # Iterate over seismometer data
    for s in range(len(seismo_data)):
        if str(sta) == str(stations[s]):
            seismometer = (seismo_utm_x_km[s], seismo_utm_y_km[s])  

            closest_p, dist_km, index= find_closest_point(flight_path, seismometer)
         
            if dist_km <= 2:
                closest_x, closest_y = closest_p
                #Calculate the time of the closest point
                flight_utm_x1, flight_utm_y1 = flight_path[index]
                flight_utm_x2, flight_utm_y2 = flight_path[index + 1]

                x_timestamp_dif_vec = flight_utm_x2 - flight_utm_x1
                y_timestamp_dif_vec = flight_utm_y2 - flight_utm_y1

                cx_timestamp_dif_vec =  closest_x - flight_utm_x1
                cy_timestamp_dif_vec = closest_y - flight_utm_y1

                line_vector = (x_timestamp_dif_vec, y_timestamp_dif_vec)
                cline_vector = (cx_timestamp_dif_vec, cy_timestamp_dif_vec)

                line_magnitude = np.sqrt(line_vector[0] ** 2 + line_vector[1] ** 2)
                cline_magnitude = np.sqrt(cline_vector[0] ** 2 + cline_vector[1] ** 2)

                length_ratio = cline_magnitude / line_magnitude
                closest_time = timestamp[index] + length_ratio*(timestamp[index+1] - timestamp[index])

                alt = (altitude[index]+altitude[index+1])/2
                sp = (speed[index]+speed[index+1])/2

                alt_m = alt * 0.3048
                elevation = elevations[s]
                speed_mps = sp * 0.514444
                height_m = alt_m - elevation 
                dist_m = dist_km * 1000
                tmid = closest_time
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
                
                middle_index =  len(times) // 2
                middle_column = spec[:, middle_index]
                vmin = 0  
                vmax = np.max(middle_column) 

                c = 329
                tprime0 = 120
                v0 = speed_mps
                l = np.sqrt(dist_m**2 + (height_m)**2)

                tf = np.arange(0, 240, 1)

                output1 = '/home/irseppi/REPOSITORIES/parkshwynodal/output/C185_data_picks/inversepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight)+'.csv'

                coords = []
                with open(output1, 'r') as file:
                    for line in file:
                        # Split the line using commas
                        pick_data = line.split(',')
                        coords.append((float(pick_data[0]), float(pick_data[1])))
                file.close()  # Close the file after reading
                    
            
                if len(coords) == 0:

                    print('No picks for: ', date, flight, sta)
                    continue
                # Convert the list of coordinates to a numpy array
                coords_array = np.array(coords)

                f0 = 116
                c = 329
                m0 = [f0, v0, l, tprime0]

                m,covm = invert_f(m0, coords_array, num_iterations=8)
                f0 = m[0]
                v0 = m[1]
                l = m[2]
                tprime0 = m[3]
                if f0 == np.nan or v0 == np.nan or l == np.nan or tprime0 == np.nan:
                    continue
                output3 = '/home/irseppi/REPOSITORIES/parkshwynodal/output/C185_data_picks/timepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight)+'.csv'
                if Path(output3).exists():
                    set_time = []
                    with open(output3, 'r') as file:
                        for line in file:
                            # Split the line using commas
                            pick_data = line.split(',')
                            set_time.append(float(pick_data[0]))
                    file.close()  # Close the file after reading

                    start_time = set_time[0]
                    end_time = set_time[1]
                else:
                    start_time = 0
                    end_time = 240
                ft = calc_ft(times, tprime0, f0, v0, l, c)
                
                if isinstance(sta, int):

                    peaks = []
                    p, _ = find_peaks(middle_column, distance = 7)
                    corridor_width = 6

                    coord_inv = []

                    for t_f in range(len(times)):
                        if times[t_f] >= start_time and times[t_f] <= end_time:
                            try:
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
                            except:
                                continue
                        else:
                            continue

                    coord_inv_array = np.array(coord_inv)

                    m,_ = invert_f(m0, coord_inv_array, num_iterations=12)
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

                    m,covm = invert_f(m0, coord_inv_array, num_iterations=12, sigma=5)
                    
                    f0 = m[0]
                    v0 = m[1]
                    l = m[2]
                    tprime0 = m[3]

                
                output2 = '/home/irseppi/REPOSITORIES/parkshwynodal/output/C185_data_picks/overtonepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight)+'.csv'


                peaks = []
                freqpeak = []
                with open(output2, 'r') as file:
                    for line in file:
                        # Split the line using commas
                        pick_data = line.split(',')
                        peaks.append(float(pick_data[1]))
                        freqpeak.append(float(pick_data[0]))
                file.close()  # Close the file after reading
                
                f0_array = []
                f0_array.append(f0)
                for pp in range(len(peaks)):
                    tprime = freqpeak[pp]
                    ft0p = peaks[pp]
                    f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
                    f0_array.append(f0)
                f0_array = np.array(f0_array)

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
                covm = np.sqrt(np.diag(covm))
                for pp in range(len(f0_array)):
                    f0 = f0_array[pp]
                    
                    ft = calc_ft(times, tprime0, f0, v0, l, c)

                    ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) 
                    tprime = tprime0
                    t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
                    ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
                    
                    ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 

                fss = 'x-small'
                f0lab = sorted(f0_array)
                for i in range(len(f0lab)):
                    f0lab[i] = (str(np.round(f0lab[i],2)))
                ax2.set_title("Final Model:\nt0'= "+str(np.round(tprime0,2)) + ' sec, v0 = '+str(np.round(v0,2)) +' m/s, l = '+str(np.round(l,2)) +' m, \n' + 'f0 = '+', '.join(f0lab) +' Hz', fontsize=fss)
                ax2.axvline(x=120, c = '#e41a1c', ls = '--',linewidth=0.5,label='Wave arrvial: 120 s')
        
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


                BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/C185_spec_1o/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'
                make_base_dir(BASE_DIR)
                fig.savefig('/scratch/irseppi/nodal_data/plane_info/C185_spec_1o/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight)+'.png')
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
                    try:
                        upper = int(ft0p + 5)
                        lower = int(ft0p - 5)
                    except:
                        continue

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
                    if isinstance(sta, int):
                        plt.text(freqp - 5, ampp + 0.8, freqp, fontsize=17, fontweight='bold')
                    else:
                        plt.text(freqp - 1, ampp + 0.8, freqp, fontsize=17, fontweight='bold')  

                plt.xlim(0, int(fs/2))
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                plt.ylim(0,vmax*1.1)
                plt.xlabel('Frequency (Hz)', fontsize=17)
                plt.ylabel('Relative Amplitude at t = {:.2f} s (dB)'.format(tprime0), fontsize=17)

                make_base_dir('/scratch/irseppi/nodal_data/plane_info/C185_specrum_1o/20190'+str(month)+str(day)+'/'+str(flight)+'/'+str(sta)+'/')

                fig.savefig('/scratch/irseppi/nodal_data/plane_info/C185_specrum_1o/20190'+str(month)+str(day)+'/'+str(flight)+'/'+str(sta)+'/'+str(sta)+'_' + str(closest_time) + '.png')
                plt.close() 

                print(tprime0,v0,l,f0lab,covm)
                C185_output.write(str(date)+','+str(flight)+','+str(sta)+','+str(closest_time)+','+str(tprime0)+','+str(v0)+','+str(l)+','+str(f0_array)+','+str(covm)+',\n')
C185_output.close()
