import numpy as np
import matplotlib.pyplot as plt
from prelude import *

##############################################################################################################################################################################################################

def remove_median(Sxx):
    """
    Remove the median from the spectrogram.

    Args:
        Sxx (array): The spectrogram data.

    Returns:
        spec: The spectrogram data with the median removed
        MDF: The median removed from the spectrogram
    """
    a, b = Sxx.shape

    MDF = np.zeros((a,b))
    for row in range(len(Sxx)):
        m = len(Sxx[row])
        p = sorted(Sxx[row])
        median = p[int(m/2)]
        for col in range(m):
            MDF[row][col] = median
    spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))
    return spec, MDF

##############################################################################################################################################################################################################

def plot_spectrgram(data, fs, torg, title, spec, times, frequencies, tprime0, v0, l, c, f0_array, F_m, arrive_time, MDF, covm, flight, middle_index, tarrive, closest_time, dir_name, plot_show=True):
    """
    Plot the spectrogram and spectrum of the given data.

    Args:
        data (array): The waveform data.
        fs (int): The sampling frequency.
        torg (array): The time array.
        title (str): The title of the plot.
        spec (array): The spectrogram data.
        times (array): The time array for the spectrogram.
        frequencies (array): The frequency array for the spectrogram.
        tprime0 (float): The estimated arrival time.
        v0 (float): The velocity.
        l (float): The distance.
        c (float): The speed of sound.
        f0_array (array): The array of frequencies.
        arrive_time (array): The arrival time array.
        MDF (array): Median removed from spectrogram.
        covm (array): The covariance matrix.
        flight (int): The flight number.
        middle_index (int): The index of the middle column.
        closest_time (float): The closest time.
        dir_name (str): The directory name.

    Returns:
        str: The user assigned quality number.
    """
    # Plot settings and calculations
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

        ax2.plot(times, ft, '#377eb8', ls = (0,(5,20)), linewidth=0.7) 
        tprime = tprime0
        t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
        ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
        
        ax2.scatter(tprime0, ft0p, color='black', marker='x', s=30) 

    fss = 'x-small'
    f0lab = sorted(f0_array)
    for i in range(len(f0lab)):
        f0lab[i] = (str(np.round(f0lab[i],2)))
    ax2.set_title("t0'= "+str(np.round(tprime0,2)) + ' sec, v0 = '+str(np.round(v0,2)) +' m/s, l = '+str(np.round(l,2)) +' m, \n' + 'f0 = '+', '.join(f0lab) +' Hz\nMisfit: ' + str(np.round(F_m,2)), fontsize=fss)
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
    if plot_show:
        plt.show()     
        qnum = input('What quality number would you give this?(first num for data quality(0-3), second for ability to fit model to data(0-1))')
    else:
        qnum = '__'
    fig.savefig(dir_name+'/'+str(closest_time)+'_'+str(flight)+'.png')
    plt.close()
    print(tprime0,v0,l,f0lab,covm)
    return qnum

##############################################################################################################################################################################################################

def plot_spectrum(spec, frequencies, tprime0, v0, l, c, f0_array, arrive_time, fs, closest_index, closest_time, sta, dir_name):
    """
    Plot the spectrum with arrival time markers.

    Args:
        spec (numpy.ndarray): The spectrum data.
        frequencies (numpy.ndarray): The frequencies.
        tprime0 (float): The reference arrival time.
        v0 (float): The velocity.
        l (float): The distance.
        c (float): The speed of light.
        f0_array (list): The list of frequencies.
        arrive_time (numpy.ndarray): The arrival times.
        fs (int): The sampling frequency.
        closest_index (int): The index of the closest time.
        closest_time (float): The closest time.
        sta (int or str): The station identifier.
        dir_name (str): The directory name.

    Returns:
        None
    """
    vmax = np.max(arrive_time)
    fig = plt.figure(figsize=(10,6))
    plt.grid()

    plt.plot(frequencies, spec[:,closest_index], c='#377eb8')
        
    for pp in range(len(f0_array)):
        f0 = f0_array[pp]
        if fs/2 < f0:
            continue
        tprime = tprime0
        t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
        ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
        if ft0p == np.nan:
            continue
        if ft0p > 250:
            continue
        upper = int(ft0p + 10)
        lower = int(ft0p - 10)
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


    fig.savefig(dir_name + '/'+str(sta)+'_' + str(closest_time) + '.png')
    plt.close()

##############################################################################################################################################################################################################

def doppler_picks(spec, times, frequencies, vmin, vmax, month, day, flight, sta, equip, closest_time, start_time, make_picks=True):
    """
    Pick the points for the doppler shift.

    Args:
        spec (numpy.ndarray): The spectrogram data.
        times (numpy.ndarray): The time array.
        frequencies (numpy.ndarray): The frequency array.
        vmin (float): The minimum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        vmax (float): The maximum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        month (int): The month of the data.
        day (int): The day of the data.
        flight (int): The flight number.
        sta (int or str): The station identifier.
        closest_time (float): The time of closest approach.

    Returns:
        list: The list of picks the user picked along the most prominent overtone.
    """
    file_name = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight)+'.csv'
                
    if Path(file_name).exists():
        coords = []
        with open(file_name, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                coords.append((float(pick_data[0]), float(pick_data[1])))
        file.close()  
        return coords
    elif make_picks:
        BASE_DIR = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/inversepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'
        make_base_dir(BASE_DIR)
        pick_again = 'y'
        while pick_again == 'y':
            r1 = open(file_name,'w')
            coords = []
            plt.figure()
            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
            def onclick(event, coords=coords):
                #global coords
                coords.append((event.xdata, event.ydata))
                plt.scatter(event.xdata, event.ydata, color='black', marker='x')  
                plt.draw() 
                print('Clicked:', event.xdata, event.ydata)  
                r1.write(str(event.xdata) + ',' + str(event.ydata) + ',' + str(start_time) + ',\n')
            plt.gcf().canvas.mpl_connect('button_press_event', onclick)

            plt.show(block=True)
            r1.close()
            pick_again = input("Do you want to repick your points? (y or n)")
        return coords
    else:
        return []

##############################################################################################################################################################################################################

def overtone_picks(spec, times, frequencies, vmin, vmax, month, day, flight, sta, equip, closest_time, start_time, tprime0, make_picks=True):
    """
    Pick the points for the overtone shift.

    Args:
        spec (numpy.ndarray): The spectrogram data.
        times (numpy.ndarray): The time array.
        frequencies (numpy.ndarray): The frequency array.
        vmin (float): The minimum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        vmax (float): The maximum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        month (int): The month of the data.
        day (int): The day of the data.
        flight (int): The flight number.
        sta (int or str): The station identifier.
        closest_time (float): The time of closest approach.
        tprime0 (float): The estimated arrival time.

    Returns:
        list: The list of peak amplitudes.
        list: The list of peak frequencies.
    """
    output2 = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/overtonepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight)+'.csv'
    if Path(output2).exists():

        peaks = []
        freqpeak = []
        with open(output2, 'r') as file:
            for line in file:
                pick_data = line.split(',')
                peaks.append(float(pick_data[1]))
                freqpeak.append(float(pick_data[0]))
        file.close()  
        return peaks, freqpeak
    elif make_picks:
        BASE_DIR = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/overtonepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'
        make_base_dir(BASE_DIR)
        pick_again = 'y'
        while pick_again == 'y':
            r2 = open(output2,'w')
            peaks = []
            freqpeak = []
            plt.figure()
            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
            plt.axvline(x=tprime0, c = '#377eb8', ls = '--')
            plt.axvline(x=120, c = '#e41a1c', ls = '--')
            def onclick(event):
                #global coords
                peaks.append(event.ydata)
                freqpeak.append(event.xdata)
                plt.scatter(event.xdata, event.ydata, color='black', marker='x')  # Add this line
                plt.draw() 
                print('Clicked:', event.xdata, event.ydata)  
                r2.write(str(event.xdata) + ',' + str(event.ydata) + ',' + str(start_time) + ',\n')
            cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)

            plt.show(block=True)
            r2.close()
            pick_again = input("Do you want to repick you points? (y or n)")
        
        return peaks, freqpeak
    else:
        return [], []

##############################################################################################################################################################################################################

def time_picks(month, day, flight, sta, equip, tobs, fobs, closest_time, start_time, spec, times, frequencies, vmin, vmax, w, peaks_assos, make_picks=True):
    """
    Pick the points for the time shift.

    Args:

        month (int): The month of the data.
        day (int): The day of the data.
        flight (int): The flight number.
        sta (int or str): The station identifier.
        tobs (list): The time array.
        fobs (list): The frequency array.
        closest_time (float): The time of closest approach.
        spec (numpy.ndarray): The spectrogram data.
        times (numpy.ndarray): The time array.
        frequencies (numpy.ndarray): The frequency array.
        vmin (float): The minimum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        vmax (float): The maximum amplitude value for the center line of the spectrogram. Used for adjusting colorbar.
        w (int): The number of peaks.

    Returns:
        list: The time array.
        list: The frequency array.
    """
    output3 = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/timepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight)+'.csv'
    if Path(output3).exists():
        set_time = []
        with open(output3, 'r') as file:
            for line in file:

                pick_data = line.split(',')
                set_time.append(float(pick_data[0]))
        file.close()  
        if len(set_time) == 0:
            set_time = [0, 250]
        start_time = set_time[0]
        end_time = set_time[1]

        ftobs = []
        ffobs = []
        if peaks_assos == False:
            for j in range(len(tobs)):
                if tobs[j] >= start_time and tobs[j] <= end_time:
                    ftobs.append(tobs[j])
                    ffobs.append(fobs[j])
            peaks_assos = np.nan
        else:
            peak_ass = []
            cum = 0

            for p in range(w):
                count = 0
                for j in range(cum,cum+peaks_assos[p]):
                    if tobs[j] >= start_time and tobs[j] <= end_time:
                        ftobs.append(tobs[j])
                        ffobs.append(fobs[j])
                        count += 1
                cum = cum + peaks_assos[p]
            
                peak_ass.append(count)
            peaks_assos = peak_ass
        tobs = ftobs
        fobs = ffobs

        return tobs, fobs, peaks_assos

    elif make_picks:
        BASE_DIR = '/home/irseppi/REPOSITORIES/parkshwynodal/output/' + equip + '_data_picks/timepicks/2019-0'+str(month)+'-'+str(day)+'/'+str(flight)+'/'+str(sta)+'/'
        make_base_dir(BASE_DIR)
        
        pick_again = 'y'
        while pick_again == 'y':
            r3 = open(output3,'w')
            set_time = []
            plt.figure()
            plt.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='pink_r', vmin=vmin, vmax=vmax)
            plt.scatter(tobs,fobs, color='black', marker='x')
            def onclick(event):
                #global coords
                set_time.append(event.xdata) 
                plt.scatter(event.xdata, event.ydata, color='red', marker='x')  # Add this line
                plt.draw() 
                print('Clicked:', event.xdata, event.ydata)  
                r3.write(str(event.xdata) + ',' + str(event.ydata) + ',' + str(start_time) + ',\n')

            cid = plt.gcf().canvas.mpl_connect('button_press_event', onclick)
            plt.show(block=True)

            r3.close()
            pick_again = input("Do you want to repick you points? (y or n)")
        if len(set_time) == 0:
            set_time = [0, 250]
        start_time = set_time[0]
        end_time = set_time[1]
        ftobs = []
        ffobs = []
        if peaks_assos == False:
            for j in range(len(tobs)):
                if tobs[j] >= start_time and tobs[j] <= end_time:
                    ftobs.append(tobs[j])
                    ffobs.append(fobs[j])
            peaks_assos = np.nan
        else:
            peak_ass = []
            cum = 0

            for p in range(w):
                count = 0
                for j in range(cum,cum+peaks_assos[p]):
                    if tobs[j] >= start_time and tobs[j] <= end_time:
                        ftobs.append(tobs[j])
                        ffobs.append(fobs[j])
                        count += 1
                cum = cum + peaks_assos[p]
            
                peak_ass.append(count)
            peaks_assos = peak_ass
        tobs = ftobs
        fobs = ffobs

        return tobs, fobs, peaks_assos
    else:
        return tobs, fobs, []
