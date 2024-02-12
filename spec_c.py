import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from scipy import signal
from obspy.geodetics import gps2dist_azimuth
import math
import scipy.io
import os
import obspy
from obspy.core import UTCDateTime
import datetime
from numpy.fft import fft, ifft
from pathlib import Path
from math import radians, sin, cos, sqrt, atan2
from geopy.distance import gps2dist_azimuth


def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()

def closest_encounter(aircraft_data, seismometer_location):
    closest_distance = float('inf')
    closest_time = None
    closest_altitude = None

    for timestamp, aircraft_location, altitude in aircraft_data:
        _, _, distance = gps2dist_azimuth(aircraft_location[0], aircraft_location[1], seismometer_location[0], seismometer_location[1])
        distance = distance / 1000  # convert distance to kilometers
        if distance < closest_distance:
            closest_distance = distance
            closest_time = timestamp
            closest_altitude = altitude

    return closest_time, closest_distance, closest_altitude


def calculate_wave_arrival(closest_time, closest_distance, closest_altitude, aircraft_speed, aircraft_heading, seismometer_heading):
    speed_of_sound = 343  # speed of sound in m/s at sea level
    aircraft_speed_mps = aircraft_speed * 0.514  # convert aircraft speed from knots to m/s

    # calculate the component of the aircraft's speed in the direction of the seismometer
    relative_heading = radians(aircraft_heading - seismometer_heading)
    relative_speed = aircraft_speed_mps * cos(relative_heading)

    # calculate the time it takes for the sound to travel from the aircraft to the seismometer
    time_for_sound_to_travel = sqrt((closest_distance * 1000)**2 + (closest_altitude * 0.3048)**2) / (speed_of_sound + relative_speed)

    wave_arrival_time = closest_time + datetime.timedelta(seconds=time_for_sound_to_travel)

    return wave_arrival_time

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
sta = seismo_data['Station']

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
	for line in range(len(tm)):
		if str(tm[line]) == str(time[n]):
			speed = flight_data['speed'][line]
			alt = flight_data['altitude'][line]
			for y in range(len(station)):
				if str(station[y]) == str(sta[n]):
					closest_encounter(aircraft_data, seismometer_location)	
						
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
					fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))     
					
					ax1.plot(t, data, 'k', linewidth=0.5)
					ax1.set_title(title)
					ax1.axvline(x=tim, c = 'r', ls = '--')
					ax1.margins(x=0)
					
					spec = 10 * np.log10(Sxx) - (10 * np.log10(MDF))

					middle_index = len(times) // 2
					middle_column = spec[:, middle_index]
					vmin = 0
					vmax = np.max(middle_index)

					# Plot spectrogram
					cax = ax2.pcolormesh(times, frequencies, spec, shading='gouraud', cmap='hot_r', vmin=vmin, vmax=vmax)
					ax2.set_xlabel('Time [s]')

					# Find the center of the trace
					center_index = len(data) // 2
					center_time = t[center_index]
					
					l = np.sqrt(dist**2+(alt*0.0003048)**2) #km
					vo = speed*0.000514444
					xloc = calc_time(center_time, l, vo)

					ax2.axvline(x=center_time, c = 'r', ls = '--')
					ax2.axvline(x=xloc, c = 'k', ls = '--')
					ax2.set_ylabel('Frequency (Hz)')
					ax2.margins(x=0)
					ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

					plt.colorbar(mappable=cax, cax=ax3)
					ax3.set_ylabel('Relative Amplitude (dB)')
					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-02-'+str(day[n])+ '/'+str(flight_num[n])+'/'+station[y]+'/'
					make_base_dir(BASE_DIR)

					fig.savefig('/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-02-'+str(day[n])+ '/'+str(flight_num[n])+'/'+station[y]+'/'+str(time[n])+'_'+str(flight_num[n])+'.png')
					plt.close()

					# Find the index of the middle frequency
					middle_index = len(times) // 2

					# Extract the middle line of the spectrogram
					middle_column = spec[:, middle_index]
					
					peaks, _ = signal.find_peaks(middle_column, prominence=20) #, distance = 10) 
					np.diff(peaks)
					fig = plt.figure(figsize=(10,6))
					plt.grid()
					
					plt.plot(frequencies, middle_column, c='c')
					plt.plot(peaks, middle_column[peaks], "x")
					for g in range(len(peaks)):
						plt.text(peaks[g], middle_column[peaks[g]], peaks[g])
					plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))

					plt.xlim(0,int(fs/2))
					plt.ylim(vmin,vmax)
					plt.xlabel('Freq [Hz]')
					plt.ylabel('Amplitude [dB]')
					plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))
					plt.show()
					#BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5spec/201902'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'/'
					#make_base_dir(BASE_DIR)
					#fig.savefig('/scratch/irseppi/nodal_data/plane_info/5spec/201902'+str(day[n])+'/'+str(flight_num[n])+'/'+station[y]+'/'+station[y]+'_'+str(time[n])+'.png')
					plt.close()
					
