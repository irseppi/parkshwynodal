import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.core import UTCDateTime
import datetime
import pytz
import obspy
import math
import matplotlib.dates as mdates
from pathlib import Path

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()


def distance(lat1, lon1, lat2, lon2):
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000
	return dist_km

def dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes):
	f = False
	for s in range(len(flight_latitudes)):
		for l in range(len(seismo_latitudes)):
			dist = distance(seismo_latitudes[l], seismo_longitudes[l], flight_latitudes[s], flight_longitudes[s])
			if dist <= 2:
				f = True
				break
			else:
				continue
	return f

flight_files=[]
filenames = []

# Load the seismometer location data
seismo_data = pd.read_csv('nodes_stations.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
seismo_stations = seismo_data['Latitude']
sta = seismo_data['Station']


#1-9 still needed plus redo 16,17,20  
for day in range(8, 10):
		day = '0' + str(day)
		# assign directory
		directory = '/scratch/irseppi/nodal_data/flightradar24/201903'+day+'_positions'
	
		# iterate over files in directory
		for filename in os.listdir(directory):
			filenames.append(filename)
			f = os.path.join(directory, filename)
			
			# checking if it is a file
			if os.path.isfile(f):
				flight_files.append(f)
		
for i, flight_file in enumerate(flight_files):
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']

	fname = filenames[i]	
	flight_num = fname[9:18]

	con = dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes)
	if con == True:	
		# Create a scatter plot for the seismometer locations
		for sd in range(len(seismo_data)):	
			for fd in range(len(flight_data)):
				dist = distance(seismo_latitudes[sd], seismo_longitudes[sd], flight_latitudes[fd], flight_longitudes[fd])
				if dist <= 2:
						station = str(sta[sd])
						ht = datetime.datetime.utcfromtimestamp(time[fd])
						
						h = ht.hour
						day = ht.day
						mins = ht.minute
						secs = ht.second
					
						day1 = '0'+str(day)
							
						if h < 23:			
							day2 = day1
							if h < 10:
								h_u = '0'+str(h+1)
								h = '0'+str(h)
							else:
								h_u = str(h+1)
								h = str(h)
							
						else:
							h_u = '00'
							h = '23'
							day2 = '0'+str(day + 1)
							if day1 == '09':
								day2 = '10'
						
						n = "/scratch/naalexeev/NODAL/2019-03-"+day1+"T"+h+":00:00.000000Z.2019-03-"+day2+"T"+h_u+":00:00.000000Z."+station+".mseed"
						
						if os.path.isfile(n):
							print(n)
							tr = obspy.read(n)
							try:
								tim = 120 
								tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs +tim)
								
								t                  = tr[2].times()
								data               = tr[2].data
								sampling_frequency = tr[2].stats.sampling_rate
								tr_filt = tr[2].copy()
								tr_filt.filter('highpass', freq=45, corners=2)

								title    = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'
								
								fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))     

								ax1.plot(t, tr_filt, 'k', linewidth=0.5)
								ax1.set_title(title)
								ax1.axvline(x=tim, c = 'r', ls = '--')
							
								s,f,tm,im = ax2.specgram(tr[2], Fs=sampling_frequency, noverlap=int(0.8*256), cmap='hsv')
								ax2.set_xlabel('Time - Seconds')
								ax2.axvline(x=tim, c = 'k', ls = '--')
								
								ax2.set_ylabel('Frequency (Hz)')

								ax3 = fig.add_axes([0.9, 0.1, 0.03, 0.37])

								plt.colorbar(mappable=im, cax=ax3)
								plt.ylabel('Relative Amplitude (dB)')
								BASE_DIR = "/scratch/irseppi/nodal_data/Plane_info/Plane_spec/2019-03-"+day1+"/"+flight_num+ '/'+station
								make_base_dir(BASE_DIR)
								fig.savefig('/scratch/irseppi/nodal_data/Plane_info/Plane_spec/2019-03-'+day1+'/'+flight_num + '/'+station+'/'+str(time[fd])+'_'+flight_num+'.png')
						
							except:
								continue		

	print((i/len(flight_files))*100, '% Done')		

