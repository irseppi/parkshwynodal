import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, windows
import obspy
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import os
import obspy
from obspy.core import UTCDateTime
import datetime
from numpy.fft import fft, ifft
from pathlib import Path

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()

text = open('input/all_station_crossing_db.txt', 'r')

for line in text.readlines():
	val = line.split(',')
	time = val[2]
	ht = datetime.datetime.utcfromtimestamp(int(val[2]))
	mins = ht.minute
	secs = ht.second
	h = ht.hour
	month1 = ht.month
	day1 = ht.day
	station = str(val[5])
	tim = 120
	flight_num = str(val[1])
	if ht.month == 2 and ht.day >= 11 or ht.month == 3:
		day_of_year = str((ht - datetime.datetime(2019, 1, 1)).days + 1)
	
		if val[5].isdigit() == False:
			try:
				n = "/aec/wf/2019/0"+day_of_year+"/"+str(val[5])+".*Z.20190"+day_of_year+"000000+"
				tr = obspy.read(n)
				
				tr[0].trim(tr[0].stats.starttime +(h *60 *60) + (mins * 60) + secs - tim, tr[0].stats.starttime +(h *60 *60) + (mins * 60) + secs + tim)
				data = tr[0][0:-1]
				fs = int(tr[0].stats.sampling_rate)
				title    = f'{tr[0].stats.network}.{tr[0].stats.station}.{tr[0].stats.location}.{tr[0].stats.channel} − starting {tr[0].stats["starttime"]}'						
				t                  = tr[0].times()
			except:
				continue
			
		else:	
			month2 = str(month1)
			if day1 < 10:
				day1 = '0'+str(day1)	
			if h < 23:			
				day2 = str(day1)
				if h < 10:
					h_u = '0'+str(h+1)
					h = '0'+str(h)
				else:
					h_u = str(h+1)
					h = str(h)
				
			else:
				h_u = '00'
				if month1 == '2' and day1 == '28':
					month2 = '03'
					day2 = '01'
				else:	
					if int(day1) >= 10:
						day2 = str(int(day1) + 1)
					else:
						day2 = '0'+str(int(day1)+1)
			try:		
				n = "/scratch/naalexeev/NODAL/2019-0"+str(month1)+"-"+str(day1)+"T"+str(h)+":00:00.000000Z.2019-0"+month2+"-"+day2+"T"+h_u+":00:00.000000Z."+station+".mseed"
				tr = obspy.read(n)
			
				tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
				data = tr[2][0:-1]
				fs = int(tr[2].stats.sampling_rate)
				title    = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'						
				t                  = tr[2].times()
			except:
				continue
	
		try:

			# Time array
			t = np.arange(len(data)) / fs
			g = fs*240

			# Compute spectrogram
			frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9) 

			fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))     

			ax1.plot(t, data, 'k', linewidth=0.5)
			ax1.set_title(title)
			ax1.axvline(x=tim, c = 'r', ls = '--')
			ax1.margins(x=0)
			vmin = np.min(10 * np.log10(Sxx))
			vmax = np.max(10 * np.log10(Sxx))

			# Plot spectrogram
			cax = ax2.pcolormesh(times, frequencies, 10 * np.log10(Sxx), shading='gouraud', cmap='hsv', vmin=vmin, vmax=vmax)
			ax2.set_xlabel('Time [s]')

			# Find the center of the trace
			center_index = len(data) // 2
			center_time = t[center_index]

			ax2.axvline(x=center_time, c = 'r', ls = '--')
			ax2.set_ylabel('Frequency (Hz)')
			ax2.margins(x=0)
			ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])

			plt.colorbar(mappable=cax, cax=ax3)

			ax3.set_ylabel('Relative Amplitude (dB)')
			BASE_DIR = "/scratch/irseppi/nodal_data/plane_info/plane_spec/2019-0"+str(month1)+"-"+str(day1)+"/"+flight_num+ '/'+station
			make_base_dir(BASE_DIR)
			fig.savefig('/scratch/irseppi/nodal_data/plane_info/plane_spec/2019-0'+str(month1)+'-'+str(day1)+'/'+flight_num + '/'+station+'/'+str(time)+'_'+flight_num+'.png')
			plt.close()

			# Compute and plot amplitude spectrum for the center of the trace
			# Find the index of the middle frequency
			middle_index = len(times) // 2

			# Extract the middle line of the spectrogram
			middle_column = Sxx[:, middle_index]

			# Plot the middle line of the spectrogram
			plt.figure(figsize=(10,6))
			plt.grid()
			plt.plot(frequencies, 10 * np.log10(middle_column))
			plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))

			plt.xlim(2,int(fs/2))

			plt.xlabel('Freq [Hz]')
			plt.ylabel('Amplitude [dB]')
			plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))
			plt.close()

			BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/spec/'+str(val[0])+'/'+str(val[1])+'/'+str(val[5])+'/'
			make_base_dir(BASE_DIR)
			plt.savefig('/scratch/irseppi/nodal_data/plane_info/spec/'+str(val[0])+'/'+str(val[1])+'/'+str(val[5])+'/'+str(val[5])+'_' + str(val[2]) + '.png')
			plt.close()

		except:
			continue	







			'''
			# Time array
			t = np.arange(len(data)) / fs
			g = fs*240
			# Compute spectrogram
			frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9) #, window=windows.boxcar(fs,sym=True))
			
			fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))     

			ax1.plot(t, data, 'k', linewidth=0.5)
			ax1.set_title(title)
			ax1.axvline(x=tim, c = 'r', ls = '--')
			ax1.margins(x=0)
			vmin = np.min(10 * np.log10(Sxx))
			vmax = np.max(10 * np.log10(Sxx))
			# Plot spectrogram
			cax = ax2.pcolormesh(times, frequencies, 10 * np.log10(Sxx), shading='gouraud', cmap='hsv', vmin=vmin, vmax=vmax)
			ax2.set_xlabel('Time [s]')
			
			# Find the center of the trace
			center_index = len(data) // 2
			center_time = t[center_index]
			
			ax2.axvline(x=center_time, c = 'r', ls = '--')
			ax2.set_ylabel('Frequency (Hz)')
			ax2.margins(x=0)
			ax3 = fig.add_axes([0.9, 0.11, 0.015, 0.35])
			
			plt.colorbar(mappable=cax, cax=ax3)
			ax3.set_ylabel('Relative Amplitude (dB)')
			plt.show()
			plt.close()
			
			BASE_DIR = "/scratch/irseppi/nodal_data/plane_info/plane_spec2/2019-0"+str(month1)+"-"+str(day1)+"/"+flight_num+ '/'+station
			make_base_dir(BASE_DIR)
			fig.savefig('/scratch/irseppi/nodal_data/plane_info/plane_spec2/2019-0'+str(month1)+'-'+str(day1)+'/'+flight_num + '/'+station+'/'+str(time)+'_'+flight_num+'.png')
			plt.close()
			
			plt.figure()
			# Compute and plot amplitude spectrum for the center of the trace
			window = np.hanning(fs)  # one second window
			center_data = data[int(center_index):int(center_index+fs)] * window
			freqs = np.fft.rfftfreq(fs, 1/fs)
			fft = 10*np.log10(np.fft.rfft(center_data))
			plt.plot(freqs, np.abs(fft))
			plt.title('Amplitude Spectrum at t = {:.2f} s'.format(center_time))
			plt.xlabel('Frequency [Hz]')
			plt.ylabel('Amplitude [m/s]')
			plt.xlim(2,int(fs/2))
			
			xmask = np.logical_and(freqs > .2, freqs < int(fs/2))
			plt.ylim(0,np.max(fft[xmask])*1.1)
			plt.grid()
			#plt.show()
			print(station, flight_num)
			BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/spec2/'+str(val[0])+'/'+str(val[1])+'/'+str(val[5])+'/'
			make_base_dir(BASE_DIR)
			plt.savefig('/scratch/irseppi/nodal_data/plane_info/spec2/'+str(val[0])+'/'+str(val[1])+'/'+str(val[5])+'/'+str(val[5])+'_' + str(val[2]) + '.png')
			plt.close()

			plt.figure()
			# Compute spectrogram
			frequencies, times, Sxx = spectrogram(data, fs, scaling='density', nperseg=fs, noverlap=fs * .9)

			# Find the index of the middle frequency
			middle_index = len(times) // 2

			# Extract the middle line of the spectrogram
			middle_line = Sxx[:, middle_index]

			# Plot the middle line of the spectrogram
			plt.figure(figsize=(10, 6))
			plt.plot(times, 10 * np.log10(middle_line))
			plt.xlabel('Time [s]')
			plt.ylabel('Amplitude [dB]')
			plt.title('Spectra vs Amplitude for the Middle Line of the Spectrogram')
			plt.show()
			plt.close()
			'''

