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

def wf_fft(wf):
    """ INPUT:
            wf - Numpy array of the data points in your trace
        OUTPUT:
            fft_amp - Numpy array of spectral amplitudes
            fft_freq - Numpy array of frequencies"""
    fNyq =250  
    NFFT = int(2**(np.ceil(np.log(len(wf))/np.log(2))))  # Next highest power of 2
    FFTX = np.fft.fft(wf,n=NFFT)                       # Take fft, padding with zeros.
    NumUniquePts = int(np.ceil((NFFT+1)/2))
    FFTX = FFTX[0:NumUniquePts]              # throw out neg frequencies
    MX = abs(FFTX)                           # Take magnitude of X
    MX = MX*2                                # Multiply by 2 
    fft_amp = MX/len(wf) 
    
    f = (np.arange(NumUniquePts))*2/NFFT            
    fft_freq = f*fNyq
    
    return fft_amp, fft_freq

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()
month = []
day = []

for d in range(11,29):
	day.append(str(d))
	month.append('02')

for d in range(1, 10):
	day.append('0' + str(d))
	month.append('03')

for d in range(10, 27):
	day.append(str(d))
	month.append('03')


for i in range(len(day)):
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/plane_spec/2019-'+month[i]+'-'+day[i]+'/'
	for flight in os.listdir(spec_dir):
		f = os.path.join(spec_dir, flight)
		for station in os.listdir(f):
			sta = os.path.join(f, station)
			for name in os.listdir(sta):
				time = name[0:10]
				
				ht = datetime.datetime.utcfromtimestamp(int(time))
		
				h = ht.hour
				month1 = month[i]
				day1 = day[i]
				mins = ht.minute
				secs = ht.second
				
				month2 = str(month1)
					
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
					if month == '02' and day == '28':
						month2 = '03'
						day2 = '01'
					else:	
						if int(day1) >= 10:
							day2 = str(int(day1) + 1)
						else:
							day2 = '0'+str(int(day1)+1)
						
				n = "/scratch/naalexeev/NODAL/2019-"+str(month1)+"-"+str(day1)+"T"+str(h)+":00:00.000000Z.2019-"+month2+"-"+day2+"T"+h_u+":00:00.000000Z."+station+".mseed"
				
				tr = obspy.read(n)
			
				tim = 0.5
				tr[2].trim(tr[2].stats.starttime + (mins * 60) + secs - tim, tr[2].stats.starttime + (mins * 60) + secs + tim)
		
				sampling_frequency = tr[2].stats.sampling_rate
				title    = f'{tr[2].stats.network}.{tr[2].stats.station}.{tr[2].stats.location}.{tr[2].stats.channel} − starting {tr[2].stats["starttime"]}'
				'''
				X =  fft(tr[2]) 
				N = len(X)  #Number of sample points
				n = np.arange(N) #array of number of points 0 to N
				T = N/sampling_frequency #Spacing of Samples
				freq = n/T 
				fft_amp = np.log10(np.abs(X))

				plt.figure(figsize = (9, 6))

				plt.stem(freq, np.log10(np.abs(X)), 'b', \
					 markerfmt=" ", basefmt="-b")
				plt.xlabel('Freq (Hz)')
				plt.ylabel('FFT Amplitude |X(freq)|')
				plt.xlim(1, 250)
				xmask = np.logical_and(freq > 0.2, freq < 250)
				plt.ylim(0,np.max(fft_amp[xmask]*1.1))
				#plt.title(title)

				BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/spec/2019-'+month[i]+'-'+day[i]+'/'+flight+'/'+station+'/'
				make_base_dir(BASE_DIR)
				plt.savefig('/scratch/irseppi/nodal_data/plane_info/spec/2019-'+month[i]+'-'+day[i]+'/'+flight+'/'+station+'/'+flight+'_' + str(time) + '.png')
				plt.close()
				'''
				# computing and creating a list of the amplitude spectra of the seismograms for the selected station location
				tr[2].detrend('constant')
				tr[2].detrend('linear')
				tr[2].taper(max_percentage=0.2, type="cosine")

				fft_amp, fft_freq = wf_fft(tr[2].data)

				plt.figure()
				plt.plot(fft_freq,np.log10(fft_amp))
				plt.xlim(0.2,250)

				xmask = np.logical_and(fft_freq > .2, fft_freq < 250)
				plt.ylim(0,np.max(np.log10(fft_amp[xmask])*1.1))
				#plt.title(title)
				plt.xlabel(f'Frequency (Hz)')
				plt.ylabel(f'Amplitude (m/s)')
				plt.savefig('/scratch/irseppi/nodal_data/plane_info/spec/2019-'+month[i]+'-'+day[i]+'/'+flight+'/'+station+'/'+flight+'_' + str(time) + '.png')
				plt.close()
