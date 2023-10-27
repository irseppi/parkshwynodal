#%matplotlib widget

import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
import os
import warnings

from obspy import read
from obspy.core import UTCDateTime
#from pysep import Pysep
#from pysep.recsec import plotw_rs
# script settings

warnings.filterwarnings("ignore")
plt.rcParams['figure.figsize'] = 12, 8
# set user parameters

# providing an explicit list of networks
networks      = 'AK,TA'


download  = dict( origin_time                  = UTCDateTime("2019,3,2,1,42,19"),
                  event_latitude               = 64.01,
                  event_longitude              = -148.76,
                  event_depth_km               = 5.51,
                  event_magnitude              = 1.62,
                  networks                     = networks, 
                  channels                     = 'BDF',                            # broadband channels
                  remove_response              = False,
                  remove_clipped               = False,
                  remove_insufficient_length   = False,
                  maxdistance                  = 200,
                  seconds_before_ref           = 50,
                  seconds_after_ref            = 200 / 0.3,                                # air wave travel time
                  log_level                    = 'CRITICAL',
                  write_files                  = 'sac',
                  plot_files                   = 'map',
                  output_dir                   = f'datawf/Example_{example_index}',
                  overwrite_event_tag          = f'Example_{example_index}' )

plot_info = dict( pysep_path                   = f'datawf/Example_{example_index}/Example_{example_index}/SAC',
                  sort_by                      = 'distance',
                  scale_by                     = 'normalize',
                  min_period_s                 = 0.2,
                  max_period_s                 = 1,
                  preprocess                   = 'st',
                  max_traces_per_rs            = None,
                  tmarks                       = [0],
                  save                         = '',
                  log_level                    = 'CRITICAL' )


# plot source station map

plt.figure()
source_station_map = img.imread(f'datawf/Example_{example_index}/Example_{example_index}/station_map.png')
plt.imshow(source_station_map);

# plot the record section using Pyseps's record section plotting tool
# the plotting below is only for example 1 for different seismogram sorting options

    
sort_by_tag = ['distance', 'absolute distance', 'azimuth', 'absolute azimuth', 'distance in reverse order',
               'absolute distance in reversed order', 'azimuth in reverse order', 'absolute azimuth in reverse order']

for i, sort_by in enumerate(['distance', 'abs_distance', 'azimuth', 'abs_azimuth', 'distance_r', 'abs_distance_r', 'azimuth_r', 'abs_azimuth_r']):

    plot_info = dict( pysep_path                   = f'datawf/Example_{example_index}/Example_{example_index}/SAC',
                      sort_by                      = sort_by,
                      scale_by                     = 'normalize',
                      time_shift_s                 = None,
                      min_period_s                 = 0.1,
                      max_period_s                 = 2,
                      preprocess                   = 'st',
                      max_traces_per_rs            = None,
                      azimuth_start_deg            = 0,
                      tmarks                       = [0],
                      save                         = '',
                      log_level                    = 'CRITICAL')

    print(f'\n\nCase {i+1}: Sort seismograms by {sort_by_tag[i]}\n\n')

    plotw_rs(**plot_info)


# plot the record section using Pyseps's record section plotting tool
# the plotting below is only for example 1 for different seismogram sorting options


sort_by_tag = ['None', 'P_from_TauP']

for i, sort_by in enumerate([None, 'p_arrival_time']):

    plot_info = dict( pysep_path                   = f'datawf/Example_{example_index}/Example_{example_index}/SAC',
                  sort_by                      = 'distance',
                  scale_by                     = 'normalize',
                  time_shift_s                 = sort_by,
                  min_period_s                 = 0.1,
                  max_period_s                 = 2,
                  preprocess                   = 'st',
                  max_traces_per_rs            = None,
                  azimuth_start_deg            = 0,
                  tmarks                       = [0],
                  save                         = '',
                  log_level                    = 'CRITICAL')

print(f'\n\nCase {i+1}: Align seismograms by {sort_by_tag[i]}\n\n')

plotw_rs(**plot_info)

# plot the record section using Pyseps's record section plotting tool

plotw_rs(**plot_info)


network  = '*'                                                                                 # select network
station  = '*'                                                                                 # select station
location = '*'                                                                                 # select location
channel  = '*'                                                                                 # select channel

sac_file = f'./datawf/Example_{example_index}/Example_{example_index}/SAC/Example_{example_index}.{network}.{station}.{location}.{channel}.sac'
st       = read(sac_file, 'SAC')

t                  = st[0].times()
data               = st[0].data
sampling_frequency = st[0].stats.sampling_rate

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,6))     
ax1.plot(t, data, 'k', linewidth=0.5)

image = ax2.specgram(st[0], Fs=sampling_frequency, noverlap=int(0.8*256), cmap="jet")
ax2.set_xlabel('Time - Seconds')
ax2.set_ylabel('Frequency (Hz)')

ax3 = fig.add_axes([0.9, 0.1, 0.03, 0.37])
plt.colorbar(mappable=image[3], cax=ax3)
plt.ylabel('Relative Amplitude (dB)')
plt.show()

title    = f'{st[0].stats.network}.{st[0].stats.station}.{st[0].stats.location}.{st[0].stats.channel} − starting {st[0].stats["starttime"]}'
fig.suptitle(title)

# note: the time axis of the spectrogram may not correspond to the time axis of the record sections
# this is because the spectrogram plotter assigns a default value of time = 0 to the first sample of the input data

fig.canvas.draw()
labels = np.arange(-40,100,20)
ax2.set_xticklabels(labels)
