from obspy import read, UTCDateTime
from obspy.geodetics import gps2dist_azimuth
import matplotlib.pyplot as plt
from prelude import *
from tqdm import trange, tqdm
import numpy as np

​event_time = UTCDateTime(2019,3,2,1,42,19)
stations = 'K24K,HDA,NEA2,CHUM,I53H?'
event_lat=64.01
event_lon = -148.76
plt.figure(figsize=(10,20))
plt.ylabel("Distance (km)")
plt.xlabel("Time (Hours)")


st = get_streams(UTCDateTime(2019,3,2,1,42,19),event_lat,event_lon,buffer_time = 2*60.0*60.0,max_dist = 13000.0)

for station in len(st):
    dist_meters = gps2dist_azimuth(st["lat"],st["lon"],event_lat,event_lon)[0]
    dist_km = dist_meters/1000.0
    tr = load_day_traces(event_time,str(station["station"].code))
    tr[0].decimate(10)
    tr[0].normalize()
    data=tr[0].data
    x = np.linspace(0,24.0,num=data.shape[0])
    plt.plot(x,3.0*data+dist_km)
plt.show()
