import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import pygmt
from prelude import *

#Load seismometer data
seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']

# Define the UTM projection with zone 6 and WGS84 ellipsoid
utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

#Convert to UTM coordinates
seismo_utm = [utm_proj(lon, lat) for lat, lon in zip(seismo_latitudes, seismo_longitudes)]
seismo_utm_x, seismo_utm_y = zip(*seismo_utm)

# Convert UTM coordinates to kilometers
seismo_utm_x_km = [x / 1000 for x in seismo_utm_x]
seismo_utm_y_km = [y / 1000 for y in seismo_utm_y]

file = open('C185data_atm_1o.txt', 'r')
file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

x_airport, y_airport = utm_proj(-150.1072713049972,62.30091781635389)

# Create a dictionary to store the color for each tail number
all_med = {}
points_lat = {}
points_lon = {}
flights = []
tail = {}
flight_lat = {}
flight_lon = {}
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = int(lines[1])
    nodes = int(lines[2])
    flight_file = '/scratch/irseppi/nodal_data/flightradar24/'+str(lines[0]) + '_positions/' + str(lines[0]) + '_' + str(flight_num) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude'] 
    flight_longitudes = flight_data['longitude']

    # Convert flight latitude and longitude to UTM coordinates
    flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
    flight_utm_x, flight_utm_y = zip(*flight_utm)

    # Convert UTM coordinates to kilometers
    flight_utm_x_km = [x / 1000 for x in flight_utm_x]
    flight_utm_y_km = [y / 1000 for y in flight_utm_y]
    flight_path = [(x,y) for x, y in zip(flight_utm_x_km, flight_utm_y_km)]

    # Find find the data for the example station
    for s in range(len(seismo_data)):
        if str(nodes) == str(stations[s]):
            seismometer = (seismo_utm_x_km[s], seismo_utm_y_km[s]) 
            break
        else:
            continue

    closest_p, dist_km, index = find_closest_point(flight_path, seismometer)
    closest_x, closest_y = closest_p
    closest_lon, closest_lat = utm_proj(closest_x * 1000, closest_y * 1000, inverse=True)

    peaks = np.array(lines[7])

    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = sorted(np.array(peaks.split(' ')).astype(float))

    f1 = []
    peak_old = 0
    for peak in peaks:
        if np.abs(float(peak) - float(peak_old)) < 10:
            continue
        if peak == peaks[0]:
            peak_old = peak
            continue

        diff = float(peak) - float(peak_old)
        if diff > 22 or diff < 18:
            peak_old = peak
            continue
        f1.append(diff)
        peak_old = float(peak)
    if len(f1) <= 1:
        continue
    for lp in range(len(flight)):   
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            break
    if tail_num != 10512184:
        continue

    if flight_num not in all_med:
        all_med[flight_num] = []
        points_lat[flight_num] = []
        points_lon[flight_num] = []
        tail[flight_num] = []
        flight_lat[flight_num] = []
        flight_lon[flight_num] = []
        flight_lat[flight_num].extend(flight_latitudes)
        flight_lon[flight_num].extend(flight_longitudes)
        tail[flight_num].extend([tail_num])
    all_med[flight_num].extend([np.nanmedian(f1)])
    points_lat[flight_num].extend([closest_lat])
    points_lon[flight_num].extend([closest_lon])
    flights.append(flight_num)

flight_num2 = 0

for flight_num in flights:
    if flight_num2 == flight_num:
        continue
    flight_num2 = flight_num
    tail_num = tail[flight_num][0]
    f_lat = flight_lat[flight_num]
    f_lon = flight_lon[flight_num]
    lat = np.array(points_lat[flight_num])
    lon = np.array(points_lon[flight_num])
    med = np.array(all_med[flight_num])
    print(lat,lon)
    fig = pygmt.Figure()

    grid = pygmt.datasets.load_earth_relief(resolution="15s", region=[-151.3, -150.05, 62.29, 63.15], registration="pixel")
    pygmt.config(MAP_FRAME_TYPE = 'plain',FORMAT_GEO_MAP="ddd.x")

    fig.grdimage(grid=grid, projection="M15c",frame="WSne",cmap="geo")
    fig.basemap(frame=["WSne", "xaf+lx-axis", "yaf+ly-axis"])
    fig.colorbar(frame=["a1000", "x+lElevation (m)"], position="JMR+o0.5c/5.5c+w10c/0.5c")
    fig.plot(x=np.array(f_lon), y=np.array(f_lat),pen="1p,black") 

    for i in range(len(f_lat) - 1):
        if i == 0 or i == len(f_lat)-2:
            angle = np.arctan2(np.array(f_lat)[i + 1] - np.array(f_lat)[i],np.array(f_lon)[i + 1] - np.array(f_lon)[i])
            angle = np.degrees(angle)
            no = np.sqrt((np.array(f_lon)[i + 1] - np.array(f_lon)[i]) ** 2 + (np.array(f_lat)[i + 1] - np.array(f_lat)[i]) ** 2)
            fig.plot(x=[np.array(f_lon)[i]],y=[np.array(f_lat)[i]],style="v0.1c+e",direction=[[angle], [0.2]],fill='black',pen="0.5p,black")

    fig.plot(x=seismo_longitudes, y=seismo_latitudes, style="x0.2c",pen="01p,black")
    fig.plot(x=-150.1072713049972, y=62.30091781635389, style="x0.3c",pen="02p,pink")
    pygmt.makecpt(cmap="gmt/seis", series=[np.min(med)-0.1,np.max(med)+0.1]) 
    yy = fig.plot(x=lon, y=lat, style="c0.3c",fill=med, pen="black", cmap=True) 
    fig.colorbar(frame=["a1", 'xaf+l\u0394'+'F (Hz)'], position="JMR+o0.5c/-5.5c+w10c/0.5c")

    pygmt.config(MAP_FRAME_TYPE = 'plain',FORMAT_GEO_MAP="ddd.xxx")
    zoom_region = [np.min(lon) - 0.01, np.max(lon) + 0.01, np.min(lat) - 0.01, np.max(lat)+ 0.01]
    grid_inset = pygmt.datasets.load_earth_relief(resolution="15s", region=zoom_region, registration="pixel")

    cmap_limits = [float(np.min(grid)), float(np.max(grid))]  # Get min and max elevation values
    pygmt.makecpt(cmap="geo", series=cmap_limits, continuous=True)

    with fig.inset(position="JML+o0.5i/0i+w4i/4i", region=zoom_region):
        fig.grdimage(grid=grid_inset, region=zoom_region, projection="T15c/-145/0", frame="a", cmap=True) #M4c
        fig.basemap(frame=["WSne", "xaf+lx-axis", "yaf+ly-axis"], region=zoom_region )
        fig.plot(x=np.array(f_lon), y=np.array(f_lat), pen="1p,black") 

        fig.plot(x=seismo_longitudes, y=seismo_latitudes, style="x0.2c",pen="01p,black")

        pygmt.makecpt(cmap="gmt/seis", series=[np.min(med) - 0.1, np.max(med) + 0.1]) 
        yy = fig.plot(x=lon, y=lat, style="c0.3c",fill=med, pen="black", cmap=True) 

    fig.savefig("output.png", crop=True, dpi=300)
    fig.show(method="external") 
    break
