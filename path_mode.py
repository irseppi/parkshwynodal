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
flight_alt = {}
UTM_km_x = {}
UTM_km_y = {}
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = int(lines[1])
    nodes = int(lines[2])
    flight_file = '/scratch/irseppi/nodal_data/flightradar24/'+str(lines[0]) + '_positions/' + str(lines[0]) + '_' + str(flight_num) + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    flight_latitudes = flight_data['latitude'] 
    flight_longitudes = flight_data['longitude']
    alt = flight_data['altitude'] * 0.3048  # Convert altitude from feet to meters

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
        flight_alt[flight_num] = []
        UTM_km_x[flight_num] = []
        UTM_km_y[flight_num] = []
        flight_lat[flight_num].extend(flight_latitudes)
        flight_lon[flight_num].extend(flight_longitudes)
        UTM_km_x[flight_num].extend(flight_utm_x_km)
        UTM_km_y[flight_num].extend(flight_utm_y_km)
        flight_alt[flight_num].extend(alt)
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
    alt_t = flight_alt[flight_num]
    lat = np.array(points_lat[flight_num])
    lon = np.array(points_lon[flight_num])
    med = np.array(all_med[flight_num])
    fxx = np.array(UTM_km_x[flight_num])
    fyy = np.array(UTM_km_y[flight_num])
    #Create array of distance along path 
    dist_p = np.zeros(len(fxx))
    dist_hold = 0
    for i in range(len(fxx)-2):
        dist_p[i] = np.sqrt((np.array(fxx)[i+1] - np.array(fxx)[i]) ** 2 + (np.array(fyy)[i + 1] - np.array(fyy)[i]) ** 2) + dist_hold
        dist_hold = dist_p[i]
    
    fig = pygmt.Figure()
    with pygmt.config(MAP_DEGREE_SYMBOL= "none"):
        with fig.subplot(
            nrows=1,
            ncols=2,
            figsize=("16c", "25c"),
            margins=["0.1c", "0.1c"],autolabel=False,
        ): 
            with fig.set_panel(panel=[0,1]):
                grid = pygmt.datasets.load_earth_relief(resolution="15s", region=[-151.2, -150.05, 62.29, 63.15], registration="pixel")
                with pygmt.config(MAP_FRAME_TYPE='plain', FORMAT_GEO_MAP="ddd.xx"):
                    proj = "M15c"
                    fig.grdimage(grid=grid, projection=proj, frame="a", cmap="geo")
                    fig.colorbar(frame=["a1000", "x+lElevation (m)"], position="JMR+o8c/6c+w11.5c/0.5c")
                    fig.plot(x=np.array(f_lon), y=np.array(f_lat), pen="1p,black", projection=proj)

                    for i in range(len(f_lat) - 1):
                        if i == 0 or i == len(f_lat) - 2:
                            angle = np.arctan2(np.array(f_lat)[i + 1] - np.array(f_lat)[i], np.array(f_lon)[i + 1] - np.array(f_lon)[i])
                            angle = np.degrees(angle)
                            no = np.sqrt((np.array(f_lon)[i + 1] - np.array(f_lon)[i]) ** 2 + (np.array(f_lat)[i + 1] - np.array(f_lat)[i]) ** 2)
                            fig.plot(x=[np.array(f_lon)[i]], y=[np.array(f_lat)[i]], style="v0.1c+e", direction=[[angle], [0.2]], fill='black', pen="0.5p,black", projection=proj)

                    fig.plot(x=seismo_longitudes, y=seismo_latitudes, style="x0.2c", pen="01p,black", projection=proj)
                    fig.plot(x=-150.1072713049972, y=62.30091781635389, style="x0.3c", pen="02p,pink", projection=proj)
                    pygmt.makecpt(cmap="gmt/seis", series=[np.min(med)-0.01, np.max(med)+0.01])
                    yy = fig.plot(x=lon, y=lat, style="c0.3c", fill=med, pen="black", cmap=True, projection=proj)
                    fig.colorbar(frame=["a0.5f0.1", 'xaf+l\u0394'+'F (Hz)'], position="JMR+o8c/-6.5c+w11.5c/0.5c")

                    zoom_region = [np.min(lon) - 0.01, np.max(lon) + 0.01, np.min(lat) - 0.01, np.max(lat) + 0.01]
                    rectangle = [[zoom_region[0], zoom_region[2], zoom_region[1], zoom_region[3]]]
                    fig.plot(data=rectangle, style="r+s", pen="0.5p,black", projection=proj)

            with fig.set_panel(panel=[0,0]):
                zoom_region = [np.min(lon) - 0.01, np.max(lon) + 0.01, np.min(lat) - 0.01, np.max(lat)+ 0.01]
                with pygmt.config(MAP_FRAME_TYPE = 'plain',FORMAT_GEO_MAP="ddd.xx"):
                    grid_inset = pygmt.datasets.load_earth_relief(resolution="15s", region=zoom_region, registration="pixel")
                    cent_lon =  str(((np.max(lon) + 0.01) - np.min(lon) - 0.01)/2 + (np.min(lon) - 0.01))
                    cent_lat = str(((np.max(lat) + 0.01) - np.min(lat) - 0.01)/2 + (np.min(lat) - 0.01))
                    proj = "M6.8c"
                    cmap_limits = [float(np.min(grid)), float(np.max(grid))]  # Get min and max elevation values
                    pygmt.makecpt(cmap="geo", series=cmap_limits, continuous=True)
                
                    fig.grdimage(grid=grid_inset, region=zoom_region, projection=proj,frame="a", cmap=True)
                    fig.plot(x=np.array(f_lon), y=np.array(f_lat), projection=proj, pen="1p,black") 

                    fig.plot(x=seismo_longitudes, y=seismo_latitudes, projection=proj, style="x0.2c", pen="01p,black")

                    pygmt.makecpt(cmap="gmt/seis", series=[np.min(med)-0.01, np.max(med)+0.01]) 
                    yy = fig.plot(x=lon, y=lat, style="c0.3c", fill=med, projection=proj, pen="black", cmap=True) 
        fig.shift_origin(yshift="-8c")
        proj = "X24.5c/6c"
        with fig.subplot(nrows=1, ncols=1, figsize=("27c", "10c"), margins=["0.1c", "0.1c"],autolabel=False):
            fig.basemap(
                region=[0, np.max(dist_p), 0,np.max(alt_t)+100],  # x_min, x_max, y_min, y_max
                projection=proj,
                frame=["WSrt", "xa20+lDistance / m", "ya1000+lElevation / m"], 
            )

            points = pd.DataFrame({'longitude': [], 'latitude': []})
            for i in range(len(f_lon[:-2]) - 1):
                lon_segment = np.linspace(f_lon[i], f_lon[i + 1], 10)  # 10 intermediate points + start and end
                lat_segment = np.linspace(f_lat[i], f_lat[i + 1], 10)
                points = pd.concat([points, pd.DataFrame({'longitude': lon_segment, 'latitude': lat_segment})], ignore_index=True)
            
            elevation_data = pygmt.grdtrack(
                grid=grid,
                points=points,
                newcolname="elevation",
            )      
            ev = elevation_data['elevation'].values   

            fig.plot(
                x=[0, np.max(dist_p), np.max(dist_p), 0],
                y=[0, 0, np.max(alt_t)+100, np.max(alt_t)+100],
                fill="lightblue",
                projection=proj,
                close=True
            )

            # Interpolate 10 evenly spaced points between each pair of dist_p and ev
            interpolated_dist_p = []
            for i in range(len(dist_p[:-2]) - 1):
                dist_segment = np.linspace(dist_p[i], dist_p[i + 1], 10)  # 10 intermediate points + start and end
                interpolated_dist_p.extend(dist_segment)

            #dist_x, elv_y = np.meshgrid(np.arange(len(interpolated_dist_p)), np.arange(len(ev)))
            distance_grid, elevation_grid = np.meshgrid(
            np.linspace(0, np.max(dist_p), len(interpolated_dist_p)),
            np.linspace(0, np.max(ev) + 100, len(interpolated_dist_p))
            )
            cmap_limits = [float(np.min(grid)), float(np.max(grid))]  # Get min and max elevation values
            pygmt.makecpt(cmap="geo", series=cmap_limits, continuous=True)
            # Fill the mesh grid with elevation values for color mapping
            color_fill = np.full_like(distance_grid, 0)  # Initialize with NaN

            for j,value_2 in enumerate(elevation_grid):
                for i,value_1 in enumerate(distance_grid):
                    if value_2[i] <= float(ev[i]):
                        color_fill[i,j] = value_2[i]
                    else:
                        color_fill[i, j] = -1000  # Leave as NaN for areas outside interpolated points
            print(color_fill)
            #trurn grids into 1D arrays
            distance_grid = distance_grid.flatten()
            elevation_grid = elevation_grid.flatten()
            color_fill = color_fill.flatten()

            c_fill = pygmt.xyz2grd(x=distance_grid, y=elevation_grid, z = color_fill, projection=proj, region=[0, np.max(dist_p), 0,np.max(alt_t)+100],spacing = (np.max(dist_p)/len(interpolated_dist_p),(np.max(alt_t)+100)/len(ev)))
            print(c_fill)
            cmap_limits = [float(np.min(grid)), float(np.max(grid))]  # Get min and max elevation values
            pygmt.makecpt(cmap="geo", series=cmap_limits, continuous=True)
            fig.grdimage(grid=c_fill, projection=proj, cmap=True)

            fig.plot(x=np.array(dist_p[:-2]), y=np.array(alt_t[:-2]), pen="1p,black", projection=proj)

            fig.plot(
                x=np.array(interpolated_dist_p),
                y=np.array(ev),
                pen="1p,black",
                projection=proj,
            )

            fig.plot(
                x=np.array(distance_grid),
                y=np.array(elevation_grid),
                fill = color_fill,
                projection=proj,
                cmap=True,
            )


    fig.savefig("output.png")
    fig.show(verbose="i") 
        
    break
