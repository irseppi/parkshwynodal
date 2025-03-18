import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import pygmt

#Load seismometer data
seismo_data = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

# Define the UTM projection with zone 6 and WGS84 ellipsoid
utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

#Convert to UTM coordinates
seismo_utm = [utm_proj(lon, lat) for lat, lon in zip(seismo_latitudes, seismo_longitudes)]
seismo_utm_x, seismo_utm_y = zip(*seismo_utm)

# Convert UTM coordinates to kilometers
seismo_utm_x_km = [x / 1000 for x in seismo_utm_x]
seismo_utm_y_km = [y / 1000 for y in seismo_utm_y]

######################################################################################################################################

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

#######################################################################################################################################

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

##################################################################################################################################################

file = open('output3.txt', 'r')
file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

x_airport, y_airport = utm_proj(-150.1072713049972,62.30091781635389)

# Create a dictionary to store the color for each tail number
color_dict = {}
all_med = {}
points = {}
points_lat = {}
points_lon = {}
path ={}
flights = []
tail = {}
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = lines[1]
    nodes = int(lines[2])
    distance = float(lines[6])
    vo = float(lines[4]) #try with vo and also with diffrence between heading direction and travel direction
    flight_file = '/scratch/irseppi/nodal_data/flightradar24/'+str(lines[0]) + '_positions/' + str(lines[0]) + '_' + flight_num+ '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    speed = flight_data['speed']
    alt = flight_data['altitude']
    flight_num = int(flight_num)
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
            elevation = elevations[s] 
            break
        else:
            continue

    closest_p, dist_km, index = find_closest_point(flight_path, seismometer)
    closest_x, closest_y = closest_p
    closest_lon, closest_lat = utm_proj(closest_x * 1000, closest_y * 1000, inverse=True)
    # Convert altitude and speed to meters per second
    alt_m = alt[index] * 0.3048
    speed_mps = speed[index] * 0.514444

    # Calculate the height and distance in meters
    height_m = alt_m - elevation 
    dist_m = dist_km * 1000

    peaks = np.array(lines[7])

    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))

    ppp = []
    f1 = []
    peak_old = 0
    for peak in peaks:
        if peak == '' or peak == ' ' or peak == '   ':
            continue
        else:
            if np.abs(float(peak) - float(peak_old))< 10:
                continue
            ppp.append(float(peak))

            if len(peaks) == 0 or peak == peaks[0]:
                peak_old = peak
                continue
            diff = float(peak) - float(peak_old)
            if diff > 22 or diff < 18:
                continue
            f1.append(diff)
            
        peak_old = float(peak)
    if len(f1) <= 1:
        continue
    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            if tail_num != 10512184:
                continue

            if flight_num not in path:
                path[flight_num] = []
                all_med[flight_num] = []
                points[flight_num] = []
                points_lat[flight_num] = []
                points_lon[flight_num] = []
                tail[flight_num] = []
            all_med[flight_num].extend([np.nanmedian(f1)])
            #points[flight_num].extend([closest_p])
            points_lat[flight_num].extend([closest_lat])
            points_lon[flight_num].extend([closest_lon])
            path[flight_num].extend(flight_path)
            tail[flight_num].extend([tail_num])
            flights.append(flight_num)
flight_num2 = 0

for flight_num in flights:
    if flight_num2 == flight_num:
        continue
    tail_num = tail[flight_num][0]
    p = np.array(path[flight_num])

    point = np.array(points[flight_num])
    lat = np.array(points_lat[flight_num])
    lon = np.array(points_lon[flight_num])
    med = all_med[flight_num]

    fig = pygmt.Figure()
    #grid = pygmt.datasets.load_earth_relief(resolution="30s", region=[-151.5, -150, 62.2, 63.3], registration="gridline") 
    grid = pygmt.datasets.load_earth_relief(resolution="15s", region=[-151.5, -150, 62.2, 63.3], registration="pixel")
    fig.grdimage(grid=grid, projection="M15c", frame="a", cmap="geo")
    #fig.colorbar(frame=["a1000", "x+lElevation", "y+lm"])
    fig.plot(x=np.array(flight_longitudes), y=np.array(flight_latitudes),pen="02p,black") 
    num_arrows = 4
    arrow_step = len(flight_latitudes) % (num_arrows)
    print(arrow_step)
    for i in range(1, len(flight_latitudes) - 1, arrow_step):
        angle = np.arctan2(np.array(flight_latitudes)[i + 1] - np.array(flight_latitudes)[i],np.array(flight_longitudes)[i + 1] - np.array(flight_longitudes)[i])
        angle = np.degrees(angle)
        no = np.sqrt((np.array(flight_longitudes)[i + 1] - np.array(flight_longitudes)[i]) ** 2 + (np.array(flight_latitudes)[i + 1] - np.array(flight_latitudes)[i]) ** 2)

        if i <= arrow_step * 2 or i >= len(flight_latitudes) - arrow_step * 2:
            fig.plot(x=[np.array(flight_longitudes)[i]],y=[np.array(flight_latitudes)[i]],style="v0.9c+e",direction=[[angle], [no]],color='black',pen="10p,black")

    pygmt.makecpt(cmap="gmt/split", series=[18,22]) 
    fig.plot(x=seismo_longitudes, y=seismo_latitudes, style="x0.2c",pen="01p,black")

    yy = fig.plot(x=lon, y=lat, style="c0.3c",fill=med, pen="black", cmap=True) 
    fig.colorbar(frame = 'xaf+l\u0394'+'F')
    fig.show()
    break
'''
    plt.figure()
    plt.plot(p[:, 0], p[:, 1], c='k')

    # Calculate the number of arrows to display
    num_arrows = 4
    arrow_step = len(p) // (num_arrows + 1)

    for i in range(1, len(p) - 1, arrow_step):
        direction = np.arctan2(p[i + 1, 0] - p[i, 0], p[i + 1, 1] - p[i, 1])
        no = np.sqrt((p[i + 1, 0] - p[i, 0]) ** 2 + (p[i + 1, 1] - p[i, 1]) ** 2)
        if i <= arrow_step * 2 or i >= len(p) - arrow_step * 2:
            yy = plt.quiver(p[i, 0], p[i, 1], (p[i + 1, 0] - p[i, 0]) / no, (p[i + 1, 1] - p[i, 1]) / no, angles='xy',
                            color='k', pivot='tail', headwidth=5, scale=100)
    yy = plt.scatter(point[:, 0], point[:, 1], c=med, zorder=10, cmap='seismic', vmin=18, vmax=22, s=100)
    plt.title(str(flight_num))
    flight_num2 = flight_num
    plt.scatter(seismo_utm_x_km, seismo_utm_y_km, c='k', marker='x')
    #plt.scatter(x_airport/1000, y_airport/1000,c = 'pink', marker='x',zorder = 10)
    plt.colorbar(yy, label= '\u0394'+'F')
    plt.axis('equal')
    plt.show()
'''