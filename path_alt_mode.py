import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj

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
    dist_lim = np.Infinity

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

    min_distance = np.Infinity
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

file = open('output4.txt', 'r')
file2 = pd.read_csv('input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

# Create a dictionary to store the color for each tail number
color_dict = {}
all_med = {}
points = {}
path ={}
flights = []

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
            # Assign a color to the tail number if it doesn't already have one
            if flight_num not in color_dict:
                color_dict[flight_num] = np.random.rand(3,)
                path[flight_num] = []
                all_med[flight_num] = []
                points[flight_num] = []
            
            all_med[flight_num].extend([np.nanmedian(f1)])
            points[flight_num].extend([closest_p])
            path[flight_num].extend(flight_path)
            flights.append(flight_num)
#for flight_num, med in all_med.items():
flight_num2 = 0
plt.figure(figsize=(10, 6))
plt.scatter(seismo_utm_x_km, seismo_utm_y_km, c='r', marker='x')
gg = 10
for flight_num in flights:
    gg += 13
    if flight_num2 == flight_num:
        continue
    color = color_dict[flight_num]
    p = np.array(path[flight_num])
    point = np.array(points[flight_num])
    med = all_med[flight_num]
    plt.plot(p[:,0], p[:,1], c='k')  
    for i in range(1, len(p)-1,gg):
        direction = np.arctan2(p[i+1,0] - p[i,0], p[i+1,1] - p[i,1])
        m = (p[i+1,1] - p[i,1])/(p[i+1,0] - p[i,0])
        b = p[i,0] - m*p[i,0]
        plt.quiver((p[i,0]-b)/m, p[i,1], np.cos(direction), np.sin(direction), angles='xy', color='k', scale=150)
    c = plt.scatter(point[:,0], point[:,1], c=med, zorder=10, cmap='seismic', vmin=18, vmax=22, s=100)
    #plt.scatter(seismo_utm_x_km, seismo_utm_y_km, c='r', marker='x')
    #plt.title(str(flight_num))
    #plt.colorbar(c, label= '\u0394'+'F')
    flight_num2 = flight_num
plt.colorbar(c, label= '\u0394'+'F')
plt.show()
