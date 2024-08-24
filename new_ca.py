import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyproj

from prelude import load_flights, dist_less

# Load flight files and filenames
flight_files,filenames = load_flights(2, 4, 11, 27)

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']

# Convert latitude and longitude to UTM coordinates
utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')
seismo_utm = [utm_proj(lon, lat) for lat, lon in zip(seismo_latitudes, seismo_longitudes)]
seismo_utm_x, seismo_utm_y = zip(*seismo_utm)

# Convert UTM coordinates to kilometers
seismo_utm_x_km = [x / 1000 for x in seismo_utm_x]
seismo_utm_y_km = [y / 1000 for y in seismo_utm_y]
seismo_utm_km = [(x, y) for x, y in zip(seismo_utm_x_km, seismo_utm_y_km)]

def closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, seismo_utm_x, seismo_utm_y):
    closest_point = None
    dist_lim = 2.01
    x = [flight_utm_x1, flight_utm_x2]
    y = [flight_utm_y1, flight_utm_y2]  
    m = (y[1]-y[0])/(x[1]-x[0])
    b = y[0] - m*x[0]
    if m >= 0:
        ggg = 10
    else:
        ggg = -10
    for point in np.arange(x[0]+ggg, x[1]-ggg, ggg):
        xx = point
        yy = m*xx + b
        dist_km = np.sqrt((seismo_utm_y-yy)**2 +(seismo_utm_x-xx)**2) 
        print(dist_km)
        if dist_km < dist_lim:
            dist_lim = dist_km
            closest_point = (xx,yy)
        else:
            continue

    return closest_point, dist_lim


def find_closest_point(flight_path, seismo_utm):
    min_distance = 2.01
    closest_point = None
    for i in range(len(flight_utm) - 1):
        flight_utm_x1, flight_utm_y1 = flight_path[i]
        flight_utm_x2, flight_utm_y2 = flight_path[i + 1]
        seismo_utm_x, seismo_utm_y = seismo_utm
        point, d = closest_point_on_segment(flight_utm_x1, flight_utm_y1, flight_utm_x2, flight_utm_y2, seismo_utm_x, seismo_utm_y)
        if point == None:
            continue
        elif d < min_distance:
            min_distance = d
            closest_point = point
        else:
            continue

    return closest_point, min_distance

# Iterate over flight files
for i, flight_file in enumerate(flight_files):
    print((i/len(flight_files))*100, '%')	# Print progress

    # Load flight data and extract relevant columns
    flight_data = pd.read_csv(flight_file, sep=",") 
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']

    if dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes) == False:
        continue
    else:
        # Convert flight latitude and longitude to UTM coordinates
        flight_utm = [utm_proj(lon, lat) for lat, lon in zip(flight_latitudes, flight_longitudes)]
        flight_utm_x, flight_utm_y = zip(*flight_utm)

        # Convert UTM coordinates to kilometers
        flight_utm_x_km = [x / 1000 for x in flight_utm_x]
        flight_utm_y_km = [y / 1000 for y in flight_utm_y]
        flight_path = [(x,y) for x, y in zip(flight_utm_x_km, flight_utm_y_km)]
        
        # Iterate over seismometer data
        for s in range(len(seismo_data)):
            seismometer = (seismo_utm_x_km[s], seismo_utm_y_km[s])  
            
            closest_p, d = find_closest_point(flight_path, seismometer)
            if closest_p == None:
                continue    
            # Extracting latitudes and longitudes for plotting

            closest_x, closest_y = closest_p

            y =[closest_y,  seismometer[1]]
            x = [closest_x, seismometer[0]]
            yy = sum(y)/len(y)
            xx = sum(x)/len(x)

            plt.figure()
            plt.scatter(flight_utm_x_km, flight_utm_y_km)
            plt.scatter(seismometer[0], seismometer[1], 'ro')
            plt.scatter(closest_x, closest_y, 'go')
            plt.plot(x, y, 'b--')
            plt.xlim(xx-10, xx+10)
            plt.ylim(yy-10, yy+10)
            plt.axis('equal')
            plt.show()


