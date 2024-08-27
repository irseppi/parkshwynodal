import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from matplotlib.patches import Rectangle
import pyproj

from prelude import load_flights, dist_less, make_base_dir

# Load flight files and filenames
flight_files,filenames = load_flights(2, 4, 11, 27)

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3
def f(y):
    return -1.0/np.cos((y)*np.pi/180)
def rf(y):
    return (-180*np.arccos(-1/y))/np.pi
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
    dist_lim = np.Infinity
    x = [flight_utm_x1, flight_utm_x2]
    y = [flight_utm_y1, flight_utm_y2]  
    m = (y[1]-y[0])/(x[1]-x[0])
    b = y[0] - m*x[0]

    if (x[1]-x[0]) <= 0:
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


def find_closest_point(flight_path, seismo_utm):
    min_distance = np.Infinity
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
            index = i
        else:
            continue

    return closest_point, min_distance, index

# Iterate over flight files
for i, flight_file in enumerate(flight_files):
    print((i/len(flight_files))*100, '%')	# Print progress

    # Load flight data and extract relevant columns
    flight_data = pd.read_csv(flight_file, sep=",") 
    flight_latitudes = flight_data['latitude']
    flight_longitudes = flight_data['longitude']
    head = flight_data['heading']
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

            closest_p, d, index= find_closest_point(flight_path, seismometer)
            # Extracting latitudes and longitudes for plotting

            closest_x, closest_y = closest_p
            y = [closest_y, seismometer[1]]
            x = [closest_x, seismometer[0]]
            yy = sum(y) / len(y)
            xx = sum(x) / len(x)

            #plt.figure()
            #plt.plot(flight_utm_x_km, flight_utm_y_km)
            #plt.scatter(flight_utm_x_km, flight_utm_y_km)
            #plt.scatter(seismometer[0], seismometer[1], c='red')
            #plt.scatter(closest_x, closest_y, c='green')
            #plt.plot(x, y, 'b--')
            #plt.axis('equal')
            #plt.xlim(xx - (x[1] - x[0]) / 1.5, xx + (x[1] - x[0]) / 1.5)
            #plt.ylim(yy - (y[1] - y[0]) / 1.5, yy + (y[1] - y[0]) / 1.5)
            #plt.title('Distance: ' + str(d) + ' km')
            #plt.grid(True, linestyle='dotted', color='gray')
            #plt.show()
            largm = 2
            if d <= 5:
                if largm == 1:
                    # Create a figure with two subplots side by side
                    fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 2]})
                    fig.subplots_adjust(wspace=0.5)  # Adjust the spacing between subplots

                    minl = (xx - (x[1] - x[0]) / 1.5)
                    maxl = (xx + (x[1] - x[0]) / 1.5)
                    minla = (yy - (y[1] - y[0]) * 2)
                    maxla = (yy + (y[1] - y[0]) * 2)
                    axs[1].set_xticks(np.arange(minl, maxl, (x[1] - x[0]) / 2))
                    axs[1].set_yticks(np.arange(minla, maxla, (y[1] - y[0]) / 2))
                    axs[0].set_xticks(np.arange(min_lon, max_lon, 1))
                    axs[0].set_yticks(np.arange(min_lat, max_lat, 1))

                    axs[0].set_xscale('function', functions=(f, rf))
                    axs[1].set_aspect('equal')
                    axs[0].grid(True, linestyle='dotted', color='gray')
                    axs[1].grid(True, linestyle='dotted', color='gray')
                    axs[0].scatter(seismo_longitudes, seismo_latitudes, c='#e41a1c', s = 3, label='seismometers')
                    axs[0].plot(flight_longitudes, flight_latitudes, '-', c='#377eb8', lw=1, ms = 1, label='flight path')
                    for i in range(1, len(flight_latitudes)-1, int(len(flight_latitudes)/5)):
                        direction = np.arctan2(flight_latitudes[i+1] - flight_latitudes[i], flight_longitudes[i+1] - flight_longitudes[i])
                        m = (flight_latitudes[i+1] - flight_latitudes[i])/(flight_longitudes[i+1] - flight_longitudes[i])
                        b = flight_latitudes[i] - m*flight_longitudes[i]
                        axs[0].quiver((flight_latitudes[i]-b)/m, flight_latitudes[i], np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', headwidth = 5)
                    # Set labels and title
                    axs[0].set_xlim(min_lon, max_lon)
                    axs[0].set_ylim(min_lat, max_lat)
                    axs[0].tick_params(axis='both', which='major', labelsize=12)
                    

                    heading = np.deg2rad(head[index])
                    # Define the UTM and latitude/longitude coordinate systems

                    minlon, minlat = utm_proj(minl*1000, minla*1000, inverse=True)
                    maxlon, maxlat = utm_proj(maxl*1000, maxla*1000, inverse=True)
                    rect = Rectangle((minlon, minlat), (maxlon-minlon), (maxlat-minlat), ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
                    axs[0].add_patch(rect)
                    axs[1].plot(x,y, '--', c='#ff7f00')

                    # Draw the zoomed in map on the second subplot
                    axs[1].plot(flight_utm_x_km, flight_utm_y_km, c='#377eb8',linestyle ='dotted')
                    axs[1].set_xlim(minl, maxl)
                    axs[1].set_ylim(minla, maxla)
                    axs[1].tick_params(axis='both', which='major', labelsize=9)
                    axs[1].ticklabel_format(useOffset=False, style='plain')
                    if len(axs[1].get_xticklabels()) > 4:
                        axs[1].set_xticklabels([round(x, 4) for x in axs[1].get_xticks()], rotation=20, fontsize=9)
                    m = (flight_latitudes[index+1] - flight_latitudes[index])/(flight_longitudes[index+1] - flight_longitudes[index])
                    b = flight_latitudes[index] - m*flight_longitudes[index]

                    direction = np.arctan2(flight_latitudes[index+1] - flight_latitudes[index], flight_longitudes[index+1] - flight_longitudes[index])
                        
                    axs[1].quiver(closest_x, closest_y, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=7)

                    axs[1].quiver(closest_x, closest_y, np.cos(heading), np.sin(heading), angles='xy', scale = 7, color='#999999')
                    axs[1].scatter(closest_x, closest_y, c='#377eb8', s=50, zorder=3)
                    axs[1].scatter(seismometer[0], seismometer[1], c='#e41a1c', s=50, zorder=3)

                    axs[1].text(xx,yy, str(round(d, 3))+' km', fontsize=15, fontweight='bold')

                    # Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
                    con = mpatch.ConnectionPatch(xyA=(minl, minla), xyB=(maxlon, minlat), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                    con = mpatch.ConnectionPatch(xyA=(minl, maxla), xyB=(maxlon, maxlat), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                    plt.show()
                    #BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/map_all/' + date + '/'+flight+'/'+station+'/'
                    #make_base_dir(BASE_DIR)
                    #plt.savefig('/scratch/irseppi/nodal_data/plane_info/map_all/'+ date + '/'+flight+'/'+station+'/map_'+flight+'_' + str(time[l]) + '.png')
                    plt.close()

                elif largm == 2:
                    # Create a figure with two subplots side by side
                    fig, axs = plt.subplots(1, 2) 
                    fig.subplots_adjust(wspace=0.5)  # Adjust the spacing between subplots

                    lxmin,lymin = utm_proj(min_lon, min_lat)
                    lxmax, lymax = utm_proj(max_lon, max_lat)

                    min_x = (xx - (x[1] - x[0]) / 1.5)
                    max_x = (xx + (x[1] - x[0]) / 1.5)
                    min_y = (yy - (y[1] - y[0]) * 4)
                    max_y = (yy + (y[1] - y[0]) * 4)

                    axs[1].set_xticks(np.arange(min_x, max_x, (x[1] - x[0]) / 2))
                    axs[1].set_yticks(np.arange(min_y, max_y, (y[1] - y[0]) / 2))

                    axs[0].set_xticks(np.arange(lxmin/1000, lxmax/1000, 50))
                    axs[0].set_yticks(np.arange(lymin/1000, lymax/1000, 50))

                    axs[0].set_aspect('equal')
                    axs[1].set_aspect('equal')
                    axs[0].grid(True, linestyle='dotted', color='gray')
                    axs[1].grid(True, linestyle='dotted', color='gray')
                    axs[0].scatter(seismo_utm_x_km, seismo_utm_y_km, c='#e41a1c', s = 3, label='seismometers')
                    axs[0].plot(flight_utm_x_km, flight_utm_y_km, '-', c='#377eb8', lw=1, ms = 1, label='flight path')
                    for i in range(1, len(flight_utm_y_km)-1, int(len(flight_utm_y_km)/5)):
                        direction = np.arctan2(flight_utm_y_km[i+1] - flight_utm_y_km[i], flight_utm_x_km[i+1] - flight_utm_x_km[i])
                        m = (flight_utm_y_km[i+1] - flight_utm_y_km[i])/(flight_utm_x_km[i+1] - flight_utm_x_km[i])
                        b = flight_utm_y_km[i] - m*flight_utm_x_km[i]
                        axs[0].quiver((flight_utm_y_km[i]-b)/m, flight_utm_y_km[i], np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', headwidth = 5)
                    # Set labels and title
                    axs[0].set_xlim(lxmin/1000, lxmax/1000)
                    axs[0].set_ylim(lymin/1000, lymax/1000)
                    axs[0].tick_params(axis='both', which='major', labelsize=12)
                    
                    heading = np.deg2rad(head[index])
                    # Define the UTM and latitude/longitude coordinate systems

                    rect = Rectangle((min_x, min_y), (max_x-min_x), (max_y-min_y), ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5)
                    axs[0].add_patch(rect)
                    axs[1].plot(x,y, '--', c='#ff7f00')

                    # Draw the zoomed in map on the second subplot
                    axs[1].plot(flight_utm_x_km, flight_utm_y_km, c='#377eb8',linestyle ='dotted')
                    axs[1].set_xlim(min_x, max_x)
                    axs[1].set_ylim(min_y, max_y)
                    axs[1].tick_params(axis='both', which='major', labelsize=9)
                    axs[1].ticklabel_format(useOffset=False, style='plain')
                    if len(axs[1].get_xticklabels()) > 4:
                        axs[1].set_xticklabels([round(x, 4) for x in axs[1].get_xticks()], rotation=20, fontsize=9)
                    m = (flight_utm_y_km[index+1] - flight_utm_y_km[index])/(flight_utm_x_km[index+1] - flight_utm_x_km[index])
                    b = flight_utm_y_km[index] - m*flight_utm_x_km[index]

                    direction = np.arctan2(flight_utm_y_km[index+1] - flight_utm_y_km[index], flight_utm_x_km[index+1] - flight_utm_x_km[index])
                        
                    axs[1].quiver(closest_x, closest_y, np.cos(direction), np.sin(direction), angles='xy', color='#377eb8', scale=7)

                    axs[1].quiver(closest_x, closest_y, np.cos(heading), np.sin(heading), angles='xy', scale = 7, color='#999999')
                    axs[1].scatter(closest_x, closest_y, c='#377eb8', s=50, zorder=3)
                    axs[1].scatter(seismometer[0], seismometer[1], c='#e41a1c', s=50, zorder=3)

                    axs[1].text(xx,yy, str(round(d, 3))+' km', fontsize=15, fontweight='bold')

                    # Draw dashed lines connecting the rectangle on the existing map to the zoomed-in map
                    con = mpatch.ConnectionPatch(xyA=(min_x, min_y), xyB=(min_x, min_y), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                    con = mpatch.ConnectionPatch(xyA=(min_x, max_y), xyB=(min_x, max_y), coordsA="data", coordsB="data", axesA=axs[1], axesB=axs[0], color="black", linestyle="--")
                    fig.add_artist(con)
                    plt.show()
                    #BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/map_all/' + date + '/'+flight+'/'+station+'/'
                    #make_base_dir(BASE_DIR)
                    #plt.savefig('/scratch/irseppi/nodal_data/plane_info/map_all/'+ date + '/'+flight+'/'+station+'/map_'+flight+'_' + str(time[l]) + '.png')
                    plt.close()
