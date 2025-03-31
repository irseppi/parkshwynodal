import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from pyproj import Proj
from prelude import *
from scipy.optimize import minimize

# Function to calculate L1 norm and produce best fit line for plot
def fit_l1_line(x, y, bounds=None):
    """
    Fits a line to the given data points (x, y) using the L1 norm (minimizing absolute deviations).

    Parameters:
        x (list or np.array): Independent variable data points.
        y (list or np.array): Dependent variable data points.
        bounds (tuple, optional): A tuple (min_val, max_val) specifying the range of values to consider for both x and y.

    Returns:
        tuple: Slope (m) and intercept (b) of the best fit line.
    """
    # Apply bounds if provided
    if bounds is not None:
        min_val, max_val = bounds
        x = np.array(x)  # Ensure x is a NumPy array
        y = np.array(y)  # Ensure y is a NumPy array
        mask = (x >= min_val) & (x <= max_val) & (y >= min_val) & (y <= max_val)
        x = np.array(x)[mask]
        y = np.array(y)[mask]

    # Define the objective function for L1 norm
    def L1(params):
        m, b = params
        x_array = np.array(x)  # Ensure x is a NumPy array
        return np.sum(np.abs(y - (m * x_array + b)))

    # Initial guess for slope and intercept
    initial_guess = [0, 0]

    # Minimize the L1 norm
    result = minimize(L1, initial_guess)

    # Extract the slope and intercept from the result
    m, b = result.x
    return m, b

seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
stations = seismo_data['Station']
elevations = seismo_data['Elevation']

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt','r')

latc = []
lonc = []
altc = []
times_list = []
speeds_list = []
dists_list = []
stat_list = []
for line in file_in.readlines():
    text = line.split(',')
    x =  float(text[2])  
    y = float(text[3])  

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)
    latc.append(lat)
    lonc.append(lon)

    altc.append(float(text[4])*0.0003048) #convert between feet and km
    times_list.append(float(text[5]))
    speeds_list.append(float(text[7]))
    dists_list.append(float(text[4]))
    stat_list.append(text[9])
file_in.close()

file_list = ['C185data_atm_1o.txt','C185data_atm_full.txt','C185data_1o.txt','C185data_full.txt']
title = ['OH/FT','OH/VT','FH/FT','FH/VT']

fig, axs = plt.subplots(4, 3, figsize=(18, 24), sharey=False)

for idx, fil in enumerate(file_list):
    data = open(fil, 'r')

    time_new = []
    v0_new = []
    distance_new = []

    times_org = []
    speeds_org = []
    dists_org = []
    
    date = []
    y=0
    temp_c = []

    for line in data.readlines():
        y += 1
        counts = []
        lines = line.split(',')
        time = float(lines[3])
        flight_num = lines[1]
        date_lab = lines[0]
        for i in range(len(latc)):
            if times_list[i] == time:
                index_UTC = i
                lat = latc[i]
                lon = lonc[i]
                alt = altc[i]
                sta = stat_list[i]
                break

        input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat'

        try:
            file =  open(input_files, 'r') 
        except:
            continue
        data = json.load(file)
        if idx == 2 or idx == 3:
            Tc = -2
            c = speed_of_sound(Tc)
        else:
            # Extract metadata
            metadata = data['metadata']
            sourcefile = metadata['sourcefile']
            datetim = metadata['time']['datetime']
            latitude = metadata['location']['latitude']
            longitude = metadata['location']['longitude']
            parameters = metadata['parameters']

            # Extract data
            data_list = data['data']

            # Convert data to a DataFrame
            data_frame = pd.DataFrame(data_list)

            # Find the "Z" parameter and extract the value at index
            z_index = None
            hold = np.inf
            for item in data_list:
                if item['parameter'] == 'Z':
                    for i in range(len(item['values'])):
                        if abs(float(item['values'][i]) - float(alt)) < hold:
                            hold = abs(float(item['values'][i]) - float(alt))
                            z_index = i
            for item in data_list:
                if item['parameter'] == 'T':
                    Tc = - 273.15 + float(item['values'][z_index])
            c = speed_of_sound(Tc)

        flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + str(date_lab) + '_positions/' + str(date_lab) + '_' + str(flight_num) + '.csv'
        flight_data = pd.read_csv(flight_file, sep=",")
        flight_latitudes = flight_data['latitude']
        flight_longitudes = flight_data['longitude']
        time = flight_data['snapshot_id']
        timestamps = flight_data['snapshot_id']
        speed = flight_data['speed']
        altitude = flight_data['altitude']

        closest_x, closest_y, dist_km, closest_time, tarrive, alt, sp, elevation, speed_mps, height_m, dist_m, tmid = closest_approach_UTM(seismo_latitudes, seismo_longitudes, flight_latitudes, flight_longitudes, timestamps, altitude, speed, stations, elevations, c, sta)
        if closest_x == None:
            continue
        
        #To set the initial window of arrival correct picks your start end Must use the tarrive time to get the correct data
        ta_old = calc_time(tmid,dist_m,height_m,343)

        temp_c.append(Tc)

        time_new.append(float(lines[4]))
        times_org.append(tarrive - (ta_old-120))

        v0_new.append(float(lines[5]))
        distance_new.append(float(lines[6]))

        speeds_org.append(speeds_list[index_UTC])
        dists_org.append(dists_list[index_UTC])
        date.append(y)
    if idx == 0 or idx == 1:
        scatter1 = axs[idx, 0].scatter(v0_new, speeds_org, c=temp_c, cmap='coolwarm', s=15)
        axs[idx, 0].set_title(f"{title[idx]}: Velocity (m/s)", fontsize=10)
        axs[idx, 0].set_xlim(50, 80)
        axs[idx, 0].axline((0, 0), slope=1, color='black', linestyle='--')
        axs[idx, 0].set_ylim(50, 80)
        axs[idx, 0].set_aspect('equal')
        axs[idx, 0].set_xticks(np.arange(50, 81, 10))
        axs[idx, 0].set_yticks(np.arange(50, 81, 10))
        axs[idx, 0].set_xlabel('Nodal Data', fontsize=8)
        axs[idx, 0].set_ylabel('Flightradar24 Data', fontsize=8)
        axs[idx, 0].tick_params(axis='both', labelsize=8)
        m, b = fit_l1_line(v0_new, speeds_org, bounds=(50, 80))
        x = np.linspace(min(v0_new), max(v0_new), 100)
        axs[idx, 0].plot(x, m * x + b, color='k')

        scatter2 = axs[idx, 1].scatter(distance_new, dists_org, c=temp_c, cmap='coolwarm', s=15)
        axs[idx, 1].set_title(f"{title[idx]}: Distance (m)", fontsize=10)
        axs[idx, 1].set_xlim(0, 2000)
        axs[idx, 1].set_ylim(0, 2000)
        axs[idx, 1].axline((0, 0), slope=1, color='black', linestyle='--')
        axs[idx, 1].set_aspect('equal', adjustable='box')
        axs[idx, 1].set_xticks(np.arange(0, 2001, 1000))
        axs[idx, 1].set_yticks(np.arange(0, 2001, 1000))
        axs[idx, 1].set_xlabel('Nodal Data', fontsize=8)
        axs[idx, 1].set_ylabel('Flightradar24 Data', fontsize=8)
        axs[idx, 1].tick_params(axis='both', labelsize=8)
        m, b = fit_l1_line(distance_new, dists_org, bounds=(0, 2000))
        x = np.linspace(min(distance_new), max(dists_org), 100)
        axs[idx, 1].plot(x, m * x + b, color='k')


        scatter3 = axs[idx, 2].scatter(np.array(time_new), np.array(times_org), c=temp_c, cmap='coolwarm', s=15)
        axs[idx, 2].set_title(f"{title[idx]}: Time (s)", fontsize=10)
        axs[idx, 2].set_xlim(110, 120)
        axs[idx, 2].set_xlabel('Nodal Data', fontsize=8)
        axs[idx, 2].set_ylabel('Flightradar24 Data', fontsize=8)
        axs[idx, 2].tick_params(axis='both', labelsize=8)

        cbar = fig.colorbar(scatter1, ax=axs[idx, 2], orientation='vertical', pad=0.1)
        cbar.set_label('Temperature (°C)')
        color_at_minus_2 = cbar.cmap(cbar.norm(-2))
        c_temp_array = scatter2
    else:
        scatter1 = axs[idx, 0].scatter(v0_new, speeds_org, c=color_at_minus_2, s=15)
        axs[idx, 0].set_title(f"{title[idx]}: Velocity (m/s)", fontsize=10)
        axs[idx, 0].set_xlim(50, 80)
        axs[idx, 0].axline((0, 0), slope=1, color='black', linestyle='--')
        axs[idx, 0].set_ylim(50, 80)
        axs[idx, 0].set_aspect('equal')
        axs[idx, 0].set_xticks(np.arange(50, 81, 10))
        axs[idx, 0].set_yticks(np.arange(50, 81, 10))
        axs[idx, 0].set_xlabel('Nodal Data', fontsize=8)
        axs[idx, 0].set_ylabel('Flightradar24 Data', fontsize=8)
        axs[idx, 0].tick_params(axis='both', labelsize=8)
        m, b = fit_l1_line(v0_new, speeds_org, bounds=(50, 80))
        x = np.linspace(min(v0_new), max(v0_new), 100)
        axs[idx, 0].plot(x, m * x + b, color='k')

        scatter2 = axs[idx, 1].scatter(distance_new, dists_org, c=color_at_minus_2, s=15)
        axs[idx, 1].set_title(f"{title[idx]}: Distance (m)", fontsize=10)
        axs[idx, 1].set_xlim(0, 2000)
        axs[idx, 1].set_ylim(0, 2000)
        axs[idx, 1].axline((0, 0), slope=1, color='black', linestyle='--')
        axs[idx, 1].set_aspect('equal', adjustable='box')
        axs[idx, 1].set_xticks(np.arange(0, 2001, 1000))
        axs[idx, 1].set_yticks(np.arange(0, 2001, 1000))
        axs[idx, 1].set_xlabel('Nodal Data', fontsize=8)
        axs[idx, 1].set_ylabel('Flightradar24 Data', fontsize=8)
        axs[idx, 1].tick_params(axis='both', labelsize=8)
        m, b = fit_l1_line(distance_new, dists_org, bounds=(0, 2000))
        x = np.linspace(min(distance_new), max(dists_org), 100)
        axs[idx, 1].plot(x, m * x + b, color='k')

        scatter3 = axs[idx, 2].scatter(np.array(time_new), np.array(times_org), c=color_at_minus_2, s=15)
        axs[idx, 2].set_title(f"{title[idx]}: Time (s)", fontsize=10)
        axs[idx, 2].set_xlim(110, 120)
        axs[idx, 2].set_xlabel('Nodal Data', fontsize=8)
        axs[idx, 2].set_ylabel('Flightradar24 Data', fontsize=8)
        axs[idx, 2].tick_params(axis='both', labelsize=8)

        # Add a single colorbar for the entire figure
        cbar = fig.colorbar(c_temp_array, ax=axs[idx, 2], orientation='vertical', pad=0.1)
        cbar.set_label('Temperature (°C)')
plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.95])  # Adjust margins for better spacing
plt.subplots_adjust(hspace=0.4, wspace=0.3)  # Increase spacing between subplots

# Ensure all plots are square, with equal aspect ratio for Velocity and Distance plots
for ax_row in axs:
    for i, ax in enumerate(ax_row):
        if i != 2:  # Velocity and Distance plots
            ax.set_aspect('equal', adjustable='box')
        else:  # Time plots
            ax.set_aspect('auto')  # Allow auto aspect ratio for Time plots
            ax.set_box_aspect(1)  # Make Time plots square
plt.show()
