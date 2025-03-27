import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from pyproj import Proj
from matplotlib.gridspec import GridSpec

# Initialize UTM projection for coordinate conversion
utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

# Read input file containing station crossing data
file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt', 'r')
latc = []  # List to store latitudes
lonc = []  # List to store longitudes
timec = []  # List to store times
altc = []  # List to store altitudes

# Parse the input file line by line
for line in file_in.readlines():
    text = line.split(',')

    # Extract time and UTM coordinates
    timec.append(float(text[5]))
    x = float(text[2])
    y = float(text[3])

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)
    latc.append(lat)
    lonc.append(lon)
    altc.append(float(text[4]) * 0.0003048)  # Convert altitude from feet to kilometers
file_in.close()

# List of files to process and their corresponding titles
file_list = ['1o_atmc_v_2c.txt', '2c_1o_v_full.txt', 'full_atmc_v_2c.txt', 'atmc_1o_v_full.txt']
title = [
    'One Harmonic/Fixed Temp v One Harmonic/Varying Temp',
    'One Harmonic/Fixed Temp v Full Harmonics/Fixed Temp',
    'Full Harmonics/Fixed Temp v Full Harmonics/Varying Temp',
    'One Harmonic/Varying Temp v Full Harmonics/Varying Temp'
]

# Process each file in the list
for gg, file in enumerate(file_list):
    file = open(file, 'r')
    date = []  # List to store indices
    y = 0  # Counter for lines
    time_new = []  # List for new times
    time_old = []  # List for old times

    v0_old = []  # List for old velocities
    v0_new = []  # List for new velocities

    distance_old = []  # List for old distances
    distance_new = []  # List for new distances

    delf0_new = []  # Placeholder for new frequency differences
    delf0_old = []  # Placeholder for old frequency differences

    pppp_old = []  # List for old peaks
    pppp_new = []  # List for new peaks
    med_old = []  # List for median old frequency differences
    med_new = []  # List for median new frequency differences
    date_old = []  # List for old peak indices
    date_new = []  # List for new peak indices
    date_all = []  # List for all indices
    temp_c = []  # List for temperatures

    # Parse the file line by line
    for line in file.readlines():
        y += 1
        counts = []

        lines = line.split(',')
        time = float(lines[3])

        # Match the time with the corresponding latitude, longitude, and altitude
        for i in range(len(latc)):
            if timec[i] == time:
                lat = latc[i]
                lon = lonc[i]
                alt = altc[i]
                break

        # Construct the input file path for atmospheric data
        input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat'
        try:
            file = open(input_files, 'r')
        except:
            continue
        data = json.load(file)

        # Extract metadata and data from the atmosphere JSON file
        metadata = data['metadata']
        sourcefile = metadata['sourcefile']
        datetim = metadata['time']['datetime']
        latitude = metadata['location']['latitude']
        longitude = metadata['location']['longitude']
        parameters = metadata['parameters']

        data_list = data['data']

        # Convert data to a pandas DataFrame
        data_frame = pd.DataFrame(data_list)

        # Find the closest altitude index
        z_index = None
        hold = np.inf
        for item in data_list:
            if item['parameter'] == 'Z':
                for i in range(len(item['values'])):
                    if abs(float(item['values'][i]) - float(alt)) < hold:
                        hold = abs(float(item['values'][i]) - float(alt))
                        z_index = i

        # Extract temperature at the closest altitude
        for item in data_list:
            if item['parameter'] == 'T':
                Tc = -273.15 + float(item['values'][z_index])  # Convert from Kelvin to Celsius
        temp_c.append(Tc)

        # Extract other parameters from the line
        time_old.append(float(lines[4]))
        time_new.append(float(lines[9]))
        v0_old.append(float(lines[5]))
        v0_new.append(float(lines[10]))
        flight_num = float(lines[1])
        distance_old.append(float(lines[6]))
        distance_new.append(float(lines[11]))
        date.append(y)

        # Process old and new peaks
        peaks_old = np.array(lines[7])
        peaks_old = str(peaks_old)
        peaks_old = np.char.replace(peaks_old, '[', '')
        peaks_old = np.char.replace(peaks_old, ']', '')
        peaks_old = str(peaks_old)
        peaks_old = np.array(peaks_old.split(' '))

        peaks_new = np.array(lines[12])
        peaks_new = str(peaks_new)
        peaks_new = np.char.replace(peaks_new, '[', '')
        peaks_new = np.char.replace(peaks_new, ']', '')
        peaks_new = str(peaks_new)
        peaks_new = np.array(peaks_new.split(' '))

        # Calculate median frequency differences for old peaks
        f1 = []
        for i in range(len(peaks_old)):
            peak_old = peaks_old[i]
            pppp_old.append(float(peak_old))
            date_old.append(y)
            if i == 0:
                continue
            diff = float(peaks_old[i]) - float(peaks_old[i - 1])
            if diff > 22 or diff < 18:
                continue
            f1.append(diff)
        med_old.append(np.nanmedian(f1))

        # Calculate median frequency differences for new peaks
        f2 = []
        for i in range(len(peaks_new)):
            peak_new = peaks_new[i]
            pppp_new.append(float(peak_new))
            date_new.append(y)
            if i == 0:
                continue
            diff = float(peaks_new[i]) - float(peaks_new[i - 1])
            if diff > 22 or diff < 18:
                continue
            f2.append(diff)
        med_new.append(np.nanmedian(f2))
        date_all.append(y)

    # Create plots for the current file
    fig = plt.figure(figsize=(10, 12))
    fig.suptitle(title[gg], fontsize=16)

    # Define grid layout for subplots
    gs = GridSpec(3, 2, figure=fig, width_ratios=[5, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)

    # Scatter plot to compare old and new peaks (ax1) and their median frequency differences (ax2)
    # ax1: Visualizes the old and new peaks with their respective indices
    # ax2: Shows the median frequency differences for old and new peaks
    ax1.margins(x=0)
    ax1.scatter(pppp_old, date_old, c='b', label='Old Peaks')
    ax1.scatter(pppp_new, date_new, c='r', label='New Peaks')

    # Scatter plot for median frequency differences
    ax2.scatter(med_old, date_all, c='b')
    ax2.scatter(med_new, date_all, c='r')

    # Configure axis labels and limits
    ax1.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
    ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
    ax1.set_xlabel('Frequency')
    ax2.set_xlabel('\u0394' + 'F')
    ax1.legend(loc='upper left', fontsize='x-small')
    ax1.set_xlim(0, 300)
    ax1.set_xticks(range(0, 251, 25))

    # Create additional subplots for velocity, distance, and time
    axs = fig.subplots(2, 3, gridspec_kw={'height_ratios': [1, 1], 'top': 0.55, 'hspace': 0.3, 'wspace': 0.15})

    # Velocity scatter plots
    axs[0, 0].scatter(v0_old, date, c='b')
    axs[0, 0].scatter(v0_new, date, c='r')
    axs[0, 0].set_title('v0')
    axs[0, 0].set_ylabel('Index')

    scatter = axs[1, 0].scatter((np.array(v0_new) - np.array(v0_old)), date, c=temp_c, cmap='coolwarm', label='Velocity Residuals')
    axs[1, 0].set_title("Velocity Residuals")
    axs[1, 0].set_ylabel('Index')

    # Distance scatter plots
    axs[0, 1].scatter(distance_old, date, c='b')
    axs[0, 1].scatter(distance_new, date, c='r')
    axs[0, 1].set_title('Distance')
    axs[0, 1].tick_params(left=False, labelleft=False)

    scatter = axs[1, 1].scatter((np.array(distance_new) - np.array(distance_old)), date, c=temp_c, cmap='coolwarm', label='Distance Residuals')
    axs[1, 1].set_title("Distance Residuals")
    axs[1, 1].tick_params(left=False, labelleft=False)

    # Time scatter plots
    axs[0, 2].scatter(time_old, date, c='b')
    axs[0, 2].scatter(time_new, date, c='r')
    axs[0, 2].set_title('Time')
    axs[0, 2].tick_params(left=False, labelleft=False)

    scatter = axs[1, 2].scatter((np.array(time_new) - np.array(time_old)), date, c=temp_c, cmap='coolwarm', label='Time Residuals')
    axs[1, 2].set_title("Time Residuals")
    axs[1, 2].tick_params(left=False, labelleft=False)

    # Add colorbar for temperature
    fig.colorbar(scatter, ax=axs[1, 2], orientation='vertical', label='Temperature (°C)')

    # Adjust layout and display the plot
    plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95], h_pad=0.5, w_pad=0.5)
    plt.show()
