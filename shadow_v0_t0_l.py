import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
from pyproj import Proj

file = open('1o_atmc_v_2c.txt','r')
#file = open('2c_1o_v_full.txt','r')

date = []
y=0
time_new = []
time_old = []

v0_old = []
v0_new = []

distance_old = []
distance_new = []

delf0_new = []
delf0_old= []
option = 1

utm_proj = Proj(proj='utm', zone='6', ellps='WGS84')

file_in = open('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_UTM.txt','r')
latc = []
lonc = []
timec = []
altc = []
for line in file_in.readlines():
    text = line.split(',')

    timec.append(float(text[5]))
    x =  float(text[2])  # Replace with your UTM x-coordinate
    y = float(text[3])  # Replace with your UTM y-coordinate

    # Convert UTM coordinates to latitude and longitude
    lon, lat = utm_proj(x, y, inverse=True)
    latc.append(lat)
    lonc.append(lon)
    altc.append(float(text[4])*0.0003048) #convert between feet and km
file_in.close()

temp_c = []
for line in file.readlines():
    y += 1
    counts = []

    lines = line.split(',')
    time = float(lines[3])
    for i in range(len(latc)):
        if timec[i] == time:
            lat = latc[i]
            lon = lonc[i]
            alt = altc[i]
            break

    input_files = '/scratch/irseppi/nodal_data/plane_info/atmosphere_data/' + str(time) + '_' + str(lat) + '_' + str(lon) + '.dat'
    try:
        file =  open(input_files, 'r') #as file:
    except:
        print('No file for: ', date, flight_num)
        continue
    data = json.load(file)

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
    temp_c.append(Tc)
    if option == 1:
        time_old.append(float(lines[4]))
        time_new.append(float(lines[9]))
        v0_old.append(float(lines[5]))
        v0_new.append(float(lines[10]))
        flight_num = float(lines[1])
        distance_old.append(float(lines[6]))
        distance_new.append(float(lines[11]))
        date.append(y)
    elif option == 2:
        peaks_new = np.array(lines[12])
        peaks_new = str(peaks_new)
        peaks_new = np.char.replace(peaks_new, '[', '')
        peaks_new = np.char.replace(peaks_new, ']', '')
        peaks_new = str(peaks_new)
        peaks_new = np.array(peaks_new.split(' '))
        f1 = []
        for p in range(len(peaks_new)):
            if p == 0:
                continue
            diff = float(peaks_new[p]) - float(peaks_new[p-1])
            if diff > 21 or diff < 18:
                continue
            f1.append(diff)
        if not np.isnan(np.median(f1)):
            delf0_new.append(np.median(f1))
            v0_new.append(float(lines[10]))

        peaks_old = np.array(lines[7])
        peaks_old = str(peaks_old)
        peaks_old = np.char.replace(peaks_old, '[', '')
        peaks_old = np.char.replace(peaks_old, ']', '')
        peaks_old = str(peaks_old)
        peaks_old = np.array(peaks_old.split(' '))
        f1 = []
        for p in range(len(peaks_old)):
            if p == 0:
                continue
            diff = float(peaks_old[p]) - float(peaks_old[p-1])
            if diff > 21 or diff < 18:
                continue
            f1.append(diff)
        if not np.isnan(np.median(f1)):
            delf0_old.append(np.median(f1))
            v0_old.append(float(lines[5]))
        
if option == 1:
    plt.figure()
    plt.scatter(v0_old, date, c='b')
    plt.scatter(v0_new, date, c='r')
    plt.title('v0')
    plt.show()

    plt.figure()
    plt.scatter((np.array(v0_new) - np.array(v0_old)), date, c=temp_c, cmap='coolwarm') 
    plt.colorbar(label='Temperature (°C)')
    plt.title("Velocity Residuals Between Fixed and Corrected Temeratures")
    plt.show()

    plt.figure()
    plt.scatter(distance_old, date, c='b') 
    plt.scatter(distance_new, date, c='r') 
    plt.title('distance')
    plt.show()

    plt.figure()
    plt.scatter((np.array(distance_new) - np.array(distance_old)), date, c=temp_c, cmap='coolwarm') 
    plt.colorbar(label='Temperature (°C)')
    plt.title("Distance Residuals Between Fixed and Corrected Temeratures")
    plt.show()

    plt.figure()
    plt.scatter(time_old, date, c='b') 
    plt.scatter(time_new, date, c='r')
    plt.title('time')
    plt.show()


    plt.figure()
    plt.scatter((np.array(time_new) - np.array(time_old)), date, c=temp_c, cmap='coolwarm') 
    plt.colorbar(label='Temperature (°C)')
    plt.title("Time Residuals Between Fixed and Corrected Temeratures")
    plt.show()

elif option == 2:
    plt.figure()
    plt.scatter(delf0_old, v0_old, c='b')
    plt.scatter(delf0_new, v0_new, c='r')
    plt.xlabel('delf0')
    plt.ylabel('v0')
    plt.show()
