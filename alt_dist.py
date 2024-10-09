import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj

file = open('output/C185data_updated.txt', 'r')
file2 = pd.read_csv('input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flights = file2['FLIGHT_NUM']
utm_proj = pyproj.Proj(proj='utm', zone='6', ellps='WGS84')

list_f = []
color_dict = {}
label_dict = {}
MODE = 2
plt.figure()
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = float(lines[1])

    peaks = np.array(lines[7])

    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))
    if MODE == 1:
        for peak in peaks:
            try:
                peak = float(peak[0:-1])
            except:
                continue
            if abs(float(peak) - 82.5) <= 1.5 or abs(float(peak) - 124.5) <= 1.5:
                title = 'Mode 1'

            elif abs(float(peak) - 76) <= 5 or abs(float(peak) - 114) <= 5 or abs(float(peak) - 152) <= 5:
                continue
    elif MODE == 2:
        for peak in peaks:
            try:
                peak = float(peak[0:-1])
            except:
                continue
            if abs(float(peak) - 82.5) <= 1.5 or abs(float(peak) - 124.5) <= 1.5:
                continue

            elif abs(float(peak) - 76) <= 5 or abs(float(peak) - 114) <= 1.5 or abs(float(peak) - 152) <= 5:
                title = 'Mode 2'

    if flight_num in list_f:
        continue
    else:
        list_f.append(flight_num)

    for lp in range(len(flights)):
        if int(flight_num) == int(flights[lp]):
            tail_num = tail_nums[lp]
            if tail_num not in color_dict:
                color_dict[tail_num] = np.random.rand(3,)
    date = lines[0]
    flight_num = str(int(flight_num))
    month = date[4:6]
    day = date[6:8]

    flight_file = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+ '_positions/2019'+month+day+ '_' + flight_num + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")

    altitude = flight_data['altitude']
    lat = flight_data['latitude']
    lon = flight_data['longitude']
    alts = []
    dists = []
    for a in range(1,len(altitude)-1):
        if a == 1:
            dist  = 0
            x1,y1 = utm_proj(lon[a], lat[a])
        else:
            x,y = utm_proj(lon[a], lat[a])
            dist = np.sqrt((x-x1)**2 + (y-y1)**2) + dist
            x1 = x
            y1 = y

        alts.append(altitude[a])
        dists.append(dist/1000)
    if tail_num not in label_dict:
        plt.plot(dists, alts, color=color_dict[tail_num], label=tail_num)
        label_dict[tail_num] = tail_num
    else:
        plt.plot(dists, alts, color=color_dict[tail_num])
plt.legend() # label=label_dict[tail_num])
plt.title(title)
plt.show()
