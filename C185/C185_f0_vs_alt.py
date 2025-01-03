import os
import matplotlib.pyplot as plt
import numpy as np

# Define the directory where your files are located
file = open('output3.csv', 'r')
input = open('input/all_station_crossing_db_C185.txt', 'r')

alt_pre = []
fly_pre = []
sta_pre = []
time_pre = []
for lp in input.readlines():
    le = lp.split(',')
    fly_pre.append(le[1])
    sta_pre.append(le[6])
    try:
        time_pre.append(float(le[3]))
        alt_pre.append(float(le[4]))
    except:
        time_pre.append(le[3])
        alt_pre.append(le[4])
plt.figure()

x = 0

altitudes = []  # Store altitudes

for line in file.readlines():
    x += 1
    lines = line.split(',')

    flight = lines[1]
    station = lines[2]

    for tp in range(len(fly_pre)):
        fly = fly_pre[tp]
        alt = alt_pre[tp]
        try:
            alt = float(alt)
        except:
            continue
        sta = sta_pre[tp]
        time = time_pre[tp]
        if fly == flight and sta == station:
            try:
                alt_mode = (alt_pre[tp+1] - alt_pre[tp-1])/(time_pre[tp+1] - time_pre[tp-1])
            except:
                alt_mode = 0
            break
        else:
            alt_mode = 0
            continue

    if tp != len(fly_pre):
        try:
            if x <= 29:
                peaks = np.array(lines[6])
            else:
                peaks = np.array(lines[7])
        except:
            continue
        peaks = str(peaks)  
        peaks = np.array(peaks.split(' '))
        for peak in peaks:
            try:
                peak = float(peak[0:-1])
            except:
                continue
            if abs(float(peak) - 82.5) <= 1.5 or abs(float(peak) - 124.5) <= 1.5:
                mode = 1
                break
            elif abs(float(peak) - 76) <= 5 or abs(float(peak) - 114) <= 6 or abs(float(peak) - 152) <= 6:
                mode = 2
                break
            else:
                mode = 0
                continue
    else:
        continue
    print(mode, alt_mode)
    plt.scatter(mode, alt_mode, color='black')
    altitudes.append(alt_mode)  # Add altitude to the list

sorted_altitudes = sorted(altitudes)  # Sort altitudes
#plt.yticks(np.arange(0,4500, 200))  # Set sorted altitudes on y-axis
plt.xticks([0, 1, 2])  # Set sorted altitudes on x-axis
plt.xlabel('Mode')
plt.ylabel('Change in Altitude (f/s)')
plt.xticks([0, 1, 2])
plt.show()


