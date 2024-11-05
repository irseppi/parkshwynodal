import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file = open('output4.txt', 'r')
file2 = pd.read_csv('input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']
time = file2['TIME']

# Create a dictionary to store the color for each tail number
color_dict = {}
y_pos_dict = {}

# Create a dictionary to store the peaks for each tail number
peaks_dict = {}
date_dict = {}

count = 0
xx = 'yes'
all_med = {}
y_med = {}

# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = lines[1]
    nodes = int(lines[2])

    peaks = np.array(lines[7])

    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))

    count += 1
    ppp = []
    date = []
    y=[]
    f1 = []

    for peak in peaks:
        if peak == '' or peak == ' ' or peak == '   ':
            continue
        else:
            ppp.append(float(peak))
            date.append(float((lines[3])))
            y.append(count)
            if len(peaks) == 0 or peak == peaks[0]:
                peak_old = peak
                continue
            diff = float(peak) - float(peak_old)
            #if diff > 21 or diff < 18:
            #    continue
            f1.append(diff)
        peak_old = peak

    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            # Assign a color to the tail number if it doesn't already have one
            if tail_num not in color_dict:
                color_dict[tail_num] = np.random.rand(3,)

                peaks_dict[tail_num] = []
                y_pos_dict[tail_num] = []
                date_dict[tail_num] = []

                all_med[tail_num] = []
                y_med[tail_num] = []

            peaks_dict[tail_num].extend(ppp)
            date_dict[tail_num].extend(date)
            y_pos_dict[tail_num].extend(y)

            if not np.isnan(np.median(f1)):
                all_med[tail_num].extend([np.median(f1)])
                y_med[tail_num].extend([count])

fig,ax1 = plt.subplots(1, 1, sharex=False, figsize=(8,6))     

ax1.margins(x=0)
ax1.grid(axis='x') 
ax2 = fig.add_axes([0.90, 0.11, 0.07, 0.77], sharey=ax1) 

for tail_num, peaks in peaks_dict.items():
    color = color_dict[tail_num]
    dates = date_dict[tail_num]
    y =  y_pos_dict[tail_num]
    ax1.scatter(peaks, y, c=color,label=tail_num) 
for tail_num, med in all_med.items():
    print(tail_num)
    y_med = y_med[tail_num]
    ax2.scatter(med, y_med, c=color)  

ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax2.grid(axis='x') 
ax1.legend(loc='upper right',fontsize = 'x-small')
ax1.set_xlim(100, 280)
plt.show()
