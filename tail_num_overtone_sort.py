import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file = open('output4.txt', 'r')
file2 = pd.read_csv('input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']
time = file2['TIME']
plt.figure()

# Create a dictionary to store the color for each tail number
color_dict = {}
y_pos_dict = {}

# Create a dictionary to store the peaks for each tail number
peaks_dict = {}
date_dict = {}
count = 0
xx = 'yes'
all_med = []
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = lines[1]
    #quality_num = int(lines[8]) #lines[9] for output2
    nodes = int(lines[2])
    #if quality_num < 21:
    #    continue

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
        print(peak)
        if peak == '' or peak == ' ' or peak == '   ':
            continue
        else:
            if len(peaks) == 0  or peak == peaks[0]:
                peak_old = peak
                continue
            diff = float(peak) - float(peak_old)
            if diff > 21 or diff < 18:
                continue
            f1.append(diff)
            ppp.append(float(peak))
            date.append(float((lines[3])))
        y.append(count)
        peak_old = peak
    peaks = ppp

    

    if not np.isnan(np.median(f1)):
        all_med.append(np.median(f1))

    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            # Assign a color to the tail number if it doesn't already have one
            if tail_num not in color_dict:
                color_dict[tail_num] = np.random.rand(3,)

            if tail_num not in peaks_dict:
                peaks_dict[tail_num] = []
                date_dict[tail_num] = []
                y_pos_dict[tail_num] = []
            peaks_dict[tail_num].extend(ppp)
            date_dict[tail_num].extend(date)
            y_pos_dict[tail_num].extend(y)
# Plot the data peaks vs their date and color code by tail number
for tail_num, peaks in peaks_dict.items():
    color = color_dict[tail_num]
    dates = date_dict[tail_num]
    y =  y_pos_dict[tail_num]
    #plt.scatter(peaks, dates, c=color,label=tail_num) #,color=color)
    #plt.scatter(peaks, y, c=color,label=tail_num) #,color=color)
    #plt.scatter(peaks, y, c=dates) #,color=color)
    plt.scatter(all_med, y, c=color,label=tail_num)
plt.legend()

plt.show()
