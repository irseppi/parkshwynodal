import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = open('C185data_atmosphere_1o.txt','r')

file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

# Create a dictionary to store the color for each tail number
color_dict = {}

counts = []
pppp = []
dates = []

count_dict = {}
date_dict = {} # Add a list to store dates
peaks_dict = {}

y=0

plt.figure()
# Iterate over each line in the file
for line in file.readlines():
    y += 1
    
    lines = line.split(',')
    flight_num = float(lines[1])
    date = lines[3]

 
    peaks = np.array(lines[7])

    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))

    for peak in peaks:
        if peak == '' or peak == ' ' or peak == '   ':
            continue
        else:
            pppp.append(float(peak))
            dates.append(float(date))
            counts.append(float(y))
    
    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            if tail_num not in color_dict:
                color_dict[tail_num] = np.random.rand(3,)
                peaks_dict[tail_num] = []
                date_dict[tail_num] = []
                count_dict[tail_num] = []
            date_dict[tail_num].extend(dates)
            peaks_dict[tail_num].extend(pppp)
            count_dict[tail_num].extend(counts)

for tail_num, peaks in peaks_dict.items():
    date = date_dict[tail_num]
    y = count_dict[tail_num]
    color = color_dict[tail_num]
    plt.scatter(peaks, y, c=color,label=tail_num)

plt.legend(loc='upper left',fontsize = 'x-small')

plt.show()
