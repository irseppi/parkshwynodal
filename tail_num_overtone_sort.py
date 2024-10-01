import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file = open('output/C185data_updated.txt', 'r')
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
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = lines[1]
    
    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            if tail_num == 10512184 or tail_num == 11013232 or tail_num == 11146232:
                xx = 'no'
                continue
    if xx == 'no':
        xx = 'yes'
        continue
    peaks = np.array(lines[7])

    peaks = str(peaks)
    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))

    count += 1
    ppp = []
    date = []

    for peak in peaks:
        if peak == '' or peak == ' ' or peak == '   ':
            continue
        else:
            ppp.append(float(peak))
            date.append(float(count))

    peaks = ppp

    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            if tail_num == 10512184 or tail_num == 11013232 or tail_num == 11146232:
                continue
            # Assign a color to the tail number if it doesn't already have one
            if tail_num not in color_dict:
                color_dict[tail_num] = np.random.rand(3,)

            if tail_num not in peaks_dict:
                peaks_dict[tail_num] = []
                date_dict[tail_num] = []

            peaks_dict[tail_num].extend(ppp)
            date_dict[tail_num].extend(date)

# Plot the data peaks vs their date and color code by tail number
for tail_num, peaks in peaks_dict.items():
    color = color_dict[tail_num]
    dates = date_dict[tail_num]

    plt.scatter(peaks, dates, label=tail_num, color=color)

plt.legend()

plt.show()
