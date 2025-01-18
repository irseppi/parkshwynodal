import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = open('old_new_invers_1o.txt','r')

file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

# Create a dictionary to store the color for each tail number
color_dict = {}

count_dict = {}
date_dict = {} # Add a list to store dates
peaks_dict_old = {}
peaks_dict_new = {}
date_old = []
date_new = []
y=0
time_new = []
time_old = []

v0_old = []
v0_new = []

pppp_old = []
pppp_new = []

distance_old = []
distance_new = []

plt.figure()
# Iterate over each line in the file
for line in file.readlines():
    y += 1
    counts = []

    lines = line.split(',')
    time_old.append(float(lines[3]))
    time_new.append(float(lines[8]))
    v0_old.append(float(lines[5]))
    v0_new.append(float(lines[10]))
    distance_old.append(float(lines[6]))
    distance_new.append(float(lines[10]))

    peaks_old = np.array(lines[7])

    peaks_old = str(peaks_old)
    peaks_old = np.char.replace(peaks_old, '[', '')
    peaks_old = np.char.replace(peaks_old, ']', '')

    peaks_old = str(peaks_old)
    peaks_old = np.array(peaks_old.split(' '))

    peaks_new = np.array(lines[12])

    peaks_new= str(peaks_new)
    peaks_new = np.char.replace(peaks_new, '[', '')
    peaks_new = np.char.replace(peaks_new, ']', '')

    peaks_new = str(peaks_new)
    peaks_new = np.array(peaks_new.split(' '))
    
    for i in range(len(peaks_old)):
        peak_old = peaks_old[i]
        pppp_old.append(float(peak_old))
        date_old.append(y)
    for i in range(len(peaks_new)):
        peak_new = peaks_new[i]
        pppp_new.append(float(peak_new))
        date_new.append(y)


plt.scatter(pppp_old, date_old, c='b') #c=color, label=tail_num)

plt.scatter(pppp_new, date_new, c='r')#, c=color)




plt.show()
