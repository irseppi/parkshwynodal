import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = open('2_ex_diff.txt','r')

file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

# Create a dictionary to store the color for each tail number
color_dict = {}

count_dict = {}
date_dict = {} # Add a list to store dates
peaks_dict_old = {}
peaks_dict_new = {}

y=0
time_new = []
time_old = []

v0_old = []
v0_new = []

distance_old = []
distance_new = []

plt.figure()
# Iterate over each line in the file
for line in file.readlines():
    y += 1
    counts = []
    pppp_old = []
    pppp_new = []
    dates = []

    lines = line.split(',')
    time_old.append(float(lines[3]))
    time_new.append(float(lines[8]))
    v0_old.append(float(lines[5]))
    v0_new.append(float(lines[10]))
    flight_num = float(lines[1])
    distance_old.append(float(lines[6]))
    distance_new.append(float(lines[10]))
    date = lines[3]

    peaks_old = np.array(lines[7])

    peaks_old = str(peaks_old)
    peaks_old = np.char.replace(peaks_old, '[', '')
    peaks_old = np.char.replace(peaks_old, ']', '')

    peaks_old = str(peaks_old)
    peaks_old = np.array(peaks_old.split(' '))

    peaks_new = np.array(lines[11])

    peaks_new= str(peaks_new)
    peaks_new = np.char.replace(peaks_new, '[', '')
    peaks_new = np.char.replace(peaks_new, ']', '')

    peaks_new = str(peaks_new)
    peaks_new = np.array(peaks_new.split(' '))
    
    for i in range(len(peaks_old)):
        peak_old = peaks_old[i]
        pppp_old.append(float(peak_old))

    for i in range(len(peaks_new)):
        peak_new = peaks_new[i]
        pppp_new.append(float(peak_new))

    if len(pppp_old) > len(pppp_new):
        # Find the two closest values
        closest_values = min(zip(pppp_old[:-1], pppp_old[1:]), key=lambda x: abs(x[0] - x[1]))

        # Calculate the average of the two closest values
        average = sum(closest_values) / 2

        # Replace the two closest values with the average
        pppp_old.remove(closest_values[0])
        pppp_old.remove(closest_values[1])
        pppp_old.append(average)
    elif len(pppp_new) > len(pppp_old):
        # Find the two closest values
        closest_values = min(zip(pppp_new[:-1], pppp_new[1:]), key=lambda x: abs(x[0] - x[1]))

        # Calculate the average of the two closest values
        average = sum(closest_values) / 2

        # Replace the two closest values with the average
        pppp_new.remove(closest_values[0])
        pppp_new.remove(closest_values[1])
        pppp_new.append(average)

    if len(pppp_old) == len(pppp_new):
        for i in range(len(pppp_new)):
            dates.append(float(date))
            counts.append(float(y))
    else:
        print(line)
        continue
    pppp_old = sorted(pppp_old)
    pppp_new = sorted(pppp_new)

    peak_diff = abs(np.array(pppp_new) - np.array(pppp_old))
    indices_to_remove = [i for i, diff in enumerate(peak_diff) if diff > 10]
    pppp_old = [num for i, num in enumerate(pppp_old) if i not in indices_to_remove]
    pppp_new = [num for i, num in enumerate(pppp_new) if i not in indices_to_remove]
    counts = [num for i, num in enumerate(counts) if i not in indices_to_remove]
    dates = [num for i, num in enumerate(dates) if i not in indices_to_remove]

    for lp in range(len(flight)):
        if int(flight_num) == int(flight[lp]):
            tail_num = tail_nums[lp]
            if tail_num not in color_dict:
                color_dict[tail_num] = np.random.rand(3,)
                peaks_dict_old[tail_num] = []
                peaks_dict_new[tail_num] = []
                date_dict[tail_num] = []
                count_dict[tail_num] = []
            date_dict[tail_num].extend(dates)
            peaks_dict_old[tail_num].extend(pppp_old)
            peaks_dict_new[tail_num].extend(pppp_new)
            count_dict[tail_num].extend(counts)
                
for tail_num, p_old in peaks_dict_old.items():
    date = date_dict[tail_num]
    y = count_dict[tail_num]
    color = color_dict[tail_num]
    plt.scatter(p_old, y, c='b') #c=color, label=tail_num)
    
    p_new = peaks_dict_new[tail_num]
    plt.scatter(p_new, y, c='r')#, c=color)
    
    # Plot a dashed line between p_old and p_new
    plt.plot([p_old, p_new], [y, y], linestyle='--',color='orange') # color=color)

plt.legend(loc='upper left',fontsize = 'x-small')

plt.show()
