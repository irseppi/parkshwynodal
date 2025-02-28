import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = open('full_atmc_v_2c.txt','r')

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

for line in file.readlines():
    y += 1
    counts = []

    lines = line.split(',')
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
    plt.scatter(np.abs(np.array(v0_new) - np.array(v0_old)), date, c='b') 
    plt.show()

    plt.figure()
    plt.scatter(distance_old, date, c='b') 
    plt.scatter(distance_new, date, c='r')
    plt.title('distance')
    plt.show()

    plt.figure()
    plt.scatter(np.abs(np.array(distance_new) - np.array(distance_old)), date, c='b') 
    plt.show()

    plt.figure()
    plt.scatter(time_old, date, c='b') 
    plt.scatter(time_new, date, c='r')
    plt.title('time')
    plt.show()


    plt.figure()
    plt.scatter(np.abs(np.array(time_new) - np.array(time_old)), date, c='b') 
    plt.title('time')
    plt.show()

elif option == 2:
    plt.figure()
    plt.scatter(delf0_old, v0_old, c='b')
    plt.scatter(delf0_new, v0_new, c='r')
    plt.xlabel('delf0')
    plt.ylabel('v0')
    plt.show()
