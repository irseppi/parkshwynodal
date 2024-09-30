
import matplotlib.pyplot as plt
import numpy as np

file = open('output/C185data_updated.csv', 'r')
file2 = open('input/all_station_crossing_db_C185.txt', 'r')
plt.figure()

for line in file.readlines():
    lines = line.split(',')

    peaks = np.array(lines[7])

    peaks = str(peaks)

    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks) 
    peaks = np.array(peaks.split(' '))
   
    for lp in file2.readlines():
        le = lp.split(',')
        
        if int(lines[1]) == int(le[1]):
            tail_num = le[2]
        else:
            continue 


    for peak in peaks:
        print(peak)
        if peak == '' or peak == ' ':
            continue
        plt.scatter(peak, tail_num)


plt.show()

