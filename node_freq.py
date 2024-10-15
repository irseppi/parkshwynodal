import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = open('output2.txt', 'r')

counts = []
pppp = []
dates = []
nodes = []

head_dict = {}
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
    flight = str(lines[1])
    node = lines[2]
    date = lines[3]
    nodes.append(node)
    month = lines[0][4:6]
    day = lines[0][6:8]
    flight_file = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+ '_positions/2019'+month+day+ '_' + flight + '.csv'
    flight_data = pd.read_csv(flight_file, sep=",")
    t = flight_data['snapshot_id']
    head = flight_data['heading']
    check = np.full((len(t),),int(float(date)))
    
    ind = np.argmin(abs(t - check))
    head = head[ind]    
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
    
    if node not in peaks_dict:
        peaks_dict[node] = []
        date_dict[node] = []
        count_dict[node] = []
    date_dict[node].extend(dates)
    peaks_dict[node].extend(pppp)
    count_dict[node].extend(counts)

for node, peaks in peaks_dict.items():
    date = date_dict[node]
    y = count_dict[node]
    #print(y[-1])
    plt.scatter(peaks, y, c=date, cmap='rainbow')
    #plt.axhline(y=(int(y[-1])+0.5), color='black', linewidth=0.5)

# Set the y-axis ticks and labels
#plt.yticks(range(len(nodes)), nodes)

plt.colorbar()

plt.show()
