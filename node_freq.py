import matplotlib.pyplot as plt
import numpy as np

file = open('output/C185data_updated.txt', 'r')

list_f = []
nodes = []
pppp = []


plt.figure()
# Iterate over each line in the file
for line in file.readlines():
    lines = line.split(',')
    flight_num = float(lines[1])
    node = lines[2]


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
            nodes.append(str(node))



plt.scatter(pppp, nodes)
plt.yticks(nodes)
plt.show()