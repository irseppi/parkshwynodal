import matplotlib.pyplot as plt
import numpy as np


file = open('output4.txt', 'r')
all_med = []

for line in file.readlines():
    lines = line.split(',')
    #quality_num = int(lines[8])

    peaks = lines[7]

    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))
    peaks_filt = []

    for peak in peaks:
        if peak == '' or peak == ' ' or peak == '   ':
            continue
        else:
            peaks_filt.append(float(peak))
    peaks_filt = sorted(peaks_filt)

    f1 = []
    for p in range(len(peaks_filt)):
        if p == 0:
            continue
        
        diff = float(peaks_filt[p]) - float(peaks_filt[p-1])
        if diff > 21 or diff < 18:
            continue
            diff= diff/20
        f1.append(diff)

    if not np.isnan(np.median(f1)):
        all_med.append(np.median(f1))

plt.figure()
plt.hist(all_med,bins=30) #, bins=50, color='blue')
plt.axvline(x=np.median(all_med), color='red', linestyle='--')


plt.text(20, 40, 'Median: '+str(np.round(np.median(all_med),2)))


#plt.xlim(0,140)
plt.show()
