import matplotlib.pyplot as plt
import numpy as np

file = open('output/C185data_updated.txt', 'r')
f1 = []
peaks_filt = []

for line in file.readlines():
    lines = line.split(',')
    quality_num = int(lines[9])

    peaks = lines[7]

    peaks = np.char.replace(peaks, '[', '')
    peaks = np.char.replace(peaks, ']', '')

    peaks = str(peaks)
    peaks = np.array(peaks.split(' '))

    for peak in peaks:
        if peak == '' or peak == ' ' or peak == '   ':
            continue
        else:
            peaks_filt.append(peak)

    for p in range(len(peaks_filt)):
        if p == 0:
            continue
        else:
            diff = float(peaks_filt[p]) - float(peaks_filt[p-1])
            f1.append(diff)


plt.figure()
plt.hist(f1, bins=500)
plt.axvline(x=np.median(f1), color='red', linestyle='--')

plt.text(20, 78000, 'Median: '+str(np.round(np.median(f1),2)))

plt.xlim(0,68)
plt.show()
