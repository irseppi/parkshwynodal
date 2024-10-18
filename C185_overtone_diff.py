import matplotlib.pyplot as plt
import numpy as np

file = open('output/C185data_updated.txt', 'r')
f1 = []
for line in file.readlines():
    lines = line.split(',')
    quality_num = int(lines[9])

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
    print(peaks_filt)
    for p in range(len(peaks_filt)):

        if p == 0:
            continue
        else:
            diff = float(peaks_filt[p]) - float(peaks_filt[p-1])
            print(float(peaks_filt[p]), float(peaks_filt[p-1]))
            print(diff)
            '''
            if 10 > diff:
                continue
            elif diff > 30:
                continue
                #missing = diff / 20
                #for m in range(int(missing)):
                #    f1.append()
            
            else:
            '''
            f1.append(diff)


plt.figure()
plt.hist(f1, bins=500)
plt.axvline(x=np.median(f1), color='red', linestyle='--')
plt.axvline(x=np.mean(f1), color='red', linestyle='--')

plt.text(20, 30, 'Median: '+str(np.round(np.median(f1),2)))
plt.text(20, 20, 'Mean: '+str(np.round(np.mean(f1),2)))

#plt.xlim(0,140)
plt.show()
