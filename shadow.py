import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#file = open('1o_atmc_v_2c.txt','r')
#file = open('temp_correction_1ovfull.txt','r')
file = open('2c_1o_v_full.txt','r')
file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM']
flight = file2['FLIGHT_NUM']

date_old = []
date_new = []
date_all = []

y=0
time_new = []
time_old = []

v0_old = []
v0_new = []

pppp_old = []
pppp_new = []

distance_old = []
distance_new = []

med_old = []
med_new = []

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

    f1 = []
    for i in range(len(peaks_old)):
        peak_old = peaks_old[i]
        pppp_old.append(float(peak_old))
        date_old.append(y)
        if i == 0:
            continue
        diff = float(peaks_old[i]) - float(peaks_old[i-1])
        if diff > 22 or diff < 18:
            continue
        f1.append(diff)
    med_old.append(np.nanmedian(f1))
    f2 = []
    for i in range(len(peaks_new)):
        peak_new = peaks_new[i]
        pppp_new.append(float(peak_new))
        date_new.append(y)
        if i == 0:
            continue
        diff = float(peaks_new[i]) - float(peaks_new[i-1])
        if diff > 22 or diff < 18:
            continue
        f2.append(diff)
    med_new.append(np.nanmedian(f2))
    date_all.append(y)

fig,ax1 = plt.subplots(1, 1, sharex=False, figsize=(8,6))     

ax1.margins(x=0)
ax2 = fig.add_axes([0.83, 0.11, 0.07, 0.77], sharey=ax1)

#ax1.scatter(pppp_old, date_old, c='b', label='-2C')
#ax1.scatter(pppp_new, date_new,c='r',label='Temperature Corrected')
ax1.scatter(pppp_old, date_old, c='b', label='One Overtone')
ax1.scatter(pppp_new, date_new, c='r',label='Full')

ax2.scatter(med_old, date_all, c='b')
ax2.scatter(med_new, date_all, c='r')

ax1.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
#ax1.set_title('Diffrence in Frequency outcome between Temperature Correction vs Constant -2 degrees C (One Overtone inversion)')
ax1.set_title('Diffrence in Frequency between One Overtone Inversion and Full Inversion (-2C)')
ax1.set_xlabel('Frequency')
ax2.set_xlabel('\u0394'+'F')
ax1.legend(loc='upper left',fontsize = 'x-small')
ax1.set_xlim(0, 300)
ax1.set_xticks(range(0, 251, 25)) 
plt.legend()
plt.show()
#
