import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file2 = pd.read_csv('/home/irseppi/REPOSITORIES/parkshwynodal/input/all_station_crossing_db_C185.csv', sep=",")
tail_nums = file2['TAIL_NUM'].tolist()
flight = file2['FLIGHT_NUM'].tolist()

data = np.genfromtxt('old_new_invers.txt', delimiter=',', skip_header=0, usecols=(1, 3, 5, 6, 8, 10, 12), filling_values=np.nan)
flight_num = data[:, 0]
time_old = data[:, 1]
v0_old = data[:, 2]
distance_old = data[:, 3]
time_new = data[:, 4]
v0_new = data[:, 5]
peaks_new = np.array([float(peak) for peak in data[:, 6]])

color_dict = {}
count_dict = {}
peaks_dict_new = {}
y = 0

for i, flight_num in enumerate(flight_num):
    if flight_num not in color_dict:
        color_dict[flight_num] = np.random.rand(3,)
        peaks_dict_new[flight_num] = []
        count_dict[flight_num] = []
    peaks_dict_new[flight_num].extend([list(peaks_new)])  # Convert peaks_new[i] to a list
    counts = []
    for i in range(len(peaks_new)):
        counts.append(y)
    count_dict[flight_num].extend(counts)  # Create a list of the same length as peaks_new[i]
    y =+ 1
fig, ax1 = plt.subplots(1, 1, sharex=False, figsize=(8, 6))
ax1.margins(x=0)
ax2 = fig.add_axes([0.83, 0.11, 0.07, 0.77], sharey=ax1)
ax3 = fig.add_axes([0.90, 0.11, 0.04, 0.77], sharey=ax1)
ax4 = fig.add_axes([0.94, 0.11, 0.04, 0.77], sharey=ax1)
ax1.set_title('Frequency Peaks')

for tail_num, peaks in peaks_dict_new.items():
    color = color_dict[tail_num]
    y = count_dict[tail_num]

    peaks = np.array(peaks)  # Convert peaks to a NumPy array
    y = np.array(y)  # Convert y to a NumPy array


    ax1.scatter(peaks, y, c=color, label=tail_num)
    rpm = 60 * (peaks / 3)
    ax2.scatter(rpm, y, c=color)
    diff = np.diff(peaks)
    valid_diff_indices = np.array(np.where(np.logical_and(diff > 16, diff < 24))[0])


    valid_diff = diff[valid_diff_indices]
    ax3.scatter(valid_diff, y[valid_diff_indices], c=color)

    #ax4.scatter(rpm[valid_diff_indices[1]] / valid_diff, y[valid_diff_indices[1]], c=color)

ax2.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax3.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax4.tick_params(left=False, right=False, labelleft=False, labelbottom=True, bottom=True)
ax2.set_title('rpm')
ax3.set_title('\u0394' + 'F')
ax1.set_xlabel('Frequency')
ax2.set_xlabel('rpm')
ax3.set_xlabel('\u0394' + 'F')
ax1.legend(loc='upper left', fontsize='x-small')
ax1.set_xlim(0, 300)
ax1.set_xticks(range(0, 251, 25))

plt.show()