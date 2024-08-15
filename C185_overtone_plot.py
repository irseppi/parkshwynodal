import os
import matplotlib.pyplot as plt
import numpy as np

# Define the directory where your files are located
file = 'output/C185data.csv'
plt.figure()
# Initialize lists to store the data from the files
mode1_x = []
mode1_y = []
mode2_x = []
mode2_y = []

with open(file, 'r') as f:
    x = 0
    # Read the data from the file and append it to the respective lists based on the mode
    for line in f.readlines():
        lines = line.split(',')
        try:
            if x <= 28:
                peaks = np.array(lines[6])

            else:
                peaks = np.array(lines[7])

        except:
            continue
        peaks = str(peaks)  # Replace "string" with "str"
        peaks = np.array(peaks.split(' '))
        for peak in peaks:
            try:
                peak = float(peak)
            except:
                continue
            if abs(float(peak) - 62.0) <= 2:
                mode_x = 1
            else:
                mode_x = 2
               
        if mode_x == 1:
            for peak in peaks:
                try:
                    peak = float(peak)
                except:
                    continue
                mode1_x.append(peak)
                mode1_y.append(x)
        else:
           
            for peak in peaks:
                try:
                    peak = float(peak)
                except:
                    continue
                mode2_x.append(peak)
                mode2_y.append(x)


        x += 1

# Plot the dots for mode 1
plt.scatter(mode1_x, mode1_y, label='Mode 1')

# Plot the dots for mode 2
plt.scatter(mode2_x, mode2_y, label='Mode 2')

# Add legend
plt.legend()


plt.show()
