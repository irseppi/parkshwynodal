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
mode3_x = []
mode3_y = []
mode1_f1 = []
mode1_f2 = []
mode1_f3 = []
mode1_f4 = []
mode1_f5 = []
mode1_f6 = []
mode1_f7 = []
mode1_f8 = []
mode1_f9 = []
mode1_f10 = []
mode2_f1 = []
mode2_f2 = []
mode2_f3 = []
mode2_f4 = []
mode2_f5 = []
mode2_f6 = []
mode2_f7 = []
mode2_f8 = []
mode2_f9 = []
with open(file, 'r') as f:
    x = 0
    xx = 0
    yy = 0
    zz = 0
    # Read the data from the file and append it to the respective lists based on the mode
    for line in f.readlines():
        lines = line.split(',')
        x += 1
        try:
            if x <= 29:
                peaks = np.array(lines[6])

            else:
                peaks = np.array(lines[7])

        except:
            continue
        peaks = str(peaks)  # Replace "string" with "str"
        peaks = np.array(peaks.split(' '))
        for peak in peaks:
            try:
                peak = float(peak[0:-1])
            except:
                continue
            if abs(float(peak) - 82.5) <= 1.5 or abs(float(peak) - 124.5) <= 1.5:
                xx += 1
                mode_x = 1
                break
            elif abs(float(peak) - 76) <= 5 or abs(float(peak) - 114) <= 5 or abs(float(peak) - 152) <= 5:
                mode_x = 2
                yy += 1
                break
            else:
                mode_x = 3
                zz += 1
                

        if mode_x == 1:
            for peak in peaks:
                try:
                    peak = float(peak[0:-1])
                except:
                    continue
                if abs(float(peak) - 62) <= 1.5:
                    mode1_f1.append(peak)
                elif abs(float(peak) - 82.5) <= 1.5:
                    mode1_f2.append(peak)
                elif abs(float(peak) - 103) <= 2:
                    mode1_f3.append(peak)
                elif abs(float(peak) - 124.5) <= 1.5:
                    mode1_f4.append(peak)
                elif abs(float(peak) - 145) <= 2:
                    mode1_f5.append(peak)
                elif abs(float(peak) - 166) <= 2:
                    mode1_f6.append(peak)
                elif abs(float(peak) - 186.5) <= 2:
                    mode1_f7.append(peak)
                elif abs(float(peak) - 208) <= 2:
                    mode1_f8.append(peak)
                elif abs(float(peak) - 228.5) <= 2:
                    mode1_f9.append(peak)
                elif abs(float(peak) - 249) <= 4:
                    mode1_f10.append(peak)
                mode1_x.append(peak)
                mode1_y.append(xx)
        elif mode_x == 2:
            for peak in peaks:
                try:
                    
                    peak = float(peak[0:-1])
                    print(peak)
                except:
                    continue
                mode2_x.append(peak)
                mode2_y.append(yy)
                print(peak)
                if abs(float(peak) - 28) <= 4:
                    print('here')
                    mode2_f1.append(peak)
                elif abs(float(peak) - 39) <= 4:
                    mode2_f2.append(peak)
                elif abs(float(peak) - 58) <= 4:
                    mode2_f3.append(peak)
                elif abs(float(peak) - 77) <= 4:
                    mode2_f4.append(peak)
                elif abs(float(peak) - 96.5) <= 4:
                    mode2_f5.append(peak)
                elif abs(float(peak) - 116) <= 5:
                    mode2_f6.append(peak)
                elif abs(float(peak) - 133) <= 6:
                    mode2_f7.append(peak)
                elif abs(float(peak) - 153.5) <= 4:
                    mode2_f8.append(peak)
                elif abs(float(peak) - 173) <= 3:
                    mode2_f9.append(peak)
        else:
            for peak in peaks:
                try:
                    peak = float(peak)
                except:
                    continue
                mode3_x.append(peak)
                mode3_y.append(zz)

extra = [[121.3,175.54,236.3,59.48],[120.54,179.63,217.25,236.04,59.02],[116.91,172.8,232.69,57.7],[114.21,170.63,227.46]]
for i in range(len(extra)):
    x = 0
    xx = 0
    yy = 0
    zz = 0
    for peak in peaks:
        try:
            peak = float(peak[0:-1])
        except:
            continue
        if abs(float(peak) - 82.5) <= 1.5 or abs(float(peak) - 124.5) <= 1.5:
            xx += 1
            mode_x = 1
            break
        elif abs(float(peak) - 76) <= 5 or abs(float(peak) - 114) <= 6 or abs(float(peak) - 152) <= 5:
            mode_x = 2
            yy += 1
            break
        else:
            mode_x = 3
            zz += 1
                
        if mode_x == 1:
            for peak in peaks:
                try:
                    peak = float(peak[0:-1])
                except:
                    continue
                if abs(float(peak) - 62) <= 1.5:
                    mode1_f1.append(peak)
                elif abs(float(peak) - 82.5) <= 1.5:
                    mode1_f2.append(peak)
                elif abs(float(peak) - 103) <= 2:
                    mode1_f3.append(peak)
                elif abs(float(peak) - 124.5) <= 1.5:
                    mode1_f4.append(peak)
                elif abs(float(peak) - 145) <= 2:
                    mode1_f5.append(peak)
                elif abs(float(peak) - 166) <= 2:
                    mode1_f6.append(peak)
                elif abs(float(peak) - 186.5) <= 2:
                    mode1_f7.append(peak)
                elif abs(float(peak) - 208) <= 2:
                    mode1_f8.append(peak)
                elif abs(float(peak) - 228.5) <= 2:
                    mode1_f9.append(peak)
                elif abs(float(peak) - 249) <= 4:
                    mode1_f10.append(peak)
                mode1_x.append(peak)
                mode1_y.append(xx)
        elif mode_x == 2:
            for peak in peaks:
                try:
                    peak = float(peak[0:-1])
                except:
                    continue
                print(peak)
                if abs(float(peak) - 28) <= 4:

                    mode2_f1.append(peak)
                elif abs(float(peak) - 39) <= 4:
                    mode2_f2.append(peak)
                elif abs(float(peak) - 58) <= 4:
                    mode2_f3.append(peak)
                elif abs(float(peak) - 77) <= 4:
                    mode2_f4.append(peak)
                elif abs(float(peak) - 96.5) <= 4:
                    mode2_f5.append(peak)
                elif abs(float(peak) - 116) <= 5:
                    mode2_f6.append(peak)
                elif abs(float(peak) - 133) <= 6:
                    mode2_f7.append(peak)
                elif abs(float(peak) - 153.5) <= 4:
                    mode2_f8.append(peak)
                elif abs(float(peak) - 172) <= 6:
                    mode2_f9.append(peak)
                mode2_x.append(peak)
                mode2_y.append(yy)
        else:
            for peak in peaks:
                try:
                    peak = float(peak)
                except:
                    continue
                mode3_x.append(peak)
                mode3_y.append(zz)


medians_mode1 = [np.median(mode1_f1), np.median(mode1_f2), np.median(mode1_f3), np.median(mode1_f4), np.median(mode1_f5), np.median(mode1_f6), np.median(mode1_f7), np.median(mode1_f8), np.median(mode1_f9), np.median(mode1_f10)]
medians_mode2 = [np.median(mode2_f1), np.median(mode2_f2), np.median(mode2_f3), np.median(mode2_f4), np.median(mode2_f5), np.median(mode2_f6), np.median(mode2_f7), np.median(mode2_f8), np.median(mode2_f9)]
medians_mode1_rounded = [round(median, 1) for median in medians_mode1]

medians_mode2_rounded = [round(median, 1) for median in medians_mode2]

for medians in medians_mode1:
    plt.axvline(x=medians, color='orange', linestyle='--')
for medians in medians_mode2:
    plt.axvline(x=medians, color='blue', linestyle='--')

# Set the x-axis labels as rounded median values
plt.xticks(medians_mode1_rounded + medians_mode2_rounded)
# Add labels for the differences between medians
for i in range(len(medians_mode1_rounded) - 1):
    diff = round(medians_mode1_rounded[i+1] - medians_mode1_rounded[i], 1)
    plt.text(medians_mode1_rounded[i]+(diff/2), 1, f'{diff}', ha='center', va='bottom')
    plt.annotate('', xy=(medians_mode1_rounded[i], 1), xytext=(medians_mode1_rounded[i+1], 1),
                 arrowprops=dict(arrowstyle='<->', color='orange'))
    #diff_mode2_mode1 = round(medians_mode1_rounded[i] - medians_mode2_rounded[i], 1)
    #plt.text(medians_mode1_rounded[i]+diff_mode2_mode1, np.max(mode1_y) + 2, f'{diff_mode2_mode1}', ha='center', va='bottom')
for i in range(len(medians_mode2_rounded) - 1):
    diff = round(medians_mode2_rounded[i+1] - medians_mode2_rounded[i], 1)
    plt.text(medians_mode2_rounded[i]+(diff/2), 3, f'{diff}', ha='center', va='bottom')
    plt.annotate('', xy=(medians_mode2_rounded[i], 3), xytext=(medians_mode2_rounded[i+1], 3),
                 arrowprops=dict(arrowstyle='<->', color='blue'))
for i in range(0,7):
    diff = round(medians_mode1_rounded[i] - medians_mode2_rounded[i+2], 1)
    plt.text(medians_mode1_rounded[i]-(diff/2), 5, f'{diff}', ha='center', va='bottom')
    plt.annotate('', xy=(medians_mode1_rounded[i], 5), xytext=(medians_mode2_rounded[i+2], 5),
                 arrowprops=dict(arrowstyle='<->', color='green'))
# Plot the dots for modes
plt.scatter(mode1_x, mode1_y, label='Mode 1', color='orange')
mode2_y_modified = [(np.max(mode1_y)) + y for y in mode2_y]
plt.scatter(mode2_x, mode2_y_modified, label='Mode 2', color='blue')

# Add legend
plt.legend()


plt.show()
