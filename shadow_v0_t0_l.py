import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = open('old_new_invers.txt','r')

date = []
y=0
time_new = []
time_old = []

v0_old = []
v0_new = []

distance_old = []
distance_new = []


# Iterate over each line in the file
for line in file.readlines():
    y += 1
    counts = []

    lines = line.split(',')
    time_old.append(float(lines[4]))
    time_new.append(float(lines[9]))
    v0_old.append(float(lines[5]))
    v0_new.append(float(lines[10]))
    flight_num = float(lines[1])
    distance_old.append(float(lines[6]))
    distance_new.append(float(lines[11]))
    date.append(y)

plt.figure()
plt.scatter(v0_old, date, c='b')
plt.scatter(v0_new, date, c='r')
plt.title('v0')
plt.show()

plt.figure()
plt.scatter(np.abs(np.array(v0_new) - np.array(v0_old)), date, c='b') 
plt.show()

plt.figure()
plt.scatter(distance_old, date, c='b') 
plt.scatter(distance_new, date, c='r')
plt.title('distance')
plt.show()

plt.figure()
plt.scatter(np.abs(np.array(distance_new) - np.array(distance_old)), date, c='b') 
plt.show()

plt.figure()
plt.scatter(time_old, date, c='b') 
plt.scatter(time_new, date, c='r')
plt.title('time')
plt.show()


plt.figure()
plt.scatter(np.abs(np.array(time_new) - np.array(time_old)), date, c='b') 
plt.title('time')
plt.show()