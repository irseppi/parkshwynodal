import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import obspy
import datetime
import pyproj
from prelude import *
from scipy.signal import find_peaks, spectrogram
from pathlib import Path


sta_f = open('output/C185data_updated.txt','r')
C185_output = open('output/C185data_1overtone.csv', 'w')

# Loop through each station in text file that we already know comes within 2km of the nodes
for line in sta_f.readlines():
    val = line.split(',')
    date = val[0]
    flight = val[1]
    sta = val[2]
    closest_time = val[3]
    tprime0 = float(val[4])
    quality_num = int(val[9])
    v0 = float(val[5])
    l = float(val[6])
    c = 343 
    output2 = '/home/irseppi/REPOSITORIES/parkshwynodal/output/C185_data_picks/overtonepicks/2019-'+str(date[4:6])+'-'+str(date[6:8])+'/'+str(flight)+'/'+str(sta)+'/'+str(closest_time)+'_'+str(flight)+'.csv'

    f0_array = []
    #if Path(output2).exists():
    #    print('File exists')
    peaks = []
    freqpeak = []
    with open(output2, 'r') as file:
        for line in file:
            # Split the line using commas
            pick_data = line.split(',')
            tprime = float(pick_data[0])
            ft0p = float(pick_data[1])
            f0 = calc_f0(tprime, tprime0, ft0p, v0, l, c)
            f0_array.append(f0)
    file.close()  # Close the file after reading
    #else:
    #ontinue
    f0_array = np.array(f0_array)
    #f0_array = sorted(f0_array)
    print(f0_array)
    C185_output.write(str(date)+','+str(flight)+','+str(sta)+','+str(closest_time)+','+str(tprime0)+','+str(v0)+','+str(l)+','+str(f0_array)+',\n')

C185_output.close()
