%matplotlib inline
import faulthandler
import numpy as np
import matplotlib.pyplot as plt
import obspy 

from obspy.signal.filter import bandpass
from prelude import *

faulthandler.enable()
tr = obspy.read("/home/irseppi/nodal_data/500sps/2019_02_11/ZE_1001_DPZ.msd")
tr_new = obspy.read("/home/irseppi/nodal_data/50sps/2019_02_11/ZE_1001_DPZ.msd")

fig1, ax1 = plt.subplots()
fig1.set_figwidth(100.0)
fig1.set_figheight(30.0)
ax1.set_xlim([0, 3500])
tr.spectrogram(axes = ax1,dbscale=True)
plt.savefig("assets/example500spec.png")
plt.show()

fig2, ax2 = plt.subplots()
fig2.set_figwidth(100.0)
fig2.set_figheight(30.0)
ax2.set_xlim([0, 3500])
tr_new.spectrogram(axes = ax2,dbscale=True)
plt.savefig("assets/example50spec.png")
plt.show()
