import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

#Open files in folder and load data to dataframe
dire = 'input/plane_overtones/'
for files in os.listdir(dire):
	try:
		f = os.path.join(dire, files)
		data = pd.read_csv(f)
		peaks = data.iloc[:, 1]
		amp = data.iloc[:, 2]
		pl = []
		for a in range(len(amp)):
			x = (amp[a]-np.min(amp))/(np.max(amp)-np.min(amp))
			pl.append(x)
		# Plot the data as a bar graph
		plt.figure()
		plt.hist(peaks, edgecolor="red", bins=50, weights=pl)
		plt.xlim(5,np.max(peaks+1))
		plt.xticks(np.arange(5, np.max(peaks), step=5))
		plt.xlabel('Freq [Hz]')
		plt.ylabel('Count')
		plt.title(files)
		plt.show()
		plt.close()
	except:
		continue
