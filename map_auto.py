import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib.patches as mpatch
from matplotlib.patches import Rectangle

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()

sta_f = open('input/all_station_crossing_db.txt','r')

min_lon = -150.7
max_lon = -147.3
min_lat = 62.2
max_lat = 65.3
f = 1.0/np.cos(62*np.pi/180)

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
sta = seismo_data['Station']

#Loop through each station in text file that we already know comes within 2km of the nodes			
for line in sta_f.readlines():
	val = line.split(',')
	date = val[0]
	flight = val[1]
	station = val[5]
	flight_file = '/scratch/irseppi/nodal_data/flightradar24/' + date + '_positions/' + date + '_' + flight + '.csv'
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	time = flight_data['snapshot_id']
	speed = flight_data['speed']
	alt = flight_data['altitude']

	for l in range(len(time)):
		if str(time[l]) == str(val[2]):
			for t  in range(len(sta)):
				if sta[t] == station:
					# Plot the main figure
					fig, ax = plt.subplots()
					plt.gca().set_aspect(f)
					
					ax.scatter(seismo_longitudes, seismo_latitudes, c='red', s = 3, label='seismometers')
					ax.plot(flight_longitudes, flight_latitudes, '-o', c='c', lw=1, ms = 1, label='flight path')
					
					# Set labels and title
					ax.set_xlim(min_lon, max_lon)
					ax.set_ylim(min_lat, max_lat)
					ax.set_yticks(fontsize = 13)
					ax.set_xticks(fontsize = 13)

					# Define the position of the zoomed portion
					sub_axes = plt.axes([1.1, 0.5, 0.4, 0.4])

					# Plot the zoomed portion
					sub_axes.scatter(seismo_longitudes, seismo_latitudes, c='red')
					sub_axes.scatter(flight_longitudes, flight_latitudes, c='c')

					#Added a box to highligh zoomed in location
					y =[flight_latitudes[l],  seismo_latitudes[t]]
					x = [flight_longitudes[l], seismo_longitudes[t]]
					yy = sum(y)/len(y)
					xx = sum(x)/len(x)
					
					minl = xx - 0.05
					maxl = xx + 0.05
					minla = yy - 0.03
					maxla = yy + 0.03

					plt.gca().add_patch(Rectangle((minl, minla), 0.1, 0.06, ls="-", lw = 1, ec = 'k', fc="none", zorder=2.5))
					print(sta[t], flight)


					#Label station
					sub_axes.text(seismo_longitudes[t], seismo_latitudes[t], sta[t], fontsize=9, fontweight='bold')
					sub_axes.scatter(seismo_longitudes[t], seismo_latitudes[t], c='pink')
					
					ht = datetime.datetime.utcfromtimestamp(time[l])

					#Label timestamp 
					sub_axes.text(flight_longitudes[l], flight_latitudes[l], ht, fontsize=9, fontweight='bold')
					sub_axes.scatter(flight_longitudes[l], flight_latitudes[l], c='orange')
					
					
					sub_axes.plot(x,y, '--', c='orange')
					sub_axes.text(xx,yy, str(round(dist, 2))+'km', fontsize=8, fontweight='bold')
					
					# Set labels and title
					sub_axes.set_xlim(minl, maxl)
					sub_axes.set_ylim(minla, maxla)
					sub_axes.tick_params(axis='both', which='major', labelsize=9)
					plt.show()

					#Save
					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/pmap_z2/' + date + '/'+flight+'/'+station+'/'
					make_base_dir(BASE_DIR)
					plt.savefig('/scratch/irseppi/nodal_data/plane_info/pmap_z2/'+ date + '/'+flight+'/'+station+'/zmap_'+flight+'_' + str(time[l]) + '.png')
					plt.close()
						
				else:
					continue
				
			else:
				continue

