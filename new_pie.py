import pandas as pd
import matplotlib.pyplot as plt
from prelude import load_flights, dist_less

# Load the seismometer location data
seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
sta = seismo_data['Station']

flight_files, filenames = load_flights(2, 4, 11, 27)

# Initialize an empty list to store DataFrames
eq_list = []
outfile = open('output.csv', 'w')

for i, flight_file in enumerate(flight_files):
	flight_data = pd.read_csv(flight_file, sep=",")
	flight_latitudes = flight_data['latitude']
	flight_longitudes = flight_data['longitude']
	fname = filenames[i]
	flight_num = fname[9:18]
	date = fname[0:8]
	con = dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes)

	if con == True:
		equip_info = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/' + date + '_flights.csv', sep=",")
		equip = equip_info['equip']
		flight_id = equip_info['flight_id']
		for j, fid in enumerate(flight_id):
			if int(fid) == int(flight_num):
				print('here')
				equip = str(equip[j])
				
				# Write the output to the file
				outfile.write(equip+'\n')
					
				break
			else:
				continue
outfile.close()

# Concatenate all the dataframes in the list
#df = pd.DataFrame(eq_list)

# Extract the 'equipment type' column and count the occurrences of each type
#equipment_counts = df['equip'].value_counts()
#equipment_percent = equipment_counts / equipment_counts.sum() * 100

# Group equipment types that occur less than 0.5% of the time into 'Other'
#other_count = equipment_counts[equipment_percent < 0.5].sum()
#equipment_counts = equipment_counts[equipment_percent >= 0.5]

# Ensure 'Other' category always appears in the pie chart
#if 'Other' not in equipment_counts.index:
#	equipment_counts = pd.concat([equipment_counts, pd.Series([0], index=['Other'])])

# Add the count of 'Other' equipment types
#equipment_counts['Other'] += other_count

# Plot the data as a pie chart
#equipment_counts.plot(kind='pie', autopct=lambda pct: f"{pct:.1f}% ({int(pct/100*equipment_counts.sum())})")
#plt.title('Occurrences of Equipment Types')
#plt.show()
