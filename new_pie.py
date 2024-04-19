import pandas as pd
import matplotlib.pyplot as plt

infile = open('input/all_station_crossing_db.txt', 'r')
outfile = open('output.csv', 'a')  # Open the file in append mode

equip_counts = {}  # Define the "equip_counts" dictionary before the loop
flight_nums = {}  # Define the "equip_counts" dictionary before the loop

for line in infile:
	if line.strip():  # Skip empty lines
		equip = line.split(',')[-1]  # Get the equipment type from the line
		flight_num = line.split(',')[1]
		if flight_num not in flight_nums:
			equip_counts[equip] = equip_counts.get(equip, 0) + 1  # Increment the count for the equipment type
			flight_nums[flight_num] = 1
		else:
			continue

equip_counts = {'Unkown': 296, 'AT73': 50, 'B190': 130, 'B738': 136, 'B737': 367, 'B739': 139, 'C185': 61, 'B77W': 67, 'PA31': 176, 'DH8A': 95, 'DHC6': 1, 'C172': 13, 'C208': 244, 'DH3T': 40, 'BE20': 37, 'E75S': 4, 'AS50': 4, 'B744': 9, 'SW4': 27, 'PC12': 52, 'B772': 33, 'B789': 30, 'DHC2': 12, 'B407': 12, 'R44': 5, 'A359': 10, 'B06': 2, 'C46': 2, 'B763': 9, 'GA8': 20, 'B732': 1, 'C182': 13, 'CH7B': 5, 'B788': 18, 'C441': 7, 'B18T': 4, 'C180': 8, 'PA34': 1, 'B77L': 6, 'B350': 1, 'PA18': 22, 'C206': 7, 'BE35': 1, 'C82S': 1, 'B733': 6, 'PA46': 3, 'A332': 1, 'PA32': 9, 'GLF5': 1, 'B748': 1, 'CRJ2': 1, 'AT8T': 3, 'BE10': 1, 'AC6L': 4, 'B412': 2, 'PA30': 2, 'BE58': 1, 'BE36': 1}

# Plotting the pie chart
labels = equip_counts.keys()
sizes = equip_counts.values()

plt.pie(sizes, labels=labels, autopct=str(sizes)) #, startangle=90)
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
plt.show()
