import sys
import fileinput
from datetime import datetime
import os
import pandas as pd
import numpy as np
from pathlib import Path
from obspy.geodetics import gps2dist_azimuth

###############################################################

def make_base_dir(base_dir):
	"""
	Create a directory and its parent directories if they don't exist.

	Args:
		base_dir (str): The path of the directory to be created.

	Returns:
		None
	"""
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path / parent
			if not current_path.exists():
				current_path.mkdir()

################################################################

def distance(lat1, lon1, lat2, lon2):
	"""
	Calculate the distance in kilometers between two sets of latitude and longitude coordinates.

	Parameters:
	lat1 (float): Latitude of the first point.
	lon1 (float): Longitude of the first point.
	lat2 (float): Latitude of the second point.
	lon2 (float): Longitude of the second point.

	Returns:
	float: The distance in kilometers between the two points.
	"""
	dist = gps2dist_azimuth(lat1, lon1, lat2, lon2)
	dist_km = dist[0]/1000

	return dist_km

#################################################################################################################################

def dist_less(flight_latitudes, flight_longitudes, seismo_latitudes, seismo_longitudes):
	"""
	Check if the distance between any flight location and any seismic location is less than or equal to 2.

	Args:
		flight_latitudes (list): List of flight latitudes.
		flight_longitudes (list): List of flight longitudes.
		seismo_latitudes (list): List of seismic latitudes.
		seismo_longitudes (list): List of seismic longitudes.

	Returns:
		bool: True if the distance is less than or equal to 2, False otherwise.
	"""
	f = False
	for s in range(len(flight_latitudes)):
		for l in range(len(seismo_latitudes)):
			dist = distance(seismo_latitudes[l], seismo_longitudes[l], flight_latitudes[s], flight_longitudes[s])
			if dist <= 2:
				f = True
				break
			else:
				continue
	return f

#################################################################################################################################

def calculate_distance(lat1, lon1, lat2, lon2):
	"""
	Calculate the distance between two GPS coordinates.

	Args:
		lat1 (float): Latitude of the first coordinate.
		lon1 (float): Longitude of the first coordinate.
		lat2 (float): Latitude of the second coordinate.
		lon2 (float): Longitude of the second coordinate.

	Returns:
		float: The distance between the two coordinates in meters.
	"""
	distance, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)  # distance in meters
	return distance

#################################################################################################################################

def calculate_projection(line_vector, station_vector):
	"""
	Calculates the projection length ratio of a station vector onto a line vector.

	Args:
		line_vector (list): The line vector represented as a list of two elements.
		station_vector (list): The station vector represented as a list of two elements.

	Returns:
		float: The projection length ratio.

	"""
	dot_product = line_vector[0] * station_vector[0] + line_vector[1] * station_vector[1]
	line_magnitude_squared = line_vector[0] ** 2 + line_vector[1] ** 2
	projection_length_ratio = dot_product / line_magnitude_squared
	return projection_length_ratio

#################################################################################################################################

def closest_encounter(flight_latitudes, flight_longitudes, index, timestamp, seismo_latitude, seismo_longitude):
	"""
	Calculates the closest distance between a flight point and a seismic station.

	Args:
		flight_latitudes (list): List of latitude values for flight points.
		flight_longitudes (list): List of longitude values for flight points.
		index (int): Index of the flight point to calculate the closest distance from.
		timestamp: Timestamp of the flight point.
		seismo_latitude (float): Latitude of the seismic station.
		seismo_longitude (float): Longitude of the seismic station.

	Returns:
		float: The closest distance between the flight point and the seismic station.
	"""

	closest_distance = float('inf')
	closest_lat = flight_latitudes[index]
	closest_lon = flight_longitudes[index]

	for i in [index-1, index+1]:
		if i >= 0 and i < len(flight_latitudes):
			distance = calculate_distance(flight_latitudes[i], flight_longitudes[i], seismo_latitude, seismo_longitude)
			if distance < closest_distance:
				closest_distance = distance
				closest_lat = flight_latitudes[i]
				closest_lon = flight_longitudes[i]

	line_vector = (flight_latitudes[index] - closest_lat, flight_longitudes[index] - closest_lon)
	station_vector = (seismo_latitude - closest_lat, seismo_longitude - closest_lon)

	projection_length_ratio = calculate_projection(line_vector, station_vector)

	closest_point_on_line_lat = closest_lat + projection_length_ratio * line_vector[0]
	closest_point_on_line_lon = closest_lon + projection_length_ratio * line_vector[1]

	closest_distance = calculate_distance(closest_point_on_line_lat, closest_point_on_line_lon, seismo_latitude, seismo_longitude)

	return closest_distance

###################################################################################################################################

def calc_time(t, l, vo):
	"""
	Calculate the time at which the acoustic wave reaches the station.

	Parameters:
	t (float): Epoch time at which the wave is generated by the aircraft (in seconds).
	l (float): Shortest distance between the station and the aircraft (in meters).
	vo (float): Velocity of the aircraft (in meters per second).

	Returns:
	float: Time at which the acoustic wave reaches the station (in seconds).
	"""
	c = 343  # Speed of acoustic wave propagation (in meters per second)
	to = t + (np.sqrt(l**2 + (vo*t)**2)) / c
	return to

####################################################################################

def calc_f(f0, t, l, v0):
	"""
	Calculates the frequency shift and time of flight for a wave generated at an aircraft and received at a station.

	Args:
		f0 (float): The initial frequency of the wave generated at the aircraft.
		t (float): The epoch time at which the wave arrives at the station (in seconds).
		l (float): The distance of closest approach between the station and the aircraft (in meters).
		v0 (float): The velocity of the aircraft (in meters per second).

	Returns:
		tuple: A tuple containing the frequency shift (f) and the time of flight (tflight).
	"""
	c = 343  # Speed of acoustic wave propagation (in m/sec)
	tflight = t - (np.sqrt(t**2 - (1 - v0**2/c**2) * (t**2 - l**2/c**2))) / (1 - v0**2/c**2)
	f = f0 * (1 / (1 + (v0/c) * (v0 * tflight / (np.sqrt(l**2 + (v0 * tflight)**2)))))

	return f, tflight

#########################################################################################################################################################


def load_flights(month1, month2, first_day, last_day):
	"""
	Load flight files based on the specified months and days.

	Args:
		month1 (int): The starting month.
		month2 (int): The ending month.
		first_day (int): The first day of the range.
		last_day (int): The last day of the range.

		for only Feb use month1 = 2 and month2 = 3
		for only March use month1 = 3 and month2 = 4
		for Fab and March use month1 = 2 and month2 = 4
		for entire deployment use month1 = 2,first_day = 11,month2 = 4, and last_day = 27

	Returns:
		tuple: A tuple containing two lists - flight_files and filenames.
			   flight_files: A list of file paths for the flight files.
			   filenames: A list of filenames for the flight files.
	"""
	flight_files = []
	filenames = []

	for month in range(month1, month2):
		if month1 == 2 and month2 == 4:
			if month == 2:
				month = '02'
				for day in range(first_day, 29):
					day = str(day)
					directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
					for filename in os.listdir(directory):
						filenames.append(filename)
						f = os.path.join(directory, filename)
						if os.path.isfile(f):
							flight_files.append(f)
			elif month == 3:
				month = '03'
				for day in range(1, last_day):
					if day < 10:
						day = '0' + str(day)
						directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
						for filename in os.listdir(directory):
							filenames.append(filename)
							f = os.path.join(directory, filename)
							if os.path.isfile(f):
								flight_files.append(f)
					else:
						day = str(day)
						directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
						for filename in os.listdir(directory):
							filenames.append(filename)
							f = os.path.join(directory, filename)
							if os.path.isfile(f):
								flight_files.append(f)
		elif month1 == 2 and month2 == 3:
			month = '02'
			for day in range(first_day, last_day):
				day = str(day)
				directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
				for filename in os.listdir(directory):
					filenames.append(filename)
					f = os.path.join(directory, filename)
					if os.path.isfile(f):
						flight_files.append(f)
		elif month1 == 2 and month2 == 4:
			month = '03'
			for day in range(first_day, last_day):
				if day < 10:
					day = '0' + str(day)
					directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
					for filename in os.listdir(directory):
						filenames.append(filename)
						f = os.path.join(directory, filename)
						if os.path.isfile(f):
							flight_files.append(f)
				else:
					day = str(day)
					directory = '/scratch/irseppi/nodal_data/flightradar24/2019' + month + day + '_positions'
					for filename in os.listdir(directory):
						filenames.append(filename)
						f = os.path.join(directory, filename)
						if os.path.isfile(f):
							flight_files.append(f)
	return flight_files, filenames

####################################################

def modify_content(content):
	"""
	Function to modify the content of the file by making all the letters uppercase. 

	Args:
		content (str): The content to be modified.

	Returns:
		str: The modified content.
	"""
	return content.upper()

#############################################################

def modify_file(input_file_name, output_file_name):
	"""
	Makes all the letters in the input file and writes this modified content to the output file.

	Args:
		input_file_name (str): The path of the input file.
		output_file_name (str): The path of the output file.

	Returns:
		None
	"""
	# Read the input file
	with open(input_file_name, 'r') as file:
		content = file.read()

	# Modify the content
	modified_content = modify_content(content)

	# Write the modified content to the output file
	with open(output_file_name, 'w') as file:
		file.write(modified_content)

####################################################################################################################

def station_subset(filename, steps, outputname):
	"""
	Subset stations from the input file based on the given steps and save the result to the output file.

	Parameters:
	- filename (str): The path of the input file.
	- steps (int): The amount of steps to subset stations by. (ie. every 4th station is in the subset)
	- outputname (str): The filename for the output file.
	"""
	
	with open(filename) as handle:
		for lineno, line in enumerate(handle):
			if lineno % steps == 0:
				print(line)

#######################################################################################


def replace(filename, old_string, new_string):
	"""
	Replace all occurrences of 'old_string' with 'new_string' in the given file.

	Args:
		filename (str): The path of the file to be modified.
		old_string (str): The string to be replaced.
		new_string (str): The string to replace the old_string with.
	"""
	for i, line in enumerate(fileinput.input(filename, inplace=1)):
		sys.stdout.write(line.replace(old_string, new_string))

	# Example usage:
	# replace('filename.site', '', ' "')
	# or replace('#', '\n#')

##############################################################################################


def round_replace(filename, col_2round, precision, m2km):
	"""
	Replace the values in a specific column of a text file with rounded values.

	Args:
		filename (str): The path of the text file.
		col_2round (int): The column index to round.
		precision (int): The number of digits to round to.
		m2km (bool): Determines whether to convert the rounded value to kilometers.
			if m2km = 0 - replace number with rounded number
			if m2km = 1 - replace number with rounded number converted to km
	
	Returns:
		None
	"""

	col_2round = int(col_2round)
	precision = int(precision)

	for i, line in enumerate(fileinput.input(filename, inplace=1)):
		val = line.split()

		if m2km == False:
			new_val = round(float(val[col_2round]), precision)
			sys.stdout.write(line.replace(str(val[col_2round]), str(new_val)))

		if m2km == True:
			new_val = round(float(val[col_2round]) / 1000, precision)
			sys.stdout.write(line.replace(str(val[col_2round]), str(new_val)))

#########################################################################################################

def rename_file(flight_name):
	"""
	Renames all files in the specified flight collection so that the flight name is appended
	to the beginning of the file name, instead of a folder label.

	Args:
		flight_name (str): The name of the flight collection.

	Returns:
		None
	"""
	os.getcwd()
	collection = flight_name + '/'
	for i, filename in enumerate(os.listdir(collection)):
		for p, fil in enumerate(os.listdir(collection+filename)):
			os.rename(collection + filename + '/' + fil, collection + filename +'_'+ fil)

########################################################################################################

def extract_col(input_file, output_file, col, split_str):
	"""
	Extracts a specific column from a text file and writes it to another file.

	Args:
		input_file (str): The path to the input text file.
		output_file (str): The path to the output file where the extracted column will be written.
		col (int): The index of the column to extract (0-based index).
		split_str (str): The string that splits the text file into columns.

	Returns:
		None
	"""
	i = int(col)
	with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
		for line in f_in:
			line = line.split(split_str)
			print(line[i])
			f_out.write(line[i])

#########################################################################################

def date_round(input_file, output_file):
	"""
	Rounds the seconds of each timestamp in the input file to remove the milliseconds
	and writes the rounded timestamps to the output file.

	Args:
		input_file (str): The path to the input file containing timestamps.
		output_file (str): The path to the output file where rounded timestamps will be written.

	Returns:
		None
	"""
	# Remove the milliseconds from the timestamp
	with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
		for line in f_in:
			dt = datetime.strptime(line.strip(), '%Y-%m-%dT%H:%M:%S.%f')
			rounded_dt = dt.replace(second=round(dt.second))
			f_out.write(rounded_dt.strftime('%Y-%m-%d %H:%M:%S') + '\n')

########################################################################################

def count_flight(input_file, col_f, output_file, designator): 
	"""
	Counts the number of flights with a specific aircraft type designator in the input file.

	Args:
		input_file (str): Path to the input file.
		col_f (int): Column index where the flight equipment is located.
		output_file (str): Path to the output file.
		designator (str): Aircraft type designator to search for.

	Returns:
		None
	"""
	text = open(input_file, 'r')
	f = open(output_file, 'w')
	i = int(col_f)

	flight_data = pd.read_csv('20231010_Aircraft_UA_Fairbanks.csv', sep=",")
	eq = flight_data['TypeDesignator']
	des = flight_data['Description']

	count = 0
	for line in text.readlines():
		val = line.split(',')
		equip = val[i]
		for l in range(len(eq)):
			if str(eq[l]) == str(equip[0:4]) and str(des[l]) == designator:
				count = count + 1
				f.write(eq[l]+'\n')
	f.write(str(count))
	f.close()

#######################################################################################

def print_eq():
	"""
	Prints the aircraft type designator for each line in the 'all_station_crossing_db.txt' file.
	"""
	text = open('all_station_crossing_db.txt', 'r')

	for line in text.readlines():
		val = line.split(',')
		equip = val[6]
		print(equip)

