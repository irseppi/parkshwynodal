import sys
import fileinput
import os
import pandas as pd
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from numpy.linalg import inv
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

#############################################################################################################################################################

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
			
def extract_flight(equipment):
	"""
	Extracts all rows from the 'all_station_crossing_db.txt' file whith the designated equipment type 
	and prints them into an infividual file labeled with the equipment type.
	
	Args:
		equipment (str): The equipment type to extract from the file.

	Returns:
		None

	"""

	input = open('input/all_station_crossing_db.txt','r')
	output = open('input/all_station_crossing_db_'+str(equipment)+'.txt','w')

	for line in input.readlines():
		val = line.split(',')
		if str(val[7][0:4]) == str(equipment):
			
			output.write(line)
		
	input.close()
	output.close()	

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

