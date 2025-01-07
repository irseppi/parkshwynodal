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
	line_magnitude = line_vector[0] ** 2 + line_vector[1] ** 2
	projection_length_ratio = dot_product / line_magnitude
	return projection_length_ratio

#################################################################################################################################

def closest_projection(flight_latitudes, flight_longitudes, index, timestamp, seismo_latitude, seismo_longitude):
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
		float: The closest distance between the flight point and the seismic station and the time this occurrs.
	"""
	
	closest_distance = float('inf')
	closest_lat = flight_latitudes[index]
	closest_lon = flight_longitudes[index]
	timestamp1 = timestamp[index]

	timestamp2 = None
	for i in [index-1, index+1]:
		if i >= 0 and i < len(flight_latitudes):
			distance = calculate_distance(flight_latitudes[i], flight_longitudes[i], seismo_latitude, seismo_longitude)
			if distance < closest_distance:
				second_closest_lat = flight_latitudes[i]
				second_closest_lon = flight_longitudes[i]
				closest_distance = distance
				timestamp2 = timestamp[i]
		else: 
			second_closest_lat = flight_latitudes[index]
			second_closest_lon = flight_longitudes[index]
			timestamp2 = timestamp[index]

	lat_timestamp_dif_vec = second_closest_lat - closest_lat
	lon_timestamp_dif_vec = (second_closest_lon - closest_lon)
	lat_seismo_dif_vec = seismo_latitude - closest_lat
	lon_seismo_dif_vec = (seismo_longitude - closest_lon)

	line_vector = (lat_timestamp_dif_vec, lon_timestamp_dif_vec)
	station_vector = (lat_seismo_dif_vec, lon_seismo_dif_vec)

	projection_length_ratio = calculate_projection(line_vector, station_vector)

	closest_point_on_line_lat = closest_lat + projection_length_ratio * lat_timestamp_dif_vec
	closest_point_on_line_lon = closest_lon + projection_length_ratio * lon_timestamp_dif_vec
	closest_distance = calculate_distance(closest_point_on_line_lat, closest_point_on_line_lon, seismo_latitude, seismo_longitude)

	closest_time = timestamp1 + projection_length_ratio*(timestamp2 - timestamp1)

	return closest_point_on_line_lat, closest_point_on_line_lon, closest_distance, closest_time

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
		float: A tuple containing the latitude and longitute of the closest point between the flight and the seismic station, the distance between the closest point and the station, and the time the closest approach occurs at.

	"""
	clat = flight_latitudes[index]
	clon = flight_longitudes[index]
		
	closest_lat = clat
	closest_lon = clon
	dist_lim = 2.01

	for tr in range(0,2):
		if  flight_latitudes[index+1] < flight_latitudes[index-1]:
			if tr == 0:
				sclat = flight_latitudes[index-1]
				sclon = flight_longitudes[index-1]  

				x = [clon, sclon]
				y = [clat, sclat]
				m = (y[1]-y[0])/(x[1]-x[0])
				b = y[0] - m*x[0]
				
				for point in np.arange(clat, sclat, 0.000001):
					lat = point
					lon = (lat - b)/m
					dist_km = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000

					if dist_km < dist_lim:
						dist_lim = dist_km
						closest_lat = lat
						closest_lon = lon
						c2lat = flight_latitudes[index-1]
						c2lon = flight_longitudes[index-1]
						index2 = index - 1
			elif tr == 1:
				sclat = flight_latitudes[index+1]
				sclon = flight_longitudes[index+1]

				x = [clon, sclon]
				y = [clat, sclat]
				m = (y[0]-y[1])/(x[0]-x[1])
				b = y[0] - m*x[0]
				
				for point in np.arange(sclat, clat, 0.000001):
					lat = point
					lon = (lat - b)/m
					dist_km = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000

					if dist_km < dist_lim:
						dist_lim = dist_km
						closest_lat = lat
						closest_lon = lon
						c2lat = flight_latitudes[index+1]
						c2lon = flight_longitudes[index+1]
						index2 = index + 1
			
		elif flight_latitudes[index+1] > flight_latitudes[index-1]:
			if tr == 0:
				sclat = flight_latitudes[index-1]
				sclon = flight_longitudes[index-1]  

				x = [clon, sclon]
				y = [clat, sclat]
				m = (y[0]-y[1])/(x[0]-x[1])
				b = y[0] - m*x[0]

				for point in np.arange(sclat, clat, 0.000001):
					lat = point
					lon = (lat - b)/m
					dist_km = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000

					if dist_km < dist_lim:
						dist_lim = dist_km
						closest_lat = lat
						closest_lon = lon
						c2lat = flight_latitudes[index-1]
						c2lon = flight_longitudes[index-1]
						index2 = index - 1
			elif tr == 1:
				sclat = flight_latitudes[index+1]
				sclon = flight_longitudes[index+1]

				x = [clon, sclon]
				y = [clat, sclat]
				m = (y[1]-y[0])/(x[1]-x[0])
				b = y[0] - m*x[0]

				for point in np.arange(clat, sclat, 0.000001):
					lat = point
					lon = (lat - b)/m
					dist_km = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000

					if dist_km < dist_lim:
						dist_lim = dist_km
						closest_lat = lat
						closest_lon = lon
						c2lat = flight_latitudes[index+1]
						c2lon = flight_longitudes[index+1]
						index2 = index + 1
		else:
			continue
	
	if dist_lim < 2.01:
		'''
		for location in np.arange((closest_lon-0.000001),(closest_lon+0.000001),0.0000000001):
			lon = location
			lat = m*lon + b
			dist = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000
			if dist < dist_lim:
				dist_lim = dist
				closest_lon = lon
				closest_lat = lat
		if dist_lim < 1:
			for location in np.arange((closest_lon-0.0000000001),(closest_lon+0.0000000001),0.0000000000001):
				lon = location
				lat = m*lon + b
				dist = calculate_distance(lat, lon, seismo_latitude, seismo_longitude)/1000
				if dist < dist_lim:
					dist_lim = dist
					closest_lon = lon
					closest_lat = lat
		'''
		dist_old_new = calculate_distance(closest_lat, closest_lon, clat, clon)/1000
		dist_old_old = calculate_distance(c2lat, c2lon, clat, clon)/1000
		ratio = dist_old_new/dist_old_old
		timestamp = timestamp[index] + (timestamp[index] - timestamp[index2])*ratio
		
		return  closest_lat, closest_lon, dist_lim, timestamp
	else:
		return None, None, None, None

###################################################################################################################################

def calc_time(t0,dist,alt,c):
	"""
	Calculate the time at which the acoustic wave reaches the station.

	Parameters:
	t0 (float): Epoch time at which the wave is generated by the aircraft (in seconds).
	dist (float): Horizontal distance between the station and the aircraft at t0 (in meters).
	alt (float): Altitude of the aircraft at t0 (in meters).

	Returns:
	float: Time at which the acoustic wave reaches the station (in seconds).
	"""
	t = t0 + (np.sqrt(dist**2 + alt**2))/c 
	return t

#####################################################################################################################

def calc_f(f0, t, l, v0,c):
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

	tflight = t - (np.sqrt(t**2 - (1 - v0**2/c**2) * (t**2 - l**2/c**2))) / (1 - v0**2/c**2)
	f = f0 * (1 / (1 + (v0/c) * (v0 * tflight / (np.sqrt(l**2 + (v0 * tflight)**2)))))

	return f, tflight

############################################################################################################################

def calc_ft(times, tprime0, f0, v0, l, c):
	"""
	Calculate the frequency at each given time using the model parameters.

	Args:
		times (list): List of time values.
		tprime0 (float): The time at which the central frequency of the overtones occur, when the aircraft is at the closest approach to the station.
		f0 (float): Fundamental frequency produced by the aircraft.
		v0 (float): Velocity of the aircraft.
		l (float): Distance between the station and the aircraft at the closest approach.
		c (float): Speed of sound.

	Returns:
		list: List of calculated frequency values.
	"""
	ft = []
	for tprime in times:
		t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
		ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
								
		ft.append(ft0p)
	return ft

###################################################################################################################################################################

def calc_f0(tprime, tprime0, ft0p, v0, l, c):
	"""
	Calculate the fundamental frequency produced by an aircraft where the wave is generated given the model parameters.

	Parameters:
	tprime (float): Time at which an aribitrary frequency (ft0p) is observed on the station.
	tprime0 (float):  The time at which the central frequency of the overtones occur, when the aircraft is at the closest approach to the station.
	ft0p (float): Frequencyrecorded on the seismometer, picked from the overtone doppler curve.
	v0 (float): Velocity of the aircraft.
	l (float): Distance between the station and the aircraft at the closest approach.
	c (float): Speed of sound..

	Returns:
	f0 (float): Fundamental frequency produced by the aircraft. (Frequency at the source.) 
	"""
	t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
	f0 = ft0p*(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
	return f0

####################################################################################################################################################################################################################################################################################################################

def df(f0,v0,l,tp0,tp,c):   
    """
	Calculate the derivatives of f with respect to f0, v0, l, and tp0.

	Parameters:
	f0 (float): Fundamental frequency produced by the aircraft.
	v0 (float): Velocity of the aircraft.
	l (float): Distance of closest approach between the station and the aircraft.
	tp0 (float): Time of that the central frequency of the overtones occur, when the aircraft is at the closest approach to the station.
	tp (float): Array of times.
	Returns:
	tuple: A tuple containing the derivatives of f with respect to f0, v0, l, and tp0.
	"""

    #derivative with respect to f0
    f_derivef0 = (1 / (1 - (c * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))) /((c**2 - v0**2) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2))))


    #derivative of f with respect to v0
    f_derivev0 = (-f0 * v0 * (-2 * l**4 * v0**4 + l**2 * (tp - tp0)**2 * v0**6 + c**6 * (tp - tp0) * (2 * l**2 + (tp - tp0)**2 * v0**2) * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) + 
    c**2 * (4 * l**4 * v0**2 - (tp - tp0)**4 * v0**6 + l**2 * (tp - tp0) * v0**4 * (5 * tp - 5 * tp0 - 3 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))) - c**4 * 
    (2 * l**4 - 3 * (tp - tp0)**3 * v0**4 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4)) - l**2 * (tp - tp0) * v0**2 * (-6 * tp + 6 * tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
    (l**2 + (tp - tp0)**2 * v0**2))/c**4)))) / (c * (c - v0) * (c + v0) * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
    (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) - c**2 * np.sqrt(l**2 + (c**4 * v0**2 * 
    (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2))**2))
    


    #derivative of f with respect to l
    f_derivel = ((f0 * l * (tp - tp0) * (c - v0) * v0**2 * (c + v0) * ((-tp + tp0) * v0**2 + c**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))) / 
    (c * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + 
    (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4) - 
    c**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + 
    (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2)) / c**4))**2) / (c**2 - v0**2)**2))**2))


    #derivative of f with respect to tp0
    f_derivetprime0 = ((f0 * l**2 * (c - v0) * v0**2 * (c + v0) * ((-tp + tp0) * v0**2 + c**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))) / 
    (c * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * 
    (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) * (c * (-tp + tp0) * v0**2 + c * v0**2 * np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4) - 
    c**2 * np.sqrt(l**2 + (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2) + v0**2 * np.sqrt(l**2 + 
    (c**4 * v0**2 * (-tp + tp0 + np.sqrt((-l**2 * v0**2 + c**2 * (l**2 + (tp - tp0)**2 * v0**2))/c**4))**2)/(c**2 - v0**2)**2))**2))


    return f_derivef0, f_derivev0, f_derivel, f_derivetprime0

#####################################################################################################################################################################################################################################################################################################################

def invert_f(m0, coords_array, c, num_iterations,sigma = 1):
	"""
	Inverts the function f using the given initial parameters and data array.

	Args:
		m0 (numpy.ndarray): Initial parameters for the function f.
		coords_array (numpy.ndarray): Data picks along overtone doppler curve.
		num_iterations (int): Number of iterations to perform.

	Returns:
		numpy.ndarray: The inverted parameters for the function f.
	"""
	w,_ = coords_array.shape
	fobs = coords_array[:,1]
	tobs = coords_array[:,0]
	m = m0
	n = 0
	while n < num_iterations:
		fnew = []
		G = np.zeros((w,4)) #partial derivative matrix of f with respect to m
		#partial derivative matrix of f with respect to m 
		for i in range(0,w):
			f0 = m[0]
			v0 = m[1]
			l = m[2]
			tprime0 = m[3]
			tprime = tobs[i]
			t = ((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2)
			ft0p = f0/(1+(v0/c)*(v0*t)/(np.sqrt(l**2+(v0*t)**2)))
			f_derivef0, f_derivev0, f_derivel, f_derivetprime0 = df(m[0], m[1], m[2], m[3], tobs[i],c)
			
			G[i,0:4] = [f_derivef0, f_derivev0, f_derivel, f_derivetprime0]

			fnew.append(ft0p) 
		try:
			covmlsq = (sigma**2)*la.inv(G.T@G)
		except:
			covmlsq = (sigma**2)*la.pinv(G.T@G)
		try:
			m = np.reshape(np.reshape(m0,(4,1))+ np.reshape(la.inv(G.T@G)@G.T@(np.reshape(fobs, (len(coords_array), 1)) - np.reshape(np.array(fnew), (len(coords_array), 1))), (4,1)), (4,))
		except:
			m = np.reshape(np.reshape(m0,(4,1))+ np.reshape(la.pinv(G.T@G)@G.T@(np.reshape(fobs, (len(coords_array), 1)) - np.reshape(np.array(fnew), (len(coords_array), 1))), (4,1)), (4,))
		print(m)
		m0 = m
		n += 1

	return m, covmlsq

####################3####################################################################################################################################################################

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

