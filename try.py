import datetime
from math import radians, sin, cos, sqrt, atan2
from geopy.distance import gps2dist_azimuth
from math import radians, cos
def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0  # radius of the Earth in kilometers

    lat1_rad = radians(lat1)
    lon1_rad = radians(lon1)
    lat2_rad = radians(lat2)
    lon2_rad = radians(lon2)

    dlon = lon2_rad - lon1_rad
    dlat = lat2_rad - lat1_rad

    a = sin(dlat / 2)**2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    return R * c  # returns distance in kilometers

    return wave_arrival_time

def closest_encounter(aircraft_data, seismometer_location):
    closest_distance = float('inf')
    closest_time = None
    closest_altitude = None

    for timestamp, aircraft_location, altitude in aircraft_data:
        _, _, distance = gps2dist_azimuth(aircraft_location[0], aircraft_location[1], seismometer_location[0], seismometer_location[1])
        distance = distance / 1000  # convert distance to kilometers
        if distance < closest_distance:
            closest_distance = distance
            closest_time = timestamp
            closest_altitude = altitude

    return closest_time, closest_distance, closest_altitude

def calculate_wave_arrival(closest_time, closest_distance, closest_altitude, aircraft_speed, aircraft_heading, seismometer_heading):
    speed_of_sound = 343  # speed of sound in m/s at sea level
    aircraft_speed_mps = aircraft_speed * 0.514  # convert aircraft speed from knots to m/s

    # calculate the component of the aircraft's speed in the direction of the seismometer
    relative_heading = radians(aircraft_heading - seismometer_heading)
    relative_speed = aircraft_speed_mps * cos(relative_heading)

    # calculate the time it takes for the sound to travel from the aircraft to the seismometer
    time_for_sound_to_travel = sqrt((closest_distance * 1000)**2 + (closest_altitude * 0.3048)**2) / (speed_of_sound + relative_speed)

    wave_arrival_time = closest_time + datetime.timedelta(seconds=time_for_sound_to_travel)

    return wave_arrival_time


seismo_data = pd.read_csv('input/all_sta.txt', sep="|")
seismo_latitudes = seismo_data['Latitude']
seismo_longitudes = seismo_data['Longitude']
station = seismo_data['Station']


flight_num = [530342801,528485724,528473220,528407493,528293430]
time = [1551066051,1550172833,1550168070,1550165577,1550089044]
sta = [1022,1272,1173,1283,1004]
day = [25,14,14,14,13]
for n in range(0,5):
	ht = datetime.datetime.utcfromtimestamp(time[n])
	mins = ht.minute
	secs = ht.second
	h = ht.hour
	tim = 120	
	h_u = str(h+1)
	if h < 23:			
		day2 = str(day[n])
		if h < 10:
			h_u = '0'+str(h+1)
			h = '0'+str(h)
		else:
			h_u = str(h+1)
			h = str(h)
	else:
		h_u = '00'
		day2 = str(day[n]+1)
	flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/201902'+str(day[n])+'_positions/201902'+str(day[n])+'_'+str(flight_num[n])+'.csv', sep=",")

