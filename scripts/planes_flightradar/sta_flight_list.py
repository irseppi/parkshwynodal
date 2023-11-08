import os

month = []
day = []

for d in range(11,29):
	day.append(str(d))
	month.append('02')

for d in range(1, 10):
	day.append('0' + str(d))
	month.append('03')

for d in range(10, 27):
	day.append(str(d))
	month.append('03')

files = open('output/flight_name_sta.txt','w')

for i in range(len(day)):
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/plane_spec/2019-'+month[i]+'-'+day[i]+'/'

	for flight in os.listdir(spec_dir):
		f = os.path.join(spec_dir, flight)
		for station in os.listdir(f):
			sta = os.path.join(f, station)
			files.write('2019'+month[i]+day[i]+','+flight+','+station+'\n')

files.close()				

