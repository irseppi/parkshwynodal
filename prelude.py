import sys
import fileinput
from datetime import datetime
import os
import pandas as pd

def modify_file(input_file_name, output_file_name):
    # Function to modify the content of the file
    def modify_content(content):
        return content.upper()

    # Read the input file
    with open(input_file_name, 'r') as file:
        content = file.read()

    # Modify the content
    modified_content = modify_content(content)

    # Write the modified content to the output file
    with open(output_file_name, 'w') as file:
        file.write(modified_content)

modify_file('input.txt', 'output.txt')
def station_subset(filename,steps,outputname):
	# Read the input file
	with open(filename) as handle:
		for lineno, line in enumerate(handle):
		if lineno % step == 0:
			print(line)
	return 


# replace all occurrences of ''' with '"'
for i, line in enumerate(fileinput.input('master_stations(_all).site', inplace=1)):
	sys.stdout.write(line.replace ('-    -     0.0000       0.0000  1687896210.00000', '-    -         0.0000    0.0000  1687896210.00000')) 


for i, line in enumerate(fileinput.input('master_stations(_all).site', inplace=1)):
	val = line.split()
	new_val = round(float(val[4]), 4) # 3 digit precision
	sys.stdout.write(line.replace(str(val[4]), str(new_val)))



for i, line in enumerate(fileinput.input('master_stations(_all).site', inplace=1)):
	val = line.split()
	new_val = round(float(val[5])/1000,4)
	sys.stdout.write(line.replace(str(val[5]), str(new_val)))



# replace all occurrences of ''' with '"'
for i, line in enumerate(fileinput.input('filename.site', inplace=1)):
    sys.stdout.write(line.replace(' '', ' "')) 
 


os.getcwd()
collection = '533595343'
for i, filename in enumerate(os.listdir(collection)):
	for p, fil in enumerate(os.listdir('533595343/'+filename)):
		
		os.rename('533595343/' + filename + '/' + fil, '533595343/' + filename +'_'+ fil)

input_file = 'color_codes.txt'
output_file = 'colors.txt'

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
       line=line.split(' 	')
       print(line[0])
       f_out.write(line[0]) 

for i, line in enumerate(fileinput.input('colors.txt', inplace=1)):
    sys.stdout.write(line.replace('#', '\n#')) 



input_file = 'dates1.txt'
output_file = 'dates2.txt'

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        dt = datetime.strptime(line.strip(), '%Y-%m-%dT%H:%M:%S.%f')
        rounded_dt = dt.replace(second=round(dt.second))
        f_out.write(rounded_dt.strftime('%Y-%m-%d %H:%M:%S') + '\n') 

for i, line in enumerate(fileinput.input('dates1.txt', inplace=1)):
        sys.stdout.write(line.replace(line[-, '')) 

f = open('coun.txt', 'w')
text = open('all_station_crossing_db.txt', 'r')
flight_data = pd.read_csv('20231010_Aircraft _UA_Fairbanks.csv', sep=",")
eq = flight_data['TypeDesignator']
des = flight_data['Description']
count = 0
for line in text.readlines():
	val = line.split(',')
	equip = val[6]
	for l in range(len(eq)):
		if str(eq[l]) == str(equip[0:4]) and str(des[l]) == 'Helicopter':
			count = count + 1
			f.write(eq[l]+'\n')
f.write(count)
f.close()
text = open('all_station_crossing_db.txt', 'r')

for line in text.readlines():
	val = line.split(',')
	equip = val[6]
	print(equip)
