import sys
import fileinput

# replace all occurrences of ''' with '"'
for i, line in enumerate(fileinput.input('master_stations(_all).site', inplace=1)):
	sys.stdout.write(line.replace ('-    -     0.0000       0.0000  1687896210.00000', '-    -         0.0000    0.0000  1687896210.00000')) 

'''
for i, line in enumerate(fileinput.input('master_stations(_all).site', inplace=1)):
	val = line.split()
	new_val = round(float(val[4]), 4) # 3 digit precision
	sys.stdout.write(line.replace(str(val[4]), str(new_val)))



for i, line in enumerate(fileinput.input('master_stations(_all).site', inplace=1)):
	val = line.split()
	new_val = round(float(val[5])/1000,4)
	sys.stdout.write(line.replace(str(val[5]), str(new_val)))
'''
