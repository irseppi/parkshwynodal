import sys
import fileinput

input_file = 'color_codes.txt'
output_file = 'colors.txt'

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
       line=line.split(' 	')
       print(line[0])
       f_out.write(line[0]) 

for i, line in enumerate(fileinput.input('colors.txt', inplace=1)):
    sys.stdout.write(line.replace('#', '\n#')) 
