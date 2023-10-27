from datetime import datetime
import sys
import fileinput

input_file = 'dates1.txt'
output_file = 'dates2.txt'

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        dt = datetime.strptime(line.strip(), '%Y-%m-%dT%H:%M:%S.%f')
        rounded_dt = dt.replace(second=round(dt.second))
        f_out.write(rounded_dt.strftime('%Y-%m-%d %H:%M:%S') + '\n') 

for i, line in enumerate(fileinput.input('dates1.txt', inplace=1)):
        sys.stdout.write(line.replace(line[-, '')) 
