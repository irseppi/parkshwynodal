import sys
import fileinput

# replace all occurrences of ''' with '"'
for i, line in enumerate(fileinput.input('filename.site', inplace=1)):
    sys.stdout.write(line.replace(' '', ' "')) 
 

