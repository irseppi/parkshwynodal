import sys
import fileinput

step = 3
 
with open("nodal_part.site") as handle:
	for lineno, line in enumerate(handle):
		if lineno % step == 0:
			print(line)
