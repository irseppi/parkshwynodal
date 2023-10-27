import os
os.getcwd()
collection = '533595343'
for i, filename in enumerate(os.listdir(collection)):
	for p, fil in enumerate(os.listdir('533595343/'+filename)):
		
		os.rename('533595343/' + filename + '/' + fil, '533595343/' + filename +'_'+ fil)

