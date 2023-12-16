import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.core import UTCDateTime
import datetime
import pytz
import obspy
import math 
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

def make_base_dir(base_dir):
	base_dir = Path(base_dir)
	if not base_dir.exists():
		current_path = Path("/")
		for parent in base_dir.parts:
			current_path = current_path/parent
			if not current_path.exists():
				current_path.mkdir()
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

text = open('input/all_station_crossing_db.txt', 'r')

airplane = 'input/20231010_Aircraft _UA_Fairbanks.csv'

plane_data = pd.read_csv(airplane, sep=",")
man = plane_data['MANUFACTURER']
model = plane_data['Model']
des = plane_data['Type Designator']
descrip = plane_data['Description']
engine = plane_data['Engine Type']
coun = plane_data['Engine Count']
turb_cat = plane_data['Wake Turbulence Category']
							
for i in range(len(day)):
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/plane_spec2/2019-'+month[i]+'-'+day[i]+'/'
	flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+'_flights.csv', sep=",")
	flight_id = flight_data['flight_id']
	equipment = flight_data['equip']
	callsign = flight_data['callsign'] 
	fly = flight_data['flight']

	for flight in os.listdir(spec_dir):
		f = os.path.join(spec_dir, flight)
		for l in range(len(flight_id)):
			if str(flight_id[l]) == str(flight):
				pla = equipment[l]
				for h in range(len(des)):
					if pla == des[h]:
						text2 = str(man[h]) + ', '+ str(model[h]) + ' (' + str(descrip[h]) + ')'
						break
					else:
						text2 = 'Callsign: ' + str(callsign[l])
					print(text2) 
				
			else:
				continue
		for station in os.listdir(f):
			sta = os.path.join(f, station)
			for image in os.listdir(sta):
				time = image[0:10]
				im = os.path.join(sta, image)
				
				#try:
				# Open images
				spectrogram = Image.open(im)
				map_img = Image.open('/scratch/irseppi/nodal_data/plane_info/pmap2/2019'+month[i]+day[i]+'/'+flight+'/map_'+station+'_'+time+'.png')
				zoom_map = Image.open('/scratch/irseppi/nodal_data/plane_info/pmap_zoom2/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/zmap_'+flight+'_' + str(time) + '.png')
				spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/spec2/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'+station+'_' + str(time) + '.png')

				# Resize images
				google_slide_width = 1280  # Width of a Google Slide in pixels
				google_slide_height = 720  # Height of a Google Slide in pixels
				
				path = '/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(pla)+'.jpg'
				if os.path.isfile(path):
					plane_img = Image.open(path)
					
				else:
					plane_img = Image.open('hold.png')
					#search text files for plane
				scale = 70/1280
				plane = plane_img.resize((int(google_slide_width * 0.27), int(google_slide_height * 0.27)))
				spec = spec_img.resize((int(google_slide_width * 0.31), int(google_slide_height * 0.36)))  
				zoom = zoom_map.resize((int(google_slide_width * 0.236), int(google_slide_height * 0.32)))
				maps = map_img.resize((int(google_slide_width *  0.256), int(google_slide_height * 0.32)))
				spectrogram = spectrogram.resize((int(google_slide_width * 0.75), int(google_slide_height)))
				#cover = Image.open('hold.png')
				#cov = cover.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.25)))

				# Create blank canvas
				canvas = Image.new('RGB', (google_slide_width, google_slide_height), 'white')

				# Paste images onto canvas
				
				#canvas.paste(plane, (google_slide_width - plane.width, 0))
				#canvas.paste(maps, (google_slide_width  - 520, plane.height))
				#canvas.paste(spec, (google_slide_width - spec.width + int(spec.width/12), google_slide_height - spec.height))
				#canvas.paste(zoom, (google_slide_width  - zoom.width + int(zoom.width/5.5), plane.height))
				#canvas.paste(spectrogram, (-40, 0))
				#width, height = zoom.size
 				
				canvas.paste(maps, (google_slide_width - int(maps.width*1.36), plane.height))
				
				canvas.paste(spec, (google_slide_width - spec.width+ int(spec.width/12), google_slide_height - spec.height))
				canvas.paste(zoom, (google_slide_width - zoom.width + int(zoom.width/4), plane.height))
				canvas.paste(plane, (google_slide_width - plane.width, 0))
				canvas.paste(spectrogram, (-40, 0))
				#canvas.paste(maps, (google_slide_width - int(maps.width*1.5), plane.height))
				
				# Draw text from files
				draw = ImageDraw.Draw(canvas)
				font = ImageFont.truetype('input/Arialn.ttf', 18) #load_default()  
				flight_file = '/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+ '_positions/2019'+month[i]+day[i]+ '_' + flight + '.csv'
				flight_data = pd.read_csv(flight_file, sep=",")
				t = flight_data['snapshot_id']
				speed = flight_data['speed']
				alt = flight_data['altitude']
				
				for l in range(len(t)):
					if int(time) == int(t[l]):
						text1 = 'Speed: '+str(round(speed[l]*0.514444,2))+'m/s\nAltitude: '+str(round(alt[l]*0.3048,2))+'m' 
					else:
						continue
				bbox = draw.textbbox((google_slide_width - plane.width, 0), text2, font=font)
				draw.rectangle(bbox, fill="white")
				draw.text((google_slide_width - plane.width, 0), text2, fill='black', font=font)
				
				draw.text((google_slide_width - 220, 430), text1, fill='black', font=font)

				BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/full_image2/'
				make_base_dir(BASE_DIR)
				name= BASE_DIR + '2019'+str(month[i])+str(day[i])+'_'+(flight)+'_'+time+'_'+str(station)+'_'+str(pla)+'_'+str(descrip[h])+'_'+str(engine[h])+str(coun[h])+'.png'

				# Save combined image
				canvas.save(name)
				#except:
					
				#continue

'''
try:
	map_img = 
except:
	map_image = Image.open('hold.png')

try:
	zoom_map = 
	
except:
	zoom_map =
try:
	spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/spec/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'+station+'_' + str(time) + '.png')
except:
	spec_img = Image.open('hold.png')

spectrogram = spectrogram.resize((int(google_slide_width * 0.75), int(google_slide_height)))
'''
						
						

											
