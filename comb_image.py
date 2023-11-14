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
								
for i in range(len(day)):
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/plane_spec/2019-'+month[i]+'-'+day[i]+'/'
	data = 'input/flight_name.txt' 

	flight_data = pd.read_csv(data, sep=",")
	flight_id = flight_data['flight_id']
	equipment = flight_data['equip']
	date = flight_data['date']
	aircraft_id = flight_data['aircraft_id']
	reg = flight_data['reg']
	callsign = flight_data['callsign'] 
	fly = flight_data['flight']
	schd_from = flight_data['schd_from']
	schd_to = flight_data['schd_to']
	real_to = flight_data['real_to']
	reserved = flight_data['reserved']

	for flight in os.listdir(spec_dir):
		f = os.path.join(spec_dir, flight)
		for l in range(len(flight_id)):
			if str(flight_id[l]) == str(flight):
				p = equipment[l]
				text2 = 'Aircraft ID: ' + str(aircraft_id[l])+ '\nEquipment: ' + str(equipment[l]) + '\nCallsign: ' + str(callsign[l]) + '\nFrom: ' + str(schd_from[l])+'\nTo: '+str(schd_to[l])+'/'+str(real_to[l])
		for station in os.listdir(f):
			sta = os.path.join(f, station)
			for image in os.listdir(sta):
				time = image[0:10]
				im = os.path.join(sta, image)

				# Open images
				spectrogram = Image.open(im)
				map_img = Image.open('/scratch/irseppi/nodal_data/plane_info/plane_map/2019-'+month[i]+'-'+day[i]+'/map_2019'+month[i]+day[i]+'_'+flight+'.png')
				zoom_map = Image.open('/scratch/irseppi/nodal_data/plane_info/map_zoom/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/zmap_'+flight+'_' + str(time) + '.png')
				spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/spec/2019-'+month[i]+'-'+day[i]+'/'+flight+'/'+station+'/'+flight+'_' + str(time) + '.png')
				
				# Resize images
				google_slide_width = 1280  # Width of a Google Slide in pixels
				google_slide_height = 720  # Height of a Google Slide in pixels
				
				path = '/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(p)+'.jpg'
				if os.path.isfile(path):
					plane_img = Image.open('/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(p)+'.jpg')
					plane = plane_img.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.25)))
				else:
					plane = Image.open('hold.png')
					#search text files for plane
				
				
				spec = spec_img.resize((int(google_slide_width * 0.28), int(google_slide_height * 0.28)))  
				zoom = zoom_map.resize((int(google_slide_width * 0.258), int(google_slide_height * 0.31)))
				maps = map_img.resize((int(google_slide_width * 0.1346), int(google_slide_height * 0.31)))
				spectrogram = spectrogram.resize((int(google_slide_width * 0.75), int(google_slide_height)))

				# Create blank canvas
				canvas = Image.new('RGB', (google_slide_width, google_slide_height), 'white')

				# Paste images onto canvas
				
				canvas.paste(plane, (google_slide_width - plane.width, 0))
				canvas.paste(spec, (google_slide_width - spec.width+ int(spec.width/12), plane.height))
				canvas.paste(zoom, (google_slide_width - zoom.width + int(zoom.width/5.5), google_slide_height - zoom.height))
				canvas.paste(spectrogram, (-40, 0))
				canvas.paste(maps, (google_slide_width - int(maps.width*2.1), google_slide_height - maps.height))
				

				# Draw text from files
				draw = ImageDraw.Draw(canvas)
				font = ImageFont.truetype('input/TimesSansDisplay-PR9E.ttf', 16) #load_default()  
				#PIL.ImageFont.ImageFont.getsize
				flight_file = '/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+ '_positions/2019'+month[i]+day[i]+ '_' + flight + '.csv'
				flight_data = pd.read_csv(flight_file, sep=",")
				t = flight_data['snapshot_id']
				speed = flight_data['speed']
				alt = flight_data['altitude']
				
				for l in range(len(t)):
					if int(time) == int(t[l]):
						text1 = 'Date: 2019-' + month[i] + '-' + day[i] + '\nFlight: ' + flight + '\nStation: ' + station + '\nSpeed: '+str(round(speed[l]*0.514444,2))+'m/s\nAltitude: '+str(round(alt[l]*0.3048,2))+'m' 
					else:
						continue
				draw.text((google_slide_width - 140, 400), text1,fill='black', font=font,)
				
				draw.text((google_slide_width - 300, 400), text2, fill='black', font=font)

				BASE_DIR = "/scratch/irseppi/nodal_data/plane_info/full_image/2019-"+str(month[i])+"-"+str(day[i])+'/'+flight+'/'
				make_base_dir(BASE_DIR)
				name= BASE_DIR + time+'_'+station+'.png'

				# Save combined image
				canvas.save(name)

