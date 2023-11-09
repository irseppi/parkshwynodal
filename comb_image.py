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
	data = 'input/flight_name.txt' #'/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+'_flights.csv'
	flight_data = pd.read_csv(data, sep=",")
	flight_id = flight_data['flight_id']
	equipment = flight_data['equip']

	for flight in os.listdir(spec_dir):
		f = os.path.join(spec_dir, flight)
		for l in range(len(flight_id)):
			if str(flight_id[l]) == str(flight):
				p = equipment[l]

		for station in os.listdir(f):
			sta = os.path.join(f, station)
			for image in os.listdir(sta):
				im = os.path.join(sta, image)

				# Open images
				spectrogram = Image.open(im)
				map_img = Image.open('/scratch/irseppi/nodal_data/plane_info/plane_map/2019-'+month[i]+'-'+day[i]+'/map_2019'+month[i]+day[i]+'_'+flight+'.png')
				blown_map = Image.open('/scratch/irseppi/nodal_data/plane_info/map_zoom/2019'+month[i]+day[i]+'/'+flight+'/'+flight+'_'+station+'_' + str(time[l]) + '.png')
				#spec_amp = 
				#if p == 'nan':
				#plane_img = Image.open('plane.png')
				#else:

				plane_img = Image.open('/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(p)+'.jpg')
				# Resize images
				google_slide_width = 1280  # Width of a Google Slide in pixels
				google_slide_height = 720  # Height of a Google Slide in pixels

				spectrogram = spectrogram.resize((int(google_slide_width * 0.75), int(google_slide_height )))
				plane = plane_img.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.25)))
				#map_img = map_img.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.375)))  change to spec_amp and make map small in blown up map
				blown_map = blown_map.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.375)))
				
				# Create blank canvas
				canvas = Image.new('RGB', (google_slide_width, google_slide_height), 'white')

				# Paste images onto canvas
				canvas.paste(spectrogram, (0, 0))
				canvas.paste(plane, (google_slide_width - plane.width, 0))
				#canvas.paste(map_img, (google_slide_width - map_img.width, plane.height))
				canvas.paste(blown_map, (google_slide_width - blown_map.width, plane.height + map_img.height))

				# Draw text from files
				draw = ImageDraw.Draw(canvas)
				font = ImageFont.load_default()  # Or replace with the font you want to use

				#with open('file1.txt', 'r') as f:
				#text1 = f.read()
				#if p != '':
				draw.text((10, google_slide_height - 50), str(p), fill='black', font=font)

				#with open('file2.txt', 'r') as f:
				#text2 = f.read()
				#draw.text((10, google_slide_height - 30), text2, fill='black', font=font)

				BASE_DIR = "/scratch/irseppi/nodal_data/Plane_info/Full_image/2019-"+str(month[i])+"-"+str(day[i])+'/'
				make_base_dir(BASE_DIR)
				name= BASE_DIR + str(p)+'_'+station+'.png'
				# Save combined image
				canvas.save(name)

