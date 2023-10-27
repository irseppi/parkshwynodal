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
for month in range(2,4):
	if month == 2:
		for day in range(11,29):
			day.append(str(day))
			month.append('02')
			
	elif month == 3:
		for day in range(1, 27):
			if day < 10:
				day.append('0' + str(day))
				month.append('03')
			else:
				day.append(str(day))
				month.append('03')
plane_pic = []
plane_map = []
plane_spec = []
b_map = []
data = []
for month in range(2,4):
	if month == 2:
		month = '02'
		for day in range(11,29):
			day = str(day)
			# assign directory
			directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_positions'
			
			# iterate over files in directory
			for filename in os.listdir(directory):
				
				f = os.path.join(directory, filename)
				
				# checking if it is a file
				if os.path.isfile(f):
					data.append(filename)
			DIR = '/scratch/irseppi/nodal_data/Plane_info/Plane_spec/2019-'+month+'-'+day
			# iterate over files in directory
			for filename in os.listdir(DIR):
				
				f = os.path.join(directory, filename)
				for fil in os.listdir(f):
					# checking if it is a file
					if os.path.isfile(f):
						plane_spec.append(filename)git clone https://GITHUBUSERNAME@github.com/irseppi/denali_nodal_set.git

			D = '/scratch/irseppi/nodal_data/Plane_info/Plane_map/2019-'+month+'-'+day
			# iterate over files in directory
			for filename in os.listdir(D):
				f = os.path.join(directory, filename)
	
				# checking if it is a file
				if os.path.isfile(f):
					plane_map.append(filename)
	elif month == 3:
		month = '03'
		for day in range(1, 27):
			if day < 10:
				day = '0' + str(day)
				# assign directory
				directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_positions'
				DIR = '/scratch/irseppi/nodal_data/Plane_info/'

				# iterate over files in directory
				for filename in os.listdir(directory):
					
					f = os.path.join(directory, filename)
					
					# checking if it is a file
					if os.path.isfile(f):
						filenames.append(filename)
			else:
				day = str(day)
				# assign directory
				directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month+day+'_positions'
				DIR = '/scratch/irseppi/nodal_data/Plane_info/'
				# iterate over files in directory
				for filename in os.listdir(directory):
					filenames.append(filename)
					f = os.path.join(directory, filename)
					
					# checking if it is a file
					if os.path.isfile(f):
						filenames.append(filename)
								
for i in range(len(day)):
	DIR = '/scratch/irseppi/nodal_data/Plane_info/'
	directory = '/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+'_positions'
	for m, f in enumerate(DIR):
		# Open images
		spectrogram = Image.open(DIR + '/Plane_spec/2019-'+month[i]+'-'+day[i]+'/'+spectrogram.jpg')
		plane = Image.open('plane.jpg')
		map_img = Image.open('map.jpg')
		blown_map = Image.open('blown_map.jpg')

		# Resize images
		google_slide_width = 1280  # Width of a Google Slide in pixels
		google_slide_height = 720  # Height of a Google Slide in pixels

		spectrogram = spectrogram.resize((google_slide_width, int(google_slide_height * 0.75)))
		plane = plane.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.25)))
		map_img = map_img.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.25)))
		blown_map = blown_map.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.25)))

		# Create blank canvas
		canvas = Image.new('RGB', (google_slide_width, google_slide_height), 'white')

		# Paste images onto canvas
		canvas.paste(spectrogram, (0, 0))
		canvas.paste(plane, (google_slide_width - plane.width, 0))
		canvas.paste(map_img, (google_slide_width - map_img.width, plane.height))
		canvas.paste(blown_map, (google_slide_width - blown_map.width, plane.height + map_img.height))

		# Draw text from files
		draw = ImageDraw.Draw(canvas)
		font = ImageFont.load_default()  # Or replace with the font you want to use

		with open('file1.txt', 'r') as f:
		    text1 = f.read()
		draw.text((10, google_slide_height - 50), text1, fill='black', font=font)

		with open('file2.txt', 'r') as f:
		    text2 = f.read()
		draw.text((10, google_slide_height - 30), text2, fill='black', font=font)

		BASE_DIR = "/scratch/irseppi/nodal_data/Plane_info/Full_image/2019-0"+str(month)+"-"+str(day)
		make_base_dir(BASE_DIR)

		# Save combined image
		canvas.save('/scratch/irseppi/nodal_data/Plane_info/Full_image/2019-0'+str(month)+"-"+str(day)+'.jpg')
