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

for d in range(27,29):
	day.append(str(d))
	month.append('02')

for d in range(1, 10):
	day.append('0' + str(d))
	month.append('03')

for d in range(10, 27):
	day.append(str(d))
	month.append('03')

text = open('all_station_crossing_db.txt', 'r')

data = 'input/flight_name.txt' 
airplane = 'input/20231010_Aircraft _UA_Fairbanks.csv'

flight_data = pd.read_csv(data, sep=",")
flight_id = flight_data['flight_id']
equipment = flight_data['equip']
callsign = flight_data['callsign'] 
fly = flight_data['flight']

plane_data = pd.read_csv(airplane, sep=",")
man = plane_data['MANUFACTURER']
model = plane_data['Model']
des = plane_data['Type Designator']
descrip = plane_data['Description']
engine = plane_data['Engine Type']
coun = plane_data['Engine Count']
turb_cat = plane_data['Wake Turbulence Category']
							
for i in range(len(day)):
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/perm_spec/2019-'+month[i]+'-'+day[i]+'/'

	for flight in os.listdir(spec_dir):
		f = os.path.join(spec_dir, flight)
		for l in range(len(flight_id)):
			if str(flight_id[l]) == str(flight):
				p = equipment[l]
				if p == float('nan'):
					text2 = 'Callsign: ' + str(callsign[l]) 
				else:
					text2 = str(man) + ', '+ str(model) + ' (' + str(descrip) + ')'
		for station in os.listdir(f):
			sta = os.path.join(f, station)
			for image in os.listdir(sta):
				time = image[0:10]
				im = os.path.join(sta, image)

				# Open images
				spectrogram = Image.open(im)
				map_img = Image.open('/scratch/irseppi/nodal_data/plane_info/plane_map/2019-'+month[i]+'-'+day[i]+'/map_2019'+month[i]+day[i]+'_'+flight+'.png')
				zoom_map = Image.open('/scratch/irseppi/nodal_data/plane_info/map_zoom/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/zmap_'+flight+'_' + str(time) + '.png')
				try:
					spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/spec/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'+station+'_' + str(time) + '.png')
				except:
					spec_img = Image.open('hold.png')
				# Resize images
				google_slide_width = 1280  # Width of a Google Slide in pixels
				google_slide_height = 720  # Height of a Google Slide in pixels
				
				path = '/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(p)+'.jpg'
				if os.path.isfile(path):
					plane_img = Image.open('/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(p)+'.jpg')
					
				else:
					plane = Image.open('hold.png')
					#search text files for plane
				
				plane = plane_img.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.25)))
				spec = spec_img.resize((int(google_slide_width * 0.31), int(google_slide_height * 0.31)))  
				zoom = zoom_map.resize((int(google_slide_width * 0.258), int(google_slide_height * 0.31)))
				maps = map_img.resize((int(google_slide_width * 0.1346), int(google_slide_height * 0.31)))
				spectrogram = spectrogram.resize((int(google_slide_width * 0.75), int(google_slide_height * 0.042)))
				cover = Image.open('hold.png')
				cov = cover.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.25)))
				# Create blank canvas
				canvas = Image.new('RGB', (google_slide_width, google_slide_height), 'white')

				# Paste images onto canvas
				
				canvas.paste(plane, (google_slide_width - plane.width, 0))
				canvas.paste(spec, (google_slide_width - spec.width+ int(spec.width/12), google_slide_height - spec.height))
				canvas.paste(zoom, (google_slide_width - zoom.width + int(zoom.width/5.5), plane.height))
				canvas.paste(spectrogram, (-40, 0))
				canvas.paste(maps, (google_slide_width - int(maps.width*2.1), plane.height))
				

				# Draw text from files
				draw = ImageDraw.Draw(canvas)
				font = ImageFont.truetype('input/Arialn.ttf', 18) #load_default()  
				#PIL.ImageFont.ImageFont.getsize
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
				draw.text((google_slide_width - plane.width, 0), text1, fill='black', font=font)
				
				draw.text((google_slide_width - 340, 400), text2, fill='black', font=font)

				BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/full_image/'
				make_base_dir(BASE_DIR)
				name= BASE_DIR + '2019'+str(month[i])+str(day[i])+'_'+flight+'_'+time+'_'+station+'.png'

				# Save combined image
				canvas.save(name)

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

#text = open('all_station_crossing_db.txt', 'r')
								
for i in range(len(day)):
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/perm_spec/2019'+month[i]+day[i]+'/'

	for flight in os.listdir(spec_dir):
		f = os.path.join(spec_dir, flight)
		for l in range(len(flight_id)):
			if str(flight_id[l]) == str(flight):
				pla = equipment[l]
				if pla == float('nan'):
					text2 = 'Callsign: ' + str(callsign[l]) 
				else:
					text2 = str(man) + ', '+ str(model) + ' (' + str(descrip) + ')'
				for station in os.listdir(f):
					sta = os.path.join(f, station)
					
					for image in os.listdir(sta):
						time = image[10:20]
						
						im = os.path.join(sta, image)
						
						# Open images
						spectrogram = Image.open(im)
						p = '/scratch/irseppi/nodal_data/plane_info/pmap/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'+'zmap_'+image
						
						try:
							map_img = Image.open(p)
						except:
							map_image = Image.open('hold.png')
						t ='/scratch/irseppi/nodal_data/plane_info/pmap_zoom/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'+'zmap_'+image
						
						try:
							zoom_map = Image.open(t)
							
							g = '/scratch/irseppi/nodal_data/plane_info/spec/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'+station+'_' + str(time) + '.png'
						except:
							zoom_map = Image.open('hold.png')
						try:
							spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/spec/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'+station+'_' + str(time) + '.png')
						except:
							spec_img = Image.open('hold.png')
						# Resize images
						google_slide_width = 1280  # Width of a Google Slide in pixels
						google_slide_height = 720  # Height of a Google Slide in pixels
						
						path = '/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(pla)+'.jpg'
						print(path)
						if os.path.isfile(path):
							plane_img = Image.open('/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(pla)+'.jpg')
							
						else:
							plane_img = Image.open('hold.png')
							#search text files for plane
						
						plane = plane_img.resize((int(google_slide_width * 0.25), int(google_slide_height * 0.25)))
						spec = spec_img.resize((int(google_slide_width * 0.31), int(google_slide_height * 0.31)))  
						zoom = zoom_map.resize((int(google_slide_width * 0.258), int(google_slide_height * 0.31)))
						maps = map_img.resize((int(google_slide_width * 0.1346), int(google_slide_height * 0.31)))
						spectrogram = spectrogram.resize((int(google_slide_width * 0.75), int(google_slide_height)))

						# Create blank canvas
						canvas = Image.new('RGB', (google_slide_width, google_slide_height), 'white')

						# Paste images onto canvas
						
						canvas.paste(plane, (google_slide_width - plane.width, 0))
						canvas.paste(spec, (google_slide_width - spec.width+ int(spec.width/12), google_slide_height - spec.height))
						canvas.paste(zoom, (google_slide_width - zoom.width + int(zoom.width/5.5), plane.height))
						canvas.paste(spectrogram, (-40, 0))
						canvas.paste(maps, (google_slide_width - int(maps.width*2.1), plane.height))
						

						# Draw text from files
						draw = ImageDraw.Draw(canvas)
						font = ImageFont.truetype('Arialn.ttf', 18) #load_default()  
						#PIL.ImageFont.ImageFont.getsize
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
						draw.text((google_slide_width - 160, 400), text1,fill='black', font=font,)
						
						draw.text((google_slide_width - 340, 400), text2, fill='black', font=font)

						BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/pfull_image/'
						make_base_dir(BASE_DIR)
						name= BASE_DIR + '2019'+str(month[i])+str(day[i])+'_'+flight+'_'+time+'_'+station+'.png'

						# Save combined image
						canvas.save(name)
						
