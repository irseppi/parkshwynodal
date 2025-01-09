import pandas as pd
import os
from PIL import Image, ImageDraw, ImageFont
from prelude import make_base_dir
import glob
import numpy as np
d = []
m = []

month = '02'
for day in range(11, 29):
	day = str(day)
	d.append(day)
	m.append(month)

month = '03'
for day in range(1, 26):
	if day < 10:
		day = '0' + str(day)
		d.append(day)
		m.append(month)
	else:
		day = str(day)
		d.append(day)
		m.append(month)
day = d
month = m

text = open('input/all_station_crossing_db.txt', 'r')

airplane = 'input/20231010_Aircraft_UA_Fairbanks.csv'

plane_data = pd.read_csv(airplane, sep=",")
man = plane_data['MANUFACTURER']
model = plane_data['Model']
des = plane_data['Type Designator']
descrip = plane_data['Description']
engine = plane_data['Engine Type']
coun = plane_data['Engine Count']
turb_cat = plane_data['Wake Turbulence Category']
							
for i in range(len(day)):
	#try:
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/C185_spec_c/2019-'+month[i]+'-'+day[i]
	if os.path.exists(spec_dir):
		flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+'_flights.csv', sep=",")
		flight_id = flight_data['flight_id']
		equipment = flight_data['equip']
		callsign = flight_data['callsign'] 
		fly = flight_data['flight']
		aircraft_id = flight_data['aircraft_id']
		for flight in os.listdir(spec_dir):
			f = os.path.join(spec_dir, flight)
			
			for station in os.listdir(f):
				sta = os.path.join(f, station)
				for image in os.listdir(sta):
					split_array = np.array(image.split('_'))
					time = str(split_array[0])
					im = os.path.join(sta, image)
					for l in range(len(flight_id)):
						if str(flight_id[l]) == str(flight):
							flight_file = '/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+ '_positions/2019'+month[i]+day[i]+ '_' + flight + '.csv'
							flight_data = pd.read_csv(flight_file, sep=",")
							times = flight_data['snapshot_id']
							speed = flight_data['speed']
							alt = flight_data['altitude']
							head = flight_data['heading']
							dist = 0
							deg = 0
							temp = 0
							wind = 0 
							sound = 0
							eff_sound = 0
							az = 0
							qnum = 0
							mnum = 0
							font2 = ImageFont.truetype('input/Arial.ttf', 25)
							diff = np.inf
							for t in range(len(times)):
								if abs(float(time) - float(times[t])) < diff:
									diff = abs(float(time) - float(times[t]))
									text1 = 'Altitude: '+str(round(alt[t]*0.3048,2))+' m ('+str(round(alt[t],2)) +' ft)\nDistance: '+str(round(dist,2))+' m\nVelocity: '+str(round(speed[t]*0.514444,2))+' m/s ('+str(round(speed[2]*1.15078,2))+' mph)\n               at '+str(round(deg,2))+ '\N{DEGREE SIGN}' + '\nHeading: '+str(round(head[2],2))+ '\N{DEGREE SIGN}'
									text2 = 'Temperature: '+str(round(temp,2))+'\N{DEGREE SIGN}'+'C\nWind: '+str(round(wind,2))+' m/s\nSound Speed: '+str(round(sound,2))+' m/s\nEffective Sound Speed:\n '+str(round(eff_sound,2))+' m/s at '+str(round(az,2))+ '\N{DEGREE SIGN}'
								else:
									continue
							#search text files for plane
							pla = equipment[l]
							id = aircraft_id[l]
							for h in range(len(des)):
								if pla == des[h]:
									text3 = 'Callsign: ' +  str(callsign[l]) + ' (' + str(des[h]) + ')'

									break
								else:
									text3 = 'Callsign: ' + str(callsign[l])
						else:
							continue
					# Open images
					spectrogram = Image.open(im)
					# Get the path of the image file using a wildcard
					image_path = glob.glob('/scratch/irseppi/nodal_data/plane_info/map_all/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/map_'+flight+'_*')[0]

					map_img = Image.open(image_path)
					spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/C185_specrum_c/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'+station+'_' + str(time) + '.png')

					# Resize images
					google_slide_width = 1280  # Width of a Google Slide in pixels
					google_slide_height = 720  # Height of a Google Slide in pixels
					
					path = '/scratch/irseppi/nodal_data/plane_info/plane_images/'+str(pla)+'.jpg'
					if os.path.isfile(path):
						plane_img = Image.open(path)
						
					else:
						plane_img = Image.open('hold.png')
						
					scale = 70/1280
					plane = plane_img.resize((int(google_slide_width * 0.26), int(google_slide_height * 0.26)))
					spec = spec_img.resize((int(google_slide_width * 0.31), int(google_slide_height * 0.35)))  
					#maps = map_img.resize((int(google_slide_width *  0.31), int(google_slide_height * 0.27)))
					maps = map_img.resize((int(google_slide_width *  0.25), int(google_slide_width * 0.25 * map_img.height / map_img.width)))
					spectrogram = spectrogram.resize((int(google_slide_width * 0.75), int(google_slide_height)))

					# Create blank canvas
					canvas = Image.new('RGB', (google_slide_width, google_slide_height), 'white')

					# Paste images onto canvas
					canvas.paste(spec, (google_slide_width - spec.width+ int(spec.width/12), google_slide_height - spec.height))
					
					canvas.paste(maps, (google_slide_width - int(maps.width*1.05), int(plane.height)-int(plane.height*0.1)))
					canvas.paste(plane, (google_slide_width - plane.width, 0))
					canvas.paste(spectrogram, (-40, 0))
					# Draw text from files
					draw = ImageDraw.Draw(canvas)
					font = ImageFont.truetype('input/Arial.ttf', 14) #load_default()  

					# Label each image
					draw.text((15, 35), '(a)', fill='black', font=font2)
					draw.text((google_slide_width - int(plane.width*1.5), 35), 'Q#: '+str(qnum), fill='black', font=font2)
					draw.text((google_slide_width - int(plane.width*1.5), google_slide_height - spec.height - spec.height/2), '[M'+str(mnum)+']', fill='black', font=font2)
					draw.text((15, 350), '(b)', fill='black', font=font2)
					draw.text((google_slide_width - int(plane.width*1.15), 20), '(c)', fill='black', font=font2)
					draw.text((google_slide_width - int(plane.width*1.15), int(plane.height) + int(plane.height*0.05)), '(d)', fill='black', font=font2)
					draw.text((google_slide_width - spec.width + int(spec.width/12) - 15, google_slide_height - spec.height + 20), '(e)', fill='black', font=font2)

					draw.text((google_slide_width - 370, 405), text1, fill='black', font=font)			
					draw.text((google_slide_width - 155, 405), text2,fill='black', font=font)
					bbox = draw.textbbox((google_slide_width - plane.width, 0), text3, font=font)
					draw.rectangle(bbox, fill="white")
					draw.text((google_slide_width - plane.width, 0), text3, fill='black', font=font)
					#show image

					BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/C185_atmosphere_correction/'
					make_base_dir(BASE_DIR)
					name= BASE_DIR + '2019'+str(month[i])+str(day[i])+'_'+(flight)+'_'+time+'_'+str(station)+'_'+str(pla)+'_'+str(id)+'_'+str(descrip[h])+'_'+str(engine[h])+str(coun[h])+'.png'

					# Save combined image
					canvas.save(name)

	else:
		continue
	#except:
	#	continue				
					

										
