import pandas as pd
import os
from PIL import Image, ImageDraw, ImageFont
from prelude import make_base_dir

day = ['11','13','14','18','21','24','25','04','04','22','22','23']
month =  ['02','02','02','02','02','02','02','03','03','02','02','02']

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
	spec_dir = '/scratch/irseppi/nodal_data/plane_info/5plane_spec/2019-'+month[i]+'-'+day[i]
	
	flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+'_flights.csv', sep=",")
	flight_id = flight_data['flight_id']
	equipment = flight_data['equip']
	callsign = flight_data['callsign'] 
	fly = flight_data['flight']

	for flight in os.listdir(spec_dir):
		f = os.path.join(spec_dir, flight)
		for l in range(len(flight_id)):
			if str(flight_id[l]) == str(flight):

				#search text files for plane
				pla = equipment[l]
				for h in range(len(des)):
					if pla == des[h]:
						text2 = str(des[h]) + ': ' + str(man[h]) + ', '+ str(model[h]) + ' (' + str(descrip[h]) + ')'
						text3 = 'Flight Designator: '+str(des[h]) + '\nEngine: '+str(engine[h])+'\nEngine Count: ' +str(coun[h])+'\nWake Turbulence Category: '+str(coun[h])
						break
					else:
						text2 = 'Callsign: ' + str(callsign[l])
						text3 = 'Flight Designator: Unknown\nEngine: Unknown\nEngine Count: Unknown\nWake Turbulence Category: Unknown'
					print(text2) 
				
			else:
				continue
		for station in os.listdir(f):
			sta = os.path.join(f, station)
			for image in os.listdir(sta):
				time = image[0:10]
				im = os.path.join(sta, image)
				
				
				# Open images
				spectrogram = Image.open(im)
				map_img = Image.open('/scratch/irseppi/nodal_data/plane_info/map_all/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/map_'+flight+'_'+time+'.png')
				spec_img = Image.open('/scratch/irseppi/nodal_data/plane_info/5spec/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'+station+'_' + str(time) + '.png')

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
				font = ImageFont.truetype('input/Arial.ttf', 15) #load_default()  
				flight_file = '/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+ '_positions/2019'+month[i]+day[i]+ '_' + flight + '.csv'
				flight_data = pd.read_csv(flight_file, sep=",")
				t = flight_data['snapshot_id']
				speed = flight_data['speed']
				alt = flight_data['altitude']
				font2 = ImageFont.truetype('input/Arial.ttf', 25)

				# Label each image
				draw.text((15, 35), '(a)', fill='black', font=font2)
				draw.text((15, 350), '(b)', fill='black', font=font2)
				draw.text((google_slide_width - int(plane.width*1.15), 20), '(c)', fill='black', font=font2)
				draw.text((google_slide_width - int(plane.width*1.15), int(plane.height) + int(plane.height*0.05)), '(d)', fill='black', font=font2)
				draw.text((google_slide_width - spec.width + int(spec.width/12) - 15, google_slide_height - spec.height + 20), '(e)', fill='black', font=font2)

				for l in range(len(t)):
					if int(time) == int(t[l]):
						text1 = 'Speed: '+str(round(speed[l]*0.514444,2))+' m/s\n           : '+str(round(speed[l]*1.15078,2))+' mph\nAltitude: '+str(round(alt[l]*0.3048,2))+' m\n            : '+str(round(alt[l],2)) +' ft'
					else:
						continue
				
				bbox = draw.textbbox((google_slide_width - plane.width, 0), text2, font=font)
				draw.rectangle(bbox, fill="white")
				draw.text((google_slide_width - plane.width, 0), text2, fill='black', font=font)			
				draw.text((google_slide_width - 150, 405), text1,fill='black', font=font)
				draw.text((google_slide_width - 365, 405), text3, fill='black', font=font)
				
				BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/5fig/'
				make_base_dir(BASE_DIR)
				name= BASE_DIR + '2019'+str(month[i])+str(day[i])+'_'+(flight)+'_'+time+'_'+str(station)+'_'+str(pla)+'_'+str(descrip[h])+'_'+str(engine[h])+str(coun[h])+'.png'

				# Save combined image
				canvas.save(name)

					
					

										
