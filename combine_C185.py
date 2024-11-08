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
    spec_dir = '/scratch/irseppi/nodal_data/plane_info/C185_spec_1o/2019-'+month[i]+'-'+day[i]
    if os.path.exists(spec_dir):

        flight_data = pd.read_csv('/scratch/irseppi/nodal_data/flightradar24/2019'+month[i]+day[i]+'_flights.csv', sep=",")
        flight_id = flight_data['flight_id']
        equipment = flight_data['equip']
        callsign = flight_data['callsign'] 
        fly = flight_data['flight']
        aircraft_id = flight_data['aircraft_id']
        for flight in os.listdir(spec_dir):
            f = os.path.join(spec_dir, flight)
            for l in range(len(flight_id)):
                if str(flight_id[l]) == str(flight):

                    #search text files for plane
                    pla = equipment[l]
                    id = aircraft_id[l]
                    for h in range(len(des)):
                        if pla == des[h]:
                            text2 = str(des[h]) + ': ' + str(man[h]) + ', '+ str(model[h]) + ' (' + str(descrip[h]) + ')'
                            text3 = 'Aircraft ID: '+str(id) + '\nEngine: '+str(engine[h])+'\nEngine Count: ' +str(coun[h])+'\nWake Turbulence Category: '+str(coun[h])
                            break
                        else:
                            text2 = 'Callsign: ' + str(callsign[l])
                            text3 = 'Aircraft ID: Unknown\nEngine: Unknown\nEngine Count: Unknown\nWake Turbulence Category: Unknown'
                        #print(text2) 
                    
                else:
                    continue
            
            for station in os.listdir(f):
                
                sta = os.path.join(f, station)
                for image in os.listdir(sta):
                    time = image[0:10]
                    im = os.path.join(sta, image)
                    
                    # Open images
                    spectrogram = Image.open(im)
                    # Get the path of the image file using a wildcard
                    image_path = glob.glob('/scratch/irseppi/nodal_data/plane_info/map_all_UTM/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/map_'+flight+'_*')[0]

                    map_img = Image.open(image_path)
                    for image in os.listdir('/scratch/irseppi/nodal_data/plane_info/C185_spectrum_1o/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/'):
                         spectrum = os.path.join('/scratch/irseppi/nodal_data/plane_info/C185_spectrum_1o/2019'+month[i]+day[i]+'/'+flight+'/'+station+'/', image)
                    spec_img = Image.open(spectrum)

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

                    time_array = np.array(t) - float(time)
                    closest_time = min(time_array, key=abs)
                    closest_time_index = np.argmin(np.abs(time_array - closest_time))
                    for l in range(len(t)):
                        if int(t[l]) == int(t[closest_time_index]):
                            text1 = 'Speed: '+str(round(speed[l]*0.514444,2))+' m/s\n           : '+str(round(speed[l]*1.15078,2))+' mph\nAltitude: '+str(round(alt[l]*0.3048,2))+' m\n            : '+str(round(alt[l],2)) +' ft'
                        else:
                            continue
                    
                    bbox = draw.textbbox((google_slide_width - plane.width, 0), text2, font=font)
                    draw.rectangle(bbox, fill="white")
                    draw.text((google_slide_width - plane.width, 0), text2, fill='black', font=font)			
                    draw.text((google_slide_width - 150, 405), text1,fill='black', font=font)
                    draw.text((google_slide_width - 365, 405), text3, fill='black', font=font)
                    
                    BASE_DIR = '/scratch/irseppi/nodal_data/plane_info/C185_fig_1o/'

                    make_base_dir(BASE_DIR)
                    name= BASE_DIR + '2019'+str(month[i])+str(day[i])+'_'+(flight)+'_'+time+'_'+str(station)+'_'+str(pla)+'_'+str(id)+'_'+str(descrip[h])+'_'+str(engine[h])+str(coun[h])+'.png'

                    # Save combined image
                    canvas.save(name)

    else:
        continue
		
                

                                    
