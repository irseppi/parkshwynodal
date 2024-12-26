import json
import pandas as pd
import numpy as np

alt = 1.1 # replace with altitude provided by flightradar24 (convert from feet to km)
# Read JSON data from a file
with open('1549843272.8099847_64.84832950193729_-149.82614954364118.dat', 'r') as file:
    data = json.load(file)

# Extract metadata
metadata = data['metadata']
sourcefile = metadata['sourcefile']
datetime = metadata['time']['datetime']
latitude = metadata['location']['latitude']
longitude = metadata['location']['longitude']
parameters = metadata['parameters']

# Extract data
data_list = data['data']

# Convert data to a DataFrame
df = pd.DataFrame(data_list)

# Print extracted information
print(f"Source file: {sourcefile}")
print(f"Datetime: {datetime}")
print(f"Location: ({latitude}, {longitude})")
#print(f"Parameters: {parameters}")
#print(df)


# Find the "Z" parameter and extract the value at index 600
z_index = None
hold = np.inf
for item in data_list:
    if item['parameter'] == 'Z':
        for i in range(len(item['values'])):
            if abs(float(item['values'][i]) - float(alt)) < hold:
                hold = abs(float(item['values'][i]) - float(alt))
                z_index = i

for item in data_list:
    if item['parameter'] != 'Z0':
        print(str(item['description'])+'('+str(item['units'])+'): '+str(item['values'][z_index]))

    if item['parameter'] == 'Z0':
        print(str(item['description'])+'('+str(item['units'])+'): '+str(item['values'][0]))