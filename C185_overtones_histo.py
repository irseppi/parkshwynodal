import os
import matplotlib.pyplot as plt
import numpy as np


# Define the directory where your files are located
directory = '/scratch/irseppi/nodal_data/plane_info/overtonepicks/' 


# Initialize a list to store the data from the files
data = []

# Iterate over all files in the directory and its subdirectories
for root, dirs, files in os.walk(directory):
    for file in files:
        # Check if the file has the desired format (e.g., .txt, .csv)
        if file.endswith('.csv'):
            # Construct the full file path
            file_path = os.path.join(root, file)
            f = open(file_path, 'r')
            # Read the data from the file and append it to the list
            for line in f.readlines():
                lines = line.split(',')
                peak = float(lines[1])
                data.append(peak)
plt.figure()
# Plot the histogram
plt.hist(np.round(data,1),125)

# Add tick marks at every value of 5
plt.xticks(np.arange(15, max(data)+1, 5))

plt.show()
