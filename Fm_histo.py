import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Initialize a list to store the data from the files
data = []

# Define the directory where your files are located
file = 'output/C185data_atmosphere_1o_updated.csv'

# Read the CSV file using pandas
df = pd.read_csv(file)

# Extract the 13th column
column = df.iloc[:, 12]
print(column)
# Append the data to the list
data.extend(column)


plt.figure()
# Plot the histogram
plt.hist(data) 
plt.axvline(np.median(data), color='k', linestyle='dashed', linewidth=1)
plt.text(np.median(data), plt.ylim()[0], str(np.round(np.median(data),3)), ha='left', va='bottom')
plt.title('Misfit of Atmosphere Corrected 1 Overtone Inversion')
plt.show()
