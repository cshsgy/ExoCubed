from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

# read the data
data_path = "/home/linfel/ExoCubedlinfel/linfel/last2_pres_hotjupiter.nc"
dataset = Dataset(data_path,'r')
time = dataset.variables['time'][:]
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
press = dataset.variables['press'][:]
temp = dataset.variables['temp'][:]


x = lon
y = press
data = temp[-1,:,45,:]
# Create a figure and a set of subplots
plt.figure(figsize=(10, 6))

# Plot contour
contour = plt.contour(x, y, data)
plt.colorbar(contour, label='Temperature')

# Plot pcolor
#pc = plt.pcolor(x, y, data, shading='auto')
#plt.colorbar(pc, label='Temperature')

# Labels and title
plt.xlabel('Longitude')
plt.ylabel('Pressure (log scale)')
plt.title('Temperature Distribution with Logarithmic Pressure Axis')

# Set the y-axis to a logarithmic scale
plt.yscale('log')

# Inverting the y-axis if pressure increases with depth/altitude
plt.gca().invert_yaxis()

# Adjust the ticks on the y-axis to be more readable
#plt.yticks(pressure, labels=np.round(np.log10(pressure), 2))

plt.savefig("images/test3.png",dpi=300)
