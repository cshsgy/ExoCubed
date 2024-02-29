from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np


"""User specify"""
# input nc file
data_path = "/home/linfel/data/hot_jupiter/hotjupiter-a2/last50_pres_hotjupiter.nc"

# time slices for averaging
#timeslices = range(50)     # average 0-49
timeslices = range(35,50)  # average 35-49
#timeslices = [5]           # Instantaneous at 5
print("time slices averaged:", list(timeslices))

# index of the vertical surface at a latitude 
lat_idx = 45   # 91 latitudes in total, 45 is the equator
"""============"""

# read the data
dataset = Dataset(data_path,'r')
time = dataset.variables['time'][:]
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
press = dataset.variables['press'][:]
temp = dataset.variables['temp'][:]
vlat = dataset.variables['vlat'][:]

# average data over time
mean_temp = np.mean(temp[timeslices,:,lat_idx,:], axis=0)
mean_vlat = np.mean(vlat[timeslices,:,lat_idx,:], axis=0)

# Downsample the data (one wind arrow per 3 longitude grid, per 3 pressure grid)
stride = 3
lon_downsampled = lon[::stride]
press_downsampled = press[::stride]
vlat_downsampled = mean_vlat[::stride, ::stride]


# Create a figure and a set of subplots
plt.figure(figsize=(10, 6))

# Plot temperature contour
#contour = plt.contour(lon, press, mean_temp)
#plt.colorbar(contour, label='Temperature')

# Plot temperature pcolor
pc = plt.pcolor(lon, press, mean_temp, shading='auto')
plt.colorbar(pc, label='Temperature')

# Plot the wind vectors with quiver
plt.quiver(lon_downsampled, press_downsampled, vlat_downsampled, 0, width=0.001, headwidth=5,headlength=5)

# Labels and title
plt.xlabel('Longitude')
plt.ylabel('Pressure')
plt.title(f'Temperature and horizontal winds at the equator')

# Set the y-axis to a logarithmic scale
plt.yscale('log')

# Inverting the y-axis if pressure increases with depth/altitude
plt.gca().invert_yaxis()

# Adjust the ticks on the y-axis to be more readable
#plt.yticks(pressure, labels=np.round(np.log10(pressure), 2))

plt.savefig("temp_wind_equator.png",dpi=300)
