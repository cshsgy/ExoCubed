from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

"""User specify"""
# read the data
data_path = "/home/linfel/data/hot_jupiter/hotjupiter-a2/last50_pres_hotjupiter.nc"

# time slices for averaging
timeslices = range(50)     # average 0-49
#timeslices = range(35,50)  # average 35-49
#timeslices = [5]           # Instantaneous at 5

# index of the vertical plane at a latitude
lat_idx = 45   # latitude index 0-90, equator is 45
"""============"""

# read the data
dataset = Dataset(data_path,'r')
time = dataset.variables['time'][:]
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
press = dataset.variables['press'][:]
temp = dataset.variables['temp'][:]

# average data over time
data = np.mean(temp[timeslices,:,lat_idx,:], axis=0)
print("# time slices averaged:", list(timeslices))

# get the latitude of the vertical plane
latitude = lat[lat_idx]
print(f"# latitude {latitude}° vertical plane")

# Create a figure and a set of subplots
plt.figure(figsize=(10, 6))

# Plot contour
contour = plt.contour(lon, press, data)
plt.colorbar(contour, label='Temperature')

# Plot pcolor
#pc = plt.pcolor(lon, press, data, shading='auto')
#plt.colorbar(pc, label='Temperature')

# Labels and title
plt.xlabel('Longitude')
plt.ylabel('Pressure / Pa')
if latitude == 0.:
	plt.title('Temperature Profile at Equator')
else:
	plt.title(f'Temperature Profile at latitude {latitude}°')

# Set the y-axis to a logarithmic scale
plt.yscale('log')

# Inverting the y-axis if pressure increases with depth/altitude
plt.gca().invert_yaxis()

# Adjust the ticks on the y-axis to be more readable
#plt.yticks(pressure, labels=np.round(np.log10(pressure), 2))

plt.savefig("images/T-P_profile_equator.png",dpi=300)
