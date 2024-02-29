from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np


"""User specify"""
# The data should be in time-press-lat-lon dimensions
data_path = "/home/linfel/data/hot_jupiter/hotjupiter-a2/last50_pres_hotjupiter.nc"

# time slices for averaging
timeslices = range(50)     # average 0-49
#timeslices = range(35,50)  # average 35-49
#timeslices = [5]           # Instantaneous at 5

# pressure level index
press_idx = 60   # 181 pressure levels in total
"""============"""

# read the data
dataset = Dataset(data_path,'r')
time = dataset.variables['time'][:]
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
press = dataset.variables['press'][:]
temp = dataset.variables['temp'][:]
vlat = dataset.variables['vlat'][:]
vlon = dataset.variables['vlon'][:]

# average data over time
data = np.mean(temp[timeslices,press_idx,:,:], axis=0)
mean_vlat = np.mean(vlat[timeslices,press_idx,:,:], axis=0)
mean_vlon = np.mean(vlon[timeslices,press_idx,:,:], axis=0)
print("# time slices averaged:", list(timeslices))

# Downsample the data
stride = 3 
lon_downsampled = lon[::stride]
lat_downsampled = lat[::stride]
vlat_downsampled = mean_vlat[::stride, ::stride]
vlon_downsampled = mean_vlon[::stride, ::stride]

# Calculate wind speed
wind_speed = (mean_vlat**2 + mean_vlon**2)**0.5

# get the pressure at the isobaric plane
pressure = press[press_idx]
print(f"# isobaric plane {pressure} Pa")

# Create a figure and a set of subplots
plt.figure(figsize=(10, 6))

# Plot contour
#contour = plt.contour(lon, lat, data)
#plt.colorbar(contour, label='Temperature')

# Plot pcolor
pc = plt.pcolor(lon, lat, data, shading='auto')
plt.colorbar(pc, label='Temperature')

# Plot the wind vectors with quiver
plt.quiver(lon_downsampled, lat_downsampled, vlat_downsampled, vlon_downsampled, width=0.001, headwidth=5, headlength=7)

# Labels and title
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title(f'Temperature and horizontal winds on the {pressure} Pa isobaric plane')

plt.savefig(f"temp_wind_{pressure}Pa_plane.png",dpi=300)
