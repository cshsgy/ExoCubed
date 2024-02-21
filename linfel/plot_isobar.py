from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

# read the data
data_path = "/home/linfel/data/hot_jupiter/hotjupiter-a2/last50_pres_hotjupiter.nc"
dataset = Dataset(data_path,'r')
time = dataset.variables['time'][:]
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
press = dataset.variables['press'][:]
temp = dataset.variables['temp'][:]
vlat = dataset.variables['vlat'][:]
vlon = dataset.variables['vlon'][:]


# specify the time slices for averaging
timeslices = range(50)
# specify the index of the isobaric plane
press_idx = 150

# average data over time
data = np.mean(temp[timeslices,press_idx,:,:], axis=0)
mean_vlat = np.mean(vlat[timeslices,press_idx,:,:], axis=0)
mean_vlon = np.mean(vlon[timeslices,press_idx,:,:], axis=0)

# Downsample the data
stride = 3 
lon_downsampled = lon[::stride]
lat_downsampled = lat[::stride]
vlat_downsampled = mean_vlat[::stride, ::stride]
vlon_downsampled = mean_vlon[::stride, ::stride]

# Calculate wind speed
wind_speed = (mean_vlat**2 + mean_vlon**2)**0.5

# get the pressure at the isobaric plane
pressure = press[150]
print(f"# isobaric plane {pressure} Pa")

# Create a figure and a set of subplots
plt.figure(figsize=(10, 6))

# Plot contour
#contour = plt.contour(x, y, data)
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

# Set the y-axis to a logarithmic scale
#plt.yscale('log')

# Inverting the y-axis if pressure increases with depth/altitude
#plt.gca().invert_yaxis()

# Adjust the ticks on the y-axis to be more readable
#plt.yticks(pressure, labels=np.round(np.log10(pressure), 2))

plt.savefig("images/test7.png",dpi=300)
