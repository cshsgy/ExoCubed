from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

"""User specify"""
# read the data
data_path = "/home/linfel/data/hjupiter/pres_xiz-0225-shj-main.nc"

# time slices for averaging
timeslices = range(54)     # average 0-53
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

# get the latitude
latitude = lat[lat_idx]
print(f"# latitude {latitude}°")

# Create a figure
plt.figure(figsize=(10, 6))

for ilon in range(0,91,6):
	plt.plot(data[:,ilon],press,label=str(lon[ilon])+'°')

# Plot contour
#contour = plt.contour(lon, press, data)
#plt.colorbar(contour, label='Temperature')

# Plot pcolor
#pc = plt.pcolor(lon, press, data, shading='auto')
#plt.colorbar(pc, label='Temperature')

# Labels and title
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')

if latitude == 0.:
	plt.title('T-P Profiles at Equator')
else:
	plt.title(f'T-P Profiles at latitude {latitude}°')

plt.legend()

# Set the y-axis to a logarithmic scale
plt.yscale('log')

# Inverting the y-axis if pressure increases with depth/altitude
plt.gca().invert_yaxis()

# Adjust the ticks on the y-axis to be more readable
#plt.yticks(pressure, labels=np.round(np.log10(pressure), 2))

plt.savefig("T-P_profile_equator.png",dpi=300)
