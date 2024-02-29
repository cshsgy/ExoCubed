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

"""============"""

# read the data
dataset = Dataset(data_path,'r')
lat = dataset.variables['lat'][:]
pressure = dataset.variables['press'][:]
vlat = dataset.variables['vlat'][:]
vlon = dataset.variables['vlon'][:]

# take the average
data = vlat
time_mean_data = np.mean(data[timeslices,:,:,:], axis=0)
zonal_mean_data = np.mean(time_mean_data, axis=2)

# Create a figure and a set of subplots
plt.figure(figsize=(10, 6))

# Create a pseudocolor plot with a non-regular rectangular grid
c = plt.pcolor(lat, pressure, zonal_mean_data, shading='auto', cmap='RdBu_r', vmin=-np.max(np.abs(data)), vmax=np.max(np.abs(data)))

# Add a colorbar to show the temperature scale
plt.colorbar(c, label='Wind speed (m/s)')

# Labels and title
plt.xlabel('Latitude')
plt.ylabel('Pressure')
plt.title('Zonal-mean zonal wind')

# Set the y-axis to a logarithmic scale
plt.yscale('log')

# Inverting the y-axis if pressure increases with depth/altitude
plt.gca().invert_yaxis()

# Adjust the ticks on the y-axis to be more readable
#plt.yticks(pressure, labels=np.round(np.log10(pressure), 2))

plt.savefig("zonal-mean_zonal_wind.png",dpi=300)
