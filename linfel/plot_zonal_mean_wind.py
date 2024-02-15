from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

# read the data
data_path = "linshi_averages_and_products.nc"
dataset = Dataset(data_path,'r')
lat = dataset.variables['lat'][:]
pressure = dataset.variables['pressure'][:]
vlon_avg = dataset.variables['vlon_avg'][:]


x = lat
y = pressure
data = vlon_avg
# Create a figure and a set of subplots
plt.figure(figsize=(10, 6))

# Create a pseudocolor plot with a non-regular rectangular grid
c = plt.pcolor(x, y, data, shading='auto', cmap='RdBu_r', vmin=-np.max(np.abs(data)), vmax=np.max(np.abs(data)))

# Add a colorbar to show the temperature scale
plt.colorbar(c, label='Wind speed')

# Labels and title
plt.xlabel('Latitude')
plt.ylabel('Pressure (log scale)')
plt.title('Zonal-mean meridional wind')

# Set the y-axis to a logarithmic scale
plt.yscale('log')

# Inverting the y-axis if pressure increases with depth/altitude
plt.gca().invert_yaxis()

# Adjust the ticks on the y-axis to be more readable
#plt.yticks(pressure, labels=np.round(np.log10(pressure), 2))

plt.savefig("images/meridwind_test1.png",dpi=300)
