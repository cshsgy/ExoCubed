import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp1d
from tqdm import tqdm

"""User define here"""
# Original input NetCDF file
inputfile = '/home/linfel/data/hjupiter/polar_xiz-0225-shj-main.nc' 
# New output NetCDF file
outputfile = '/home/linfel/data/hjupiter/pres_xiz-0225-shj-main.nc'

# 'new_press_levels' is provided as a 1D numpy array
new_log_press = np.linspace(7.4, 0., 149)
new_press_levels = 10.**new_log_press
"""================"""

# Open the original input NetCDF file
with nc.Dataset(inputfile, 'r') as src:
    # Read the dimensions
    t_dim = src.dimensions['time']
    x2_dim = src.dimensions['lat']
    x3_dim = src.dimensions['lon']
    x1 = src.variables['x1'][:]  # Original vertical coordinates
    press_var = src.variables['press'][:]  # Original pressure variable

    
    # Create a new output NetCDF file
    with nc.Dataset(outputfile, 'w') as dst:
        # Copy dimensions from the source to destination, except for 'x1' which is replaced by 'press'
        dst.createDimension('time', len(t_dim))
        dst.createDimension('lat', len(x2_dim))
        dst.createDimension('lon', len(x3_dim))
        dst.createDimension('press', len(new_press_levels))

        # Copy all variables from the source to destination, except those that depend on 'x1'
        for name, variable in src.variables.items():
            if 'x1' not in variable.dimensions:
                dst.createVariable(name, variable.datatype, variable.dimensions)
                dst.variables[name][:] = src.variables[name][:]

        # Define the new pressure coordinate variable
        new_press_var = dst.createVariable('press', new_press_levels.dtype, ('press',))
        new_press_var[:] = new_press_levels

        # Create the new variables that depend on 'x1'
        for name, variable in src.variables.items():
            print(name,flush=True)
            if name == 'time' or name == 'lat' or name == 'lon' or name == 'press' or name == 'x1':
                continue
            if 'x1' in variable.dimensions:
                # Create the new variable in the destination file
                new_dimensions = tuple('press' if dim == 'x1' else dim for dim in variable.dimensions)
                dst.createVariable(name, variable.datatype, new_dimensions)
               

    # Interpolate the variables that depend on 'x1'
    # Loop over time slices
    for t in tqdm(range(len(t_dim)), desc="Processing time steps"):
        print(t,'......',flush=True)
        for name, variable in src.variables.items():
            if name == 'time' or name == 'lat' or name == 'lon' or name == 'press' or name == 'x1':
                continue
            if 'x1' in variable.dimensions:
                print("Interpolating variable: " + name, flush=True)
                for x2 in range(len(x2_dim)):
                    for x3 in range(len(x3_dim)):
                        # Extract the slice of pressure values for the current point
                        original_log_press = np.log10(press_var[t, :, x2, x3])
                        # Extract the original data slice
                        original_data = variable[t, :, x2, x3]
                        # Create an interpolation function based on the original pressure and data
                        f = interp1d(original_log_press, original_data, kind='linear', bounds_error=False, fill_value="extrapolate")
                        # Interpolate to the new press levels
                        interp_data = f(new_log_press)
                        # Insert the interpolated data into the new NetCDF file
                        with nc.Dataset(outputfile, 'a') as dst: 
                            dst.variables[name][t, :, x2, x3] = interp_data
