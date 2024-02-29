import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp1d
from tqdm import tqdm


def generate_uneven_seq(low_limit,up_limit,step_order):
    """ example:
        generate_uneven_seq(0,2,1):   [1, 1.1, 1.2, 1.3, 1.4,..., 9.9, 10, 11, 12, 13, 14,..., 99, 100]
        generate_uneven_seq(-1,1,2):  [0.1, 0.101, 0.102, 0.103,..., 0.999, 1, 1.01, 1.02, 1.03,..., 9.99, 10]
    """
    i = low_limit
    seq = []
    while i<up_limit:
        seq_i = list(np.arange(10.**i, 10.**(i+1), 10.**(i-step_order)))
        seq = seq + seq_i
        i += 1
    seq.append(10.**up_limit)
    return np.array(seq)

"""User define here"""
# Original input NetCDF file
inputfile = '/home/linfel/data/hot_jupiter/hotjupiter-a2/last50_polar_hotjupiter-a2-main.nc' 
# New output NetCDF file
outputfile = '/home/linfel/data/hot_jupiter/hotjupiter-a2/last50_pres_hotjupiter.nc'

# 'new_press_levels' is provided as a 1D numpy array of the new pressure levels
#new_press_levels = np.linspace(1E5, 0.01E5, 100)
new_press_levels = np.flip(generate_uneven_seq(3,5,1))
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
                        press_slice = press_var[t, :, x2, x3]
                        # Extract the original data slice
                        original_data = variable[t, :, x2, x3]
                        # Create an interpolation function based on the original pressure and data
                        f = interp1d(press_slice, original_data, kind='linear', bounds_error=False, fill_value="extrapolate")
                        # Interpolate to the new press levels
                        interp_data = f(new_press_levels)
                        # Insert the interpolated data into the new NetCDF file
                        with nc.Dataset(outputfile, 'a') as dst: 
                            dst.variables[name][t, :, x2, x3] = interp_data
