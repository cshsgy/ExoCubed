# This util is used to transform the affine coordinate netcdf output 
# to be compared with the cartesian outputs
#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse, os, shutil
import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp2d
from tqdm import tqdm

fig = plt.figure()

# Set constants
transform_angle = np.pi/3
n = 500 # Final resolution

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
    default = '../build/bin/affine_results/sw.out2.00000.nc',
    help = 'file to be transformed and plotted'
    )
parser.add_argument('-c', '--compare',
    default = '../build/bin/sw.out2.00000.nc',
    help = 'file to be compared with'
    )
args = vars(parser.parse_args())
number = args['input'][-8:-3]
tmp_fign = 'tmp_affine_to_cart_' + number + '.png'
print('Saving to ')
print(tmp_fign)

# Input file
fn = args['input']
ds = nc.Dataset(fn)
# print(ds['rho'].shape)
# print(ds['x1'].shape)
# rho shape: (1, 2, 1024, 512)
# x1 shape: (2,)
x2 = ds['x2']
x3 = ds['x3']
y = np.transpose(np.squeeze(ds['rho'][0,0,:,:]))

# setup interpolation
f = interp2d(x2,x3,y,kind='linear')

# Calculate the limits of the de-transformed field, setup coordinates
x2t_lim = np.max(x2)+np.max(x3)*np.cos(transform_angle)
x3t_lim = np.max(x3)*np.sin(transform_angle)
x2t = np.linspace(-x2t_lim,x2t_lim,n)
x3t = np.linspace(-x3t_lim,x3t_lim,n)
rhot = np.zeros((n,n))

print("Deprojecting the affine coordinates...")

for i in tqdm(range(n)):
    for j in range(n):
        # Calculate the affine coord
        x2_now = x2t[i]-x3t[j]/np.tan(transform_angle)
        x3_now = x3t[j]/np.sin(transform_angle)
        # Do interpolation
        rhot[j,i] = f(x2_now,x3_now)

ax = fig.add_subplot(121)
ax.contourf(x2t,x3t,rhot)
ax.set_xlim([-6,6])
ax.set_ylim([-6,6])
ax.set_clim([10,12])
ax.set_aspect(1)

# Compare file
fn = args['compare']
ds = nc.Dataset(fn)
# print(ds['rho'].shape)
# print(ds['x1'].shape)
# rho shape: (1, 2, 1024, 512)
# x1 shape: (2,)
x2 = ds['x2']
x3 = ds['x3']
y = np.transpose(np.squeeze(ds['rho'][0,0,:,:]))

# setup interpolation
f = interp2d(x2,x3,y,kind='linear')

# Calculate the limits of the de-transformed field, setup coordinates
x2t_lim = np.max(x2)
x3t_lim = np.max(x3)
x2t = np.linspace(-x2t_lim,x2t_lim,n)
x3t = np.linspace(-x3t_lim,x3t_lim,n)
rhot = np.zeros((n,n))

print("Resampling the cartesian coordinates...")

for i in tqdm(range(n)):
    for j in range(n):
        # Calculate the affine coord
        x2_now = x2t[i]
        x3_now = x3t[j]
        # Do interpolation
        rhot[j,i] = f(x2_now,x3_now)

ax = fig.add_subplot(122)
ax.contourf(x2t,x3t,rhot)
ax.set_xlim([-6,6])
ax.set_ylim([-6,6])
ax.set_clim([10,12])
ax.set_aspect(1)

fig.savefig(tmp_fign)