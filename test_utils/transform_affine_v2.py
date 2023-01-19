#! /usr/bin/env python3.6
import matplotlib as mpl
#mpl.use('Agg')
from pylab import *
from matplotlib import pyplot as plt
import argparse, os, shutil
import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp2d
from tqdm import tqdm

fig = plt.figure()
fig.set_size_inches(14.5, 6.5)

# Set constants
transform_angle = np.pi/3
#transform_angle = np.pi/2
n = 200 # Final resolution

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
    default = '../build_affine/bin/sw-log13-main.nc',
    help = 'file to be transformed and plotted'
    )
parser.add_argument('-c', '--compare',
    default = '../build_cart/bin/sw-log3-main.nc',
    help = 'file to be compared with'
    )
parser.add_argument('-d', '--old',
    default = '../build_old/bin/sw-cart-main.nc',
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
y = np.transpose(np.squeeze(ds['rho'][-1,0,:,:]))

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

ax = fig.add_subplot(131)
ax.contourf(x2t,x3t,rhot,cmap='Reds')
ax.set_xlim([-4,4])
ax.set_ylim([-4,4])
ax.set_xlabel("x_2")
ax.set_ylabel("x_3")
ax.set_aspect(1)

ax3 = fig.add_subplot(133)
ax3.plot(x2t, rhot[int(rhot.shape[0]/2),:])

# Compare file
fn = args['compare']
ds = nc.Dataset(fn)
# print(ds['rho'].shape)
# print(ds['x1'].shape)
# rho shape: (1, 2, 1024, 512)
# x1 shape: (2,)
x2 = ds['x2']
x3 = ds['x3']
y = np.transpose(np.squeeze(ds['rho'][-1,0,:,:]))

# old file
fn = args['old']
ds = nc.Dataset(fn)
x2 = ds['x1']
x3 = ds['x2']
z = np.transpose(np.squeeze(ds['rho'][-1,:,:,0]))

# setup interpolation
#f = interp2d(x2,x3,y,kind='linear')
f = interp2d(x2,x3,z,kind='linear')

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

ax = fig.add_subplot(132)
ax.contourf(x2t,x3t,rhot,cmap='Reds')
ax.set_xlim([-4,4])
ax.set_ylim([-4,4])
ax.set_xlabel("x_2")
ax.set_ylabel("x_3")
ax.set_aspect(1)

ax3.plot(x2t, rhot[int(rhot.shape[0]/2),:])
ax3.set_xlim([-4,4])
ax3.set_xlabel("x_2")
ax3.set_ylabel("H")
ax3.set_aspect(4)

show()

#fig.savefig(tmp_fign)
