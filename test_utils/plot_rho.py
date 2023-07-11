#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import glob
import re
import netCDF4 as nc
import matplotlib.pyplot as plt
from tqdm import tqdm

def main():
    file_pattern = 'cubed.out2.*.nc'
    files = sorted(glob.glob(file_pattern), key=lambda f: int(re.findall(r'cubed\.out2\.(\d+)\.nc', f)[0]))

    orders = []
    rho_averages = []

    for file in tqdm(files, desc='Processing files'):
        order = int(re.findall(r'cubed\.out2\.(\d+)\.nc', file)[0])
        rho_avg = compute_rho_average(file)

        orders.append(order)
        rho_averages.append(rho_avg)

    plot_rho_vs_order(orders, rho_averages)

def compute_rho_average(file):
    with nc.Dataset(file, 'r') as dataset:
        rho = dataset.variables['rho'][:]
        rho_average = rho.mean()

    return rho_average

def plot_rho_vs_order(orders, rho_averages):
    plt.plot(orders, rho_averages, marker='o', linestyle='-')
    plt.xlabel('Order')
    plt.ylabel('Average rho')
    plt.title('Average rho vs. Order')
    plt.grid()
    plt.savefig('rho.png', dpi=300)
    plt.close()

if __name__ == '__main__':
    main()
