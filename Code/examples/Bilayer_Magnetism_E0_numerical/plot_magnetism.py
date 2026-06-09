import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import sys
import os

# ============ for latex fonts ============
from matplotlib import rc #, font_manager
rc('text.latex', preamble=r'\usepackage{lmodern}')# this helps use the plots in tex files
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'xtick.labelsize': 14,
		  'ytick.labelsize': 14,
		  'xtick.major.pad': 6,
		  'ytick.major.pad': 6,
          'axes.titlesize': 14,
		  'font.serif': 'Computer Modern Roman',
		  'axes.formatter.use_mathtext': True,
		  'axes.labelpad': 6.0 })
# ==========================================

file = "results/n_map.csv"
if(len(sys.argv) > 1):
    file = sys.argv[1]

def read_csv(path: str, delimiter = ",") -> np.ndarray:
    data = pd.read_csv(path, delimiter = delimiter, header=None)
    data = np.array(data)
    if(data.shape[0] == 1):
        data = data.flatten()
    return data

def plot(fig, ax, n, label):
    # Tworzenie siatki
    grid_x, grid_y = np.meshgrid(
        np.linspace(min(E_vals), max(E_vals), 200),
        np.linspace(min(B_vals), max(B_vals), 200)
    )
    # Interpolacja do siatki
    grid_n = griddata(
        (E_vals, B_vals), n,
        (grid_x, grid_y),
        method='cubic'
    )
    # Kontury wypełnione
    contour = ax.contourf(grid_x, grid_y, grid_n, levels=50, cmap='RdBu_r')
    # Izolinie
    ax.contour(grid_x, grid_y, grid_n, levels=50, colors='black', linewidths=0.5, alpha=0.5)
    # # Oryginalne punkty (dla porównania)
    # plt.scatter(E_vals, B_vals, c=n_vals, cmap='viridis', edgecolors='k', s=30)
    # Pasek kolorów
    cbar = fig.colorbar(contour)
    cbar.set_label(label)
    # Opisy osi
    ax.set_xlabel("E0 [eV]")
    ax.set_ylabel("B [T]")

def plot1(fig, ax, n, label):
    ax.plot(E_vals, n)
    ax.plot(E_vals, np.zeros_like(n) + 3)
    ax.set_xlabel("$E0$ [eV]")
    ax.set_ylabel(label)
    ax.grid()

data = read_csv(file)

n_vals = data[0] * 1e-11
B_vals = data[1]
E_vals = data[2]
n0_vals = data[3] * 1e-11


if len(np.unique(B_vals)) != 1:
    fig, ax = plt.subplots(1,2, figsize=(14, 5))
    print("n")
    plot(fig, ax[0], n_vals, r"n [$10^{11} \text{cm}^{-2}]$")
    print("n0")
    plot(fig, ax[1], n0_vals, r"n0 [$10^{11} \text{cm}^{-2}]$")

    plt.tight_layout()
    fig.savefig("results/n_map_izolinie.png", dpi = 350)
else:
    fig, ax = plt.subplots(1,2, figsize=(14, 5))
    plot1(fig, ax[0], n_vals, r"n [$10^{11} \text{cm}^{-2}]$")
    plot1(fig, ax[1], n0_vals, r"n0 [$10^{11} \text{cm}^{-2}]$")
    plt.tight_layout()
    print()
    fig.savefig(
        file.replace(".csv", ".png"),
        dpi = 350)
