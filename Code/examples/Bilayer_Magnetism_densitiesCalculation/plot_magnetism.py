import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import sys

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

SHOW = False
if(len(sys.argv) > 1):
    if(sys.argv[1] == str(0)):
        SHOW = False
    else:
        SHOW = bool(sys.argv[1])

def read_csv(path: str, delimiter = ",") -> np.ndarray:
    data = pd.read_csv(path, delimiter = delimiter, header=None)
    data = np.array(data)
    if(data.shape[0] == 1):
        data = data.flatten()
    if(SHOW):
        print(data)
    return data

data = read_csv("results/Magnetism/n_map.csv")

n_vals = data[0]
B_vals = data[1]
E_vals = data[2]

# ------------------------
# Tworzenie siatki
# ------------------------
grid_x, grid_y = np.meshgrid(
    np.linspace(min(E_vals), max(E_vals), 200),
    np.linspace(min(B_vals), max(B_vals), 200)
)

# Interpolacja do siatki
grid_n = griddata(
    (E_vals, B_vals), n_vals,
    (grid_x, grid_y),
    method='cubic'
)

fig, ax = plt.subplots()

# Kontury wypełnione
contour = ax.contourf(grid_x, grid_y, grid_n, levels=150, cmap='RdBu')

# Izolinie
ax.contour(grid_x, grid_y, grid_n, levels=150, colors='black', linewidths=0.5, alpha=0.5)

# Oryginalne punkty (dla porównania)
# plt.scatter(E_vals, B_vals, c=n_vals, cmap='viridis', edgecolors='k', s=30)

# Pasek kolorów
cbar = fig.colorbar(contour)
cbar.set_label("n")

# Opisy osi
ax.set_xlabel("E0 [eV]")
ax.set_ylabel("B [T]")

plt.tight_layout()
fig.savefig("results/Magnetism/n_map_izolinie.png", dpi = 350)
plt.show()