import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
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

# dG = 0.34 nm
dG = 0.34
if(len(sys.argv) > 2):
    dG = float(sys.argv[2])
print("dG = " + str(dG))

nt_dir = "results/bilayer_nt/"
nb_dir = "results/bilayer_nb/"
nt = read_csv(nt_dir + "dG_" + str(dG) + ".csv")
nb = read_csv(nb_dir + "dG_" + str(dG) + ".csv")

Vt_min = -6
Vt_max = 6
Vb_min = -10
Vb_max = 10

def dens(arr, cmap, dir, title):
    fig, ax = plt.subplots()
    im = ax.imshow(arr, cmap = cmap, extent=(Vt_min, Vt_max, Vb_min, Vb_max), aspect="auto")
    ax.set_title(title)
    ax.set_xlabel("$V_t$ [V]")
    ax.set_ylabel("$V_b$ [V]")
    fig.colorbar(im, ax = ax)
    plt.tight_layout()
    fig.savefig(dir + "dG_" + str(dG) + ".png")

    if SHOW:
        plt.show()
    
dens(nt, cmap = "RdBu", dir = nt_dir, title = "$dG = " + str(dG) + "$ $n_t$")
dens(nb, cmap = "RdBu", dir = nb_dir, title = "$dG = " + str(dG) + "$ $n_b$")

nt_cut = np.where(np.abs(nt) < 1e8, 1, 0)
nb_cut = np.where(np.abs(nb) < 1e8, 2, 0)

dens(nt_cut + nb_cut, cmap = "viridis", dir = "results/zeros/", title = "$dG = " + str(dG) + "$ zeros, nt = 1, nb = 2")
