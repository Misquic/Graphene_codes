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

