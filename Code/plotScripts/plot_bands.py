#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from utils import *
from matplotlib.ticker import AutoMinorLocator

if (len(sys.argv) > 1):
    dir = sys.argv[1]
else:
    dir = "./results/"

# file = dir + "bands1.dat"

files, _ = getFiles(dir, "dat")

for f in files:
    if not "bands" in f:
        continue
    # data = np.loadtxt(file)
    data = np.loadtxt(f)
    print(data.shape)
    no_lines = np.size(data[0,:]) - 1
    print(no_lines)
    x = data[:,0]

    fig, ax = plt.subplots()
    for i in range(no_lines):
        ax.plot(x, data[:,i+1],ls='-')

    ax.set_xlabel("k [1/unit size]")
    ax.set_ylabel("Energy [some units]")

    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid()
    ax.grid(which="minor")

    B, Vb, _ = getParamsFromDir(dir)
    name = os.path.basename(f).split('.')[0].strip(' ')
    ax.set_title(f"{name} B={B} Vb={Vb}")

    plt.tight_layout()
    plt.savefig(os.path.join(dir, name + ".pdf"))
