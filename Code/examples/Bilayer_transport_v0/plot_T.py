#!/usr/bin/python
"""
Created on Thu Mar  5 14:16:21 2015

@author: Krzysztof Kolasinski
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

if (len(sys.argv) > 1):
    dir = sys.argv[1]
else:
    dir = "./results/"

file = dir + "T.dat"

data = np.loadtxt(file, delimiter=',')
no_lines = np.size(data[0,:])
x = data[:,0]
plt.plot(x,data[:,1],c='k',ls='-')

plt.xlabel("Ef [energy units]")
plt.ylabel("Transmission")
plt.savefig(dir + "T.png")