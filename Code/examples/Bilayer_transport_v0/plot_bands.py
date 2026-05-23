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

file = dir + "bands.dat"

# data = np.loadtxt(file)
data = np.loadtxt(file)
print(data)
no_lines = np.size(data[0,:]) - 1
print(no_lines)
x = data[:,0]

fig, ax = plt.subplots()
for i in range(no_lines):
    ax.plot(x, data[:,i+1],ls='-')

ax.set_xlabel("k [1/unit size]")
ax.set_ylabel("Energy [some units]")

plt.savefig(dir + "bands.png")
