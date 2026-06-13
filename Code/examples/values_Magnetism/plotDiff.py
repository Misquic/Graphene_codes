#!/net/home/plgrid/plgkamilsocko/.conda/envs/normal/bin/python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

eV2au = 0.03674932587122423
au2eV = 1/eV2au

IND_VB = 0
IND_B = 1
IND_NT = 2
IND_NB = 3
IND_OT = 4
IND_OB = 5
IND_ET = 6
# IND_EB = 7

print("works")
dataF = np.loadtxt("outF.csv")
data = np.loadtxt("out.csv", delimiter=',')

print(dataF.shape)
print(data.shape)


B = np.unique(data[:, IND_B])
BF = np.unique(dataF[:, IND_B])
print(B  - BF)
VB = np.unique(data[:, IND_VB])
VBF = np.unique(dataF[:, IND_VB])
print(VB - VBF)


diffE0 = (data[:, IND_ET] - dataF[:, IND_ET]) / (dataF[:, IND_ET]) * 100
diffE0 = (data[:, IND_ET] - dataF[:, IND_ET]) / (dataF[:, IND_ET]) * 100

print(diffE0.shape)
diffE0 = diffE0.reshape((len(B), len(VB)))

print(diffE0.shape)

# diffOT = data[:, IND_OT] - dataF[:, IND_OT]
# diffOT = (data[:, IND_OB] - dataF[:, IND_OB]) / (dataF[:, IND_OB]) * 100
diffOT = (data[:, IND_OT] - dataF[:, IND_OT]) / (dataF[:, IND_OT]) * 100
diffOT = diffOT.reshape((len(B), len(VB)))

fig, ax = plt.subplots(1,2, figsize = (14, 8))
ax[0].plot(B, diffE0)
ax[0].set_xlabel("B [T]")
ax[0].set_ylabel("diff(E0) [%]")

ax[1].plot(VB, diffOT[::5].transpose(), label = B[::5])
ax[1].legend()
ax[1].set_xlabel("Vb [V]")
ax[1].set_ylabel("diff(onsite Top) [%]")

fig.tight_layout()
fig.savefig("diff.pdf")