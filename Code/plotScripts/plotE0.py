import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

eV2au = 0.03674932587122423
au2eV = 1/eV2au
if len(sys.argv) != 4:
    print(f"python {sys.argv[0]} dir B Vb")
    exit(1)

directory = sys.argv[1]
E0t = np.loadtxt(f"{directory}/E0t.csv", delimiter=',') * au2eV
E0b = np.loadtxt(f"{directory}/E0b.csv", delimiter=',') * au2eV
Vgt = np.loadtxt(f"{directory}/Vgt.csv", delimiter=',') * au2eV
Vgb = np.loadtxt(f"{directory}/Vgb.csv", delimiter=',') * au2eV
T   = np.loadtxt(f"{directory}/T.csv",   delimiter=',')
B   = np.loadtxt(f"{directory}/B.csv",   delimiter=',')
Vb  = np.loadtxt(f"{directory}/Vb.csv",  delimiter=',')


fig, ax = plt.subplots(3,2, figsize=(14, 14))

# cross B
Bmin = np.min(B)
Bmax = np.max(B)
crossB = float(sys.argv[2])
print(f"B {crossB}")
indexB = int(len(B) * (crossB - Bmin) / (Bmax - Bmin))
ax[0, 0].plot(Vb, -Vgb[indexB, :], label = f"B = {crossB}")
ax[0, 0].legend()
ax[0, 0].set_ylabel("-Vgb")
ax[0, 0].set_xlabel("Vb")

ax[1, 0].plot(Vb, -E0b[indexB, :] - Vgb[indexB, :], label = f"B = {crossB}")
ax[1, 0].legend()
ax[1, 0].set_ylabel("-E0b -Vgb")
ax[1, 0].set_xlabel("Vb")

ax[2, 0].plot(Vb, T[indexB, :], label = f"B = {crossB}")
ax[2, 0].legend()
ax[2, 0].set_ylabel("T")
ax[2, 0].set_xlabel("Vb")
# cross Vb
Vbmin = np.min(Vb)
Vbmax = np.max(Vb)
crossVb = float(sys.argv[3])
print(f"Vb {crossVb}")
indexVb = int(len(Vb) * (crossVb - Vbmin) / (Vbmax - Vbmin))
ax[0, 1].plot(B, -Vgb[:, indexVb], label = f"Vb = {crossVb}")
ax[0, 1].legend()
ax[0, 1].set_ylabel("- Vgb")
ax[0, 1].set_xlabel("B")

ax[1, 1].plot(B, -E0b[:, indexVb] -Vgb[:, indexVb], label = f"Vb = {crossVb}")
ax[1, 1].legend()
ax[1, 1].set_ylabel("-E0b - Vgb")
ax[1, 1].set_xlabel("B")

ax[2, 1].plot(B, T[:, indexVb], label = f"B = {crossB}")
ax[2, 1].legend()
ax[2, 1].set_ylabel("T")
ax[2, 1].set_xlabel("Vb")

fig.tight_layout()
fig.savefig("CrossVgb.png")
