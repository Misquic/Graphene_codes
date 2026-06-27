import numpy
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import colors, cm
from scipy.interpolate import griddata

from utils import *

cmap=cm.Purples

def getParamsFromDir(dir: str) -> tuple[float, float, float]:
    if dir[-1] == "/":
        dir = dir[:-1]
    baseDirName = os.path.basename(dir)
    split = baseDirName.split('_')
    B = float(split[1])
    Vb = float(split[3])
    Vt = float(split[5])

    return B, Vb, Vt

directory = "./results/"
if len(sys.argv) >= 2:
    directory = sys.argv[1]

outDir = directory
if len(sys.argv) >= 3:
    outDir = sys.argv[2]

absPath = os.path.abspath(directory)
sf = int(absPath.split("sf")[-1].split("/")[0])

files, _ = getFiles(directory, ext = "txt")
files = [f for f in files if "current" in getNameOfFile(f)]
print(files)

fig, axes = plt.subplots(2, len(files), figsize=(len(files) * 6, 12))
# try:
B, Vb, _ = getParamsFromDir(directory)
# except:
    # B = 0
    # Vb = 0

fig.suptitle(f"B = {B} T, Vb = {float(Vb)} V", fontsize = 16)
for file, ax0, ax1 in zip(files, axes[0,:], axes[1,:]):
    print(file)
    data  = np.loadtxt(file, usecols=(0, 1, 3, 4, 5))

    y1 = data[:,1] / nm2au
    x1 = data[:,0] / nm2au

    xmin = np.min(x1)
    xmax = np.max(x1)
    ymin = np.min(y1)
    ymax = np.max(y1)

    xi = np.unique(x1)
    yi = np.unique(y1)
    gx,gy = np.meshgrid(xi,yi)
    nx = len(xi)
    ny = len(yi)

    u = data[:,2] # x component
    v = data[:,3] # y component
    current = np.sqrt(u**2 + v**2) # current density

    gridr = griddata((x1, y1), current, (gx,gy), method="linear", fill_value=np.nan)
    gridu = griddata((x1, y1), u, (gx,gy), method="linear", fill_value=np.nan)
    gridv = griddata((x1, y1), v, (gx,gy), method="linear", fill_value=np.nan)

    gridu = gridu / gridr
    gridv = gridv / gridr

    # wektory pradu (co 4-ty zeby strzalki nie za gesto, mozna dobrac do rysunku)
    skip = 8
    ax0.quiver(xi[1::skip],
               yi[1::skip],
               gridu[1::skip,
               1::skip],
               gridv[1::skip,
               1::skip],
               color='k',
               scale = 30)
    ax0.set_title(getNameOfFile(file).replace("current_", "lead "))
    # mapa gestosci
    im = ax1.scatter(data[::skip,0,] / nm2au,
                     data[::skip,1] / nm2au,
                     c=current[::skip],
                     s=skip*1.2*(sf/8),
                     alpha = 1,
                     marker="H",
                     edgecolors=None,
                     linewidths=0,
                     cmap=cmap)
    # im = ax1.imshow(gridr, extent=(xi[0], xi[nx-1], yi[0], yi[ny-1]), origin='lower', interpolation='none', cmap=cmap)#, vmin=0, vmax=vmax)
    ax0.set_xlabel("x [nm]")
    ax1.set_xlabel("x [nm]")
    ax0.set_ylabel("y [nm]")
    ax1.set_ylabel("y [nm]")
    cbar = fig.colorbar(im, orientation = 'vertical')
    cbar.set_label("$|J|$")

plt.tight_layout()
print("saving")
fig.savefig(os.path.join(outDir, f"currents_B_{B}_Vb_{Vb}.png"), dpi = 400)
