import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
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

def read_csv(path: str, delimiter = ',') -> np.ndarray:
    data = pd.read_csv(path, delimiter = delimiter, header=None)
    data = np.array(data)
    if(data.shape[0] == 1):
        data = data.flatten()
    return data


def getFiles(dir, ext = None) -> tuple[list[str], list[str]]:
    absPath = os.path.abspath(dir)
    files = []
    dirs = []

    fileOrFolderList = os.listdir(absPath)
    print(f"files in folder {dir}: {fileOrFolderList}")
    for fileOrFolder in fileOrFolderList:
        abs = os.path.join(absPath, fileOrFolder)
        if os.path.isfile(abs):
            files.append(abs)
        else:
            dirs.append(abs)

    if (ext is not None):
        files = [f for f in files if f.split('.')[-1] == ext]

    files.sort()

    return files, dirs

def getNameOfFile(file: str) -> str:
    nameWithExtension = file.split('/')[-1]
    splitDot = nameWithExtension.split('.')

    l = len(splitDot)
    if (l == 1 or l == 2):
        return splitDot[0]
    elif (len(splitDot) > 2):
        name = splitDot[0]
        for part in splitDot[1:-1]:
            name += '.'
            name += part
        return name

    print(f"Error in getNameOfFile")
    assert False


def assert_mess(cond: bool, mess:str = "") -> None:
    if not cond:
        print(mess)

    assert(cond)

################################################################################

def plotIm(data, B, Vb, title, nameSave):
    fig, ax = plt.subplots()

    im = ax.imshow(data, origin="lower", extent=(Vb[0], Vb[-1], B[0], B[-1]))
    ax.set_xlabel("Vb [V]")
    ax.set_ylabel("B [T]")
    ax.set_title(title)
    fig.colorbar(im)

    fig.tight_layout()
    fig.savefig(nameSave)

    return fig, ax


def filterFiles(files: list[str], containList: list[str]) -> list[str]:
    newFiles = [f for f in files if (f.split('/')[-1] in containList)]
    return newFiles


def plotSingle(files: list[str], dir: str) -> None:
    for f in files:

        data = read_csv(f)
        B = read_csv(f"{dir}/B.csv")
        Vb = read_csv(f"{dir}/Vb.csv")
        print(data)

        title = f"Vt = {round(float(f.split('Vt_')[-1].split('/')[0]))} [V], {getNameOfFile(f)}"
        fig, ax = plotIm(data, B, Vb, title, f.replace(".csv", ".pdf"))


def plotDiff(dir):
    Vgt = read_csv(f"{dir}/Vgt.csv")
    Vgb = read_csv(f"{dir}/Vgb.csv")
    diff = Vgt - Vgb

    B = read_csv(f"{dir}/B.csv")
    Vb = read_csv(f"{dir}/Vb.csv")

    title = f"Vgt - Vgb, Vt = {round(float(dir.split('Vt_')[-1].split('/')[0]))} [V]"
    plotIm(diff, B, Vb, title, f"{dir}/diff.pdf")

################################################################################

dir = "./results"
if (len(sys.argv) > 1):
    dir = sys.argv[1]

files, _ = getFiles(dir, "csv")

plotSingle(filterFiles(files, ["Vgb.csv", "Vgt.csv"]), dir)

plotDiff(dir)