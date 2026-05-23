'''
This file serves as a way tu run multiple simulations for different B, Vt, Vb and then plot results
It menages files and directories, use help to see usage
'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
import time
import asyncio
# from types import List
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

from utils import *

# ============ for latex fonts ============
from matplotlib import rc #, font_manager
rc("text.latex", preamble=r"\usepackage{lmodern}")# this helps use the plots in tex files
plt.rcParams.update({"font.size": 14})
plt.rcParams.update({"xtick.labelsize": 14,
		             "ytick.labelsize": 14,
		             "xtick.major.pad": 6,
		             "ytick.major.pad": 6,
                     "axes.titlesize": 14,
		             "font.serif": "Computer Modern Roman",
		             "axes.formatter.use_mathtext": True,
		             "axes.labelpad": 6.0 })
# ==========================================

# np co 0.2 V
args = {
    "BMin": 1., "BMax": 8.,   # min and max induction in T
    "VbMin": 1., "VbMax": 10., # min and max bottom gate voltages in V
    "VtMin": 0., "VtMax": 0., # min and max top gate voltages in V
    "maxParallel": 2, # maximum number of parralel simulations
    "numB": 5, "numVb": 5, "numVt":1, # number of B/Vt/Vb values to run
    "clearDir": 0, # if 1 then moves "allResultsDir" to "allResultsDir"Old
    "plotAll": 0, # if 1 then plots results after simulation
    "runSim": 0, # if 1 then runs simulation
    "allResultsDir": "./results/", # directory to store all results
    "saveStdout": 0, # if 1 then saves outputs from simulations to file else >dev/null
    "saveSystem": 0, # if 1 then saves created system to file
    "processFiles": 0, # if 1 then forces processing files even when previously processed
    "show": 0, # if 1 then shows plots at the end of plotting
    "runTransport": 1, # if 1 then runTransport = 1
    "sf": 8, # scaling factor
    }

def createTab(min, max, num = 3) -> np.ndarray:
    if min == max:
        return np.array([min])
    else:
        return np.linspace(min, max, num)

def execCommand(command: str) -> int:
    print(f"----------------------------------\n"
          f"Executing: {command}\n")
    result = os.system(command)
    print(f"Command executed with result={result}\n"
          f"----------------------------------\n")
    return result

def prepareCommandsAndDirs()-> list[str]:
    saveSystem = args["saveSystem"]
    runTransport = 1 if args["runTransport"] == 1 else 0
    runEnergyScan = 0
    plotResults = 0
    saveDensities = 0
    saveBands = 0
    sf = args["sf"]

    BTab = createTab(args["BMin"], args["BMax"], args["numB"])
    VbTab = createTab(args["VbMin"], args["VbMax"], args["numVb"])
    VtTab = createTab(args["VtMin"], args["VtMax"], args["numVt"])
    print(f"B = {BTab}")
    print(f"Vb = {VbTab}")
    print(f"Vt = {VtTab}")

    commands = []
    allResultsDir = args["allResultsDir"]

    if (args["clearDir"] == 1):
        newDir = allResultsDir[:-1] + "Old/"
        print(newDir)

        if os.path.exists(newDir):
            execCommand(f"rm -r {newDir}")
        execCommand(f"mkdir {newDir}")
        if os.path.exists(allResultsDir):
            execCommand(f"mv {allResultsDir}* {newDir}")
            execCommand(f"rm -r {allResultsDir}")
        execCommand(f"mkdir {allResultsDir}")
        execCommand(f"mkdir {allResultsDir}dirs/")

    # command to save system
    command = f"./Transport2D {allResultsDir} {np.max(BTab)} -60 0 1 0 0 0 0 0 {sf}"
    commands.append(command)

    # Vt is outside because it is likely to be single value
    for Vt in VtTab:
        for B in BTab:
            for Vb in VbTab:

                resultsDir = f"{allResultsDir}dirs/B_{B}_Vb_{Vb}_Vt_{Vt}/"
                if os.path.isdir(resultsDir) and os.path.exists(resultsDir):
                    execCommand(f"rm -r {resultsDir}") # clear resultsDir
                if not os.path.isdir(resultsDir):
                    os.makedirs(resultsDir)

                # "usage: ./Transport2D <resultsDir> <B in T> <Vb> <Vt> &
                #  <save_system> <run_transport> <plot_results> &
                #  <save_densities> <save_bands>"
                commandArgs = [ resultsDir,
                                str(B), str(Vb), str(Vt),
                                str(saveSystem), str(runTransport), str(runEnergyScan),
                                str(plotResults), str(saveDensities), str(saveBands), str(sf)]

                command = "./Transport2D "
                for a in commandArgs:
                    command += str(a) + " "
                # print(f"running with args: {commandArgs}")
                if (args["saveStdout"] == 1):
                    command += f"1>{resultsDir}stdout.txt "
                else:
                    command += f"1>/dev/null "
                command += f"2>{resultsDir}stderr.ansi "
                commands.append(command)

    with open(allResultsDir + "commands.txt", 'w') as f:
        for idx, command in enumerate(commands):
            print(f"[{idx}] {command}", file=f)

    return commands

async def runCommandAsync(command: str, id: int, maxId: int, timeStart: float) -> int:
    proc = await asyncio.create_subprocess_shell(command)
    await proc.wait()
    if (proc.returncode != 0):
        print(f"Command {id} failed, result = {proc.returncode}", flush=True)
    else:
        progressBar(id+1, 1, maxId, timeStart)
    return proc.returncode if proc.returncode is not None else 0

async def runCommandsInternal(commands: list[str], maxParallel: int, timeStart: float) -> list[tuple[int, int]]:
    semaphore = asyncio.Semaphore(maxParallel)
    maxId = len(commands)

    async def runWithLimit(command: str, id: int) -> tuple[int, int]:
        async with semaphore:
            returnCode = await runCommandAsync(command, id, maxId, timeStart)
            return returnCode, id

    results = await asyncio.gather(*[runWithLimit(cmd, id) for id, cmd in enumerate(commands)])
    return results

def runCommands(commands: list[str], maxParallel: int = 2) -> None:
    timeStart = time.time()
    print(f"starting {len(commands)} commands", flush=True)
    results = asyncio.run(runCommandsInternal(commands, maxParallel, timeStart))

    failed = [cmd_id for returncode, cmd_id in results if returncode != 0]
    if failed:
        print(f"{len(failed)} commands failed: {failed}")
    else:
        print(f"All {len(commands)} commands completed successfully")

def runSim() -> None:
    commands = prepareCommandsAndDirs()

    timeStart = time.time()
    runCommands(commands, maxParallel = args["maxParallel"]) # simulation itself uses multithreading
    timeEnd = time.time()
    totalTime = timeEnd - timeStart
    minutes = totalTime // 60
    print(f"All runs took {minutes} min {round(totalTime - minutes * 60, 3)} s\n"
          f"Avg time/sim = {round((timeEnd - timeStart)/len(commands), 3)}")

def getParamsFromDir(dir: str) -> tuple[float, float, float]:
    baseDirName = os.path.basename(dir)
    split = baseDirName.split('_')
    B = float(split[1])
    Vb = float(split[3])
    Vt = float(split[5])

    return B, Vb, Vt

def plotIm(fig, ax, x, Vb_unique, B_unique, Vt, cbar_label):
    im = ax.imshow(x, extent=(Vb_unique[0], Vb_unique[-1], B_unique[0], B_unique[-1]),
                   aspect=(Vb_unique[-1] - Vb_unique[0]) / (B_unique[-1] - B_unique[0]),
                #    origin="lower", interpolation="bilinear", cmap="viridis")
                   origin="lower", interpolation="nearest", cmap="viridis")
    cbar = fig.colorbar(im)
    ax.set_title(f"Vt = {Vt} V")
    ax.set_xlabel("Vb [V]")
    ax.set_ylabel("B [T]")
    cbar.set_label(cbar_label)


def plotVgtVgb(Vgt: np.ndarray,
               Vgb: np.ndarray,
               Vb_unique: np.ndarray,
               B_unique: np.ndarray,
               Vt: float) -> None:
    fig, ax = plt.subplots(2, 1,figsize=(6, 8))

    plotIm(fig, ax[0], Vgt, Vb_unique, B_unique, Vt, "Vgt [V]")
    plotIm(fig, ax[1], Vgb, Vb_unique, B_unique, Vt, "Vgb [V]")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}VgtVgb.pdf")

def plotConductance(T_2D: np.ndarray,
                    Vb_unique: np.ndarray,
                    B_unique: np.ndarray,
                    Vt: float) -> None:
    G = T2Gau(T_2D)
    fig, ax = plt.subplots(2, 1,figsize=(10, 8), height_ratios=[2.5,1])

    plotIm(fig, ax[0], G, Vb_unique, B_unique, Vt, r"G [$\frac{e^2}{h}$]")
    # im = ax[0].imshow(G, extent=(Vb_unique[0], Vb_unique[-1], B_unique[0], B_unique[-1]),
    #                   aspect='auto', origin="lower", interpolation="nearest", cmap="viridis")
    #             #    aspect='auto', origin="lower", interpolation="bilinear", cmap="viridis")
    # cbar = fig.colorbar(im)

    # ax[0].set_title(f"Vt={Vt} [V]")

    # ax[0].set_ylabel("B [T]")
    # ax[0].set_xlabel("Vb [V]")
    # cbar.set_label(r"G [$\frac{e^2}{h}$]")

    G_middle = G[G.shape[0]//2]

    # Get colors from the same normalization as the 2D plot
    norm = Normalize(vmin=G.min(), vmax=G.max())
    cmap = plt.get_cmap("viridis")

    # Create line segments for LineCollection
    points = np.array([Vb_unique, G_middle]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Plot with LineCollection
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(G_middle)
    lc.set_linewidth(2)
    ax[1].add_collection(lc)
    ax[1].autoscale()

    ax[1].set_title(f"B = {round(B_unique[G.shape[0]//2], 2)} [T], Vt={Vt} [V]")
    ax[1].set_xlabel("Vb [V]")
    ax[1].grid()
    ax[1].set_ylabel(r"G [$\frac{e^2}{h}$]")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}Conductance.pdf")

def T2Gau(T: np.ndarray) -> np.ndarray:
    return (2*e*e/h)*T

def T2GSi(T: np.ndarray) -> np.ndarray:
    return (2*eSi*eSi/hSi)*T

def plotdGdV(T_2D: np.ndarray,
             Vb_unique: np.ndarray,
             B_unique: np.ndarray,
             Vt: float) -> None:
    fig, ax = plt.subplots(figsize=(10, 8))

    dV = np.diff(Vb_unique)
    BAxis = 0
    VAxis = 1
    dG = np.diff(T2Gau(T_2D), 1, axis = VAxis)
    dGdV = dG/dV

    im = ax.imshow(dGdV, extent=(Vb_unique[0], Vb_unique[-1], B_unique[0], B_unique[-1]),
                   aspect='auto', origin="lower", interpolation="nearest", cmap='nipy_spectral')
                    #  aspect='auto', origin="lower", interpolation="bilinear", cmap='nipy_spectral')

    cbar = fig.colorbar(im)
    ax.set_title(f"Vt={Vt}")

    ax.set_ylabel("B [G]")
    ax.set_xlabel("Vb [V]")
    cbar.set_label(r"$\frac{dG}{dV}$")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}dGdV.pdf")

def saveProcessed(T_2D: np.ndarray,
                  Vgt_2D: np.ndarray,
                  Vgb_2D: np.ndarray,
                  Vb_unique: np.ndarray,
                  B_unique: np.ndarray):
    np.savetxt(args["allResultsDir"] + "Vb.csv", Vb_unique, delimiter=',')
    np.savetxt(args["allResultsDir"] + "B.csv", B_unique, delimiter=',')
    np.savetxt(args["allResultsDir"] + "T.csv", T_2D, delimiter=',')
    np.savetxt(args["allResultsDir"] + "Vgt.csv", Vgt_2D, delimiter=',')
    np.savetxt(args["allResultsDir"] + "Vgb.csv", Vgb_2D, delimiter=',')

def processFiles(plotForVt: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    print("Processing Files")

    _, dirs = getFiles(args["allResultsDir"] + "dirs/")

    l = len(dirs)

    B_list = [] * l
    Vb_list = [] * l
    T_list = [] * l
    Vgt_list = [] * l
    Vgb_list = [] * l

    timeStart = time.time()

    for idx, dir in enumerate(dirs):
        B, Vb, Vt = getParamsFromDir(dir)
        if Vt != plotForVt:
            continue

        if (idx % 10 == 0) or (idx+1 == l):
            progressBar(idx+1, 1, l, timeStart)

        data = read_csv(os.path.join(dir, "single_T.dat"), ',', header = [0])
        T   = data[1]
        Vgt = data[4]
        Vgb = data[5]

        B_list.append(B)
        Vb_list.append(Vb)
        T_list.append(T)
        Vgt_list.append(Vgt)
        Vgb_list.append(Vgt)

    if len(B_list) == 0:
        print(f"Nothing to plot for Vt = {plotForVt}")
        return np.zeros(0), np.zeros(0), np.zeros(0)

    print(len(B_list))
    print(len(dirs))

    # Get unique sorted values for grid
    B_unique = np.sort(np.unique(B_list))
    Vb_unique = np.sort(np.unique(Vb_list))

    # Create 2D array (Vb rows, B columns)
    T_2D = np.zeros((len(B_unique), len(Vb_unique)))
    Vgt_2D = np.zeros((len(B_unique), len(Vb_unique)))
    Vgb_2D = np.zeros((len(B_unique), len(Vb_unique)))

    # Fill the 2D array
    for B, Vb, T, Vgt, Vgb in zip(B_list, Vb_list, T_list, Vgt_list, Vgb_list):
        i = np.where(B_unique == B)[0][0]
        j = np.where(Vb_unique == Vb)[0][0]
        T_2D[i, j] = T
        Vgt_2D[i, j] = Vgt
        Vgb_2D[i, j] = Vgb

    saveProcessed(T_2D, Vgt_2D, Vgb_2D, Vb_unique, B_unique)

    return T_2D, Vgt_2D, Vgb_2D, Vb_unique, B_unique

def readFiles(plotForVt: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    Vb = np.loadtxt(args["allResultsDir"] + "Vb.csv", delimiter=',')
    B = np.loadtxt(args["allResultsDir"] + "B.csv", delimiter=',')
    T = np.loadtxt(args["allResultsDir"] + "T.csv", delimiter=',')
    Vgt = np.loadtxt(args["allResultsDir"] + "Vgt.csv", delimiter=',')
    Vgb = np.loadtxt(args["allResultsDir"] + "Vgb.csv", delimiter=',')

    return T, Vgt, Vgb, Vb, B

def plotAll(plotForVt) -> None:

    if (not os.path.exists(args["allResultsDir"] + "T.csv") or \
        not os.path.exists(args["allResultsDir"] + "Vgt.csv") or \
        not os.path.exists(args["allResultsDir"] + "Vgb.csv") or
        args["processFiles"]):
        T_2D, Vgt, Vgb, Vb_unique, B_unique = processFiles(plotForVt)
    else:
        T_2D, Vgt, Vgb, Vb_unique, B_unique = readFiles(plotForVt)

    if len(T_2D) == 0: return

    plotConductance(T_2D, Vb_unique, B_unique, plotForVt)
    plotVgtVgb(Vgt, Vgb, Vb_unique, B_unique, plotForVt)
    plotdGdV(T_2D, Vb_unique, B_unique, plotForVt)

################################################################################

if __name__ == "__main__":

    if (len(sys.argv) == 2):
        if (sys.argv[1] == "help" or sys.argv[1] == "--help" or sys.argv[1] == "-help"):
            printUsage(args)
            exit(0)

    parseArgs(args)
    printArgs(args)

    # exit(0)

    if (args["runSim"] == 1):
        runSim()

    if (args["plotAll"] == 1):
        VtTab = createTab(args["VtMin"], args["VtMax"], args["numVt"])
        for Vt in VtTab:
            plotAll(Vt)
        if (args["show"] == 1):
            plt.show()

#TODO run simulations, create plotting script, parallelize