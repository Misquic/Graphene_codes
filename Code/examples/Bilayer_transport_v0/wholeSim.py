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
    "BMin"         : 1.,           "BMax" : 8.,                # min and max induction in T
    "VbMin"        : 1.,           "VbMax": 10.,               # min and max bottom gate voltages in V
    "VtMin"        : 0.,           "VtMax": 0.,                # min and max top gate voltages in V
    "maxParallel"  : 2,                                        # maximum number of parralel simulations
    "numB"         : 5,            "numVb": 5,   "numVt": 1,   # number of B/Vt/Vb values to run
    "dB"           : -1.,          "dVb"  : -1.,               # delta used instead of num, default is negative so it is not used
    "clearDir"     : 0,                                        # if 1 then clears dir if it exists
    "saveOld"      : 0,                                        # if 1 then saves previous folder under the same name into old
    "plotAll"      : 0,                                        # if 1 then plots results after simulation
    "runSim"       : 0,                                        # if 1 then runs simulation
    "allResultsDir": "./results/",                             # directory to store all results
    "saveStdout"   : 0,                                        # if 1 then saves outputs from simulations to file else >dev/null
    "saveSystem"   : 0,                                        # if 1 then saves created system to file
    "processFiles" : 0,                                        # if 1 then forces processing files even when previously processed
    "show"         : 0,                                        # if 1 then shows plots at the end of plotting
    "runTransport" : 1,                                        # if 1 then runTransport = 1
    "sf"           : 4,                                        # scaling factor
    "filter"       : 0,                                        # if 1 then finters results before plot
    "prepCmdsOnly" : 0,                                        # if 1 then only prepares commands and doesn't run sims
    "Executable"   : "./Transport2D",                          # exetucable to simulation
    "cmap"         : "viridis",                                # cmap for plots
    "cut"          : 0,                                        # cut plots to Bmin Bmax Vbmin Vbmax
}

def createTab(min, max, num = 3) -> np.ndarray:
    if min == max:
        return np.array([min])
    else:
        return np.linspace(min, max, num)

def execCommand(command: str) -> int:
    result = os.system(command)
    if result != 0:
        print(f"Command {command} executed with result={result}\n")
    return result

def prepareCommandsAndDirs()-> list[str]:
    saveSystem    = args["saveSystem"]
    runTransport  = 1 if args["runTransport"] == 1 else 0
    runEnergyScan = 0
    plotResults   = 0
    saveDensities = 0
    saveBands     = 0
    sf            = args["sf"]

    if (args["dB"] <= 0):
        BTab = createTab(args["BMin"], args["BMax"], args["numB"])
    else:
        BTab = createTab(args["BMin"], args["BMax"],
                         int((args["BMax"] - args["BMin"]) / args["dB"] + 1))

    if (args["dVb"] <= 0):
        VbTab = createTab(args["VbMin"], args["VbMax"], args["numVb"])
    else:
        VbTab = createTab(args["VbMin"], args["VbMax"],
                          int((args["VbMax"] - args["VbMin"]) / args["dVb"] + 1))

    VtTab = createTab(args["VtMin"], args["VtMax"], args["numVt"])
    print(f"B = {BTab}")
    print(f"Vb = {VbTab}")
    print(f"Vt = {VtTab}")

    commands = []
    allResultsDir = args["allResultsDir"]

    if (args["clearDir"] == 1):
        temp = allResultsDir[:-1] + "Old/"

        if (args["saveOld"] == 1):
            # rm temp if exists
            if os.path.exists(temp):
                execCommand(f"rm -r {temp}")
            # make temp
            execCommand(f"mkdir -p {temp}")

        if os.path.exists(allResultsDir):
            # move to previous to temp
            if (args["saveOld"] == 1):
                execCommand(f"mv {allResultsDir}* {temp}")
            # rm previous
            execCommand(f"rm -r {allResultsDir}")

        if (args["saveOld"] == 1):
            # remove previous old
            execCommand(f"rm -r ./results/Old/{temp.split('/')[-2]}")
            # mv temp to old
            execCommand(f"mv {temp} ./results/Old/")

    # create new results dir and dirs inside
    execCommand(f"mkdir --p {allResultsDir}")
    execCommand(f"mkdir --p {allResultsDir}dirs/")

    # # command to save system
    # command = f"{args["Executable"]} {allResultsDir} {np.max(BTab)} -60 0 1 0 0 0 0 0 {sf}"
    # commands.append(command)

    # Vt is outside because it is likely to be single value
    maxIdx = len(VtTab) * len(BTab) * len(VbTab)
    idx = 0
    timeStart = time.time()
    for Vt in VtTab:
        for B in BTab:
            for Vb in VbTab:
                idx += 1
                if (idx % 100 == 0):
                    progressBar(idx, 0, maxIdx, timeStart, end = '\r')
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

                command = f"{args["Executable"]} "
                for a in commandArgs:
                    command += str(a) + " "
                # print(f"running with args: {commandArgs}")
                if (args["saveStdout"] == 1):
                    command += f"1>{resultsDir}stdout.txt "
                else:
                    command += f"1>/dev/null "
                command += f"2>{resultsDir}stderr.ansi "
                commands.append(command)

    progressBar(maxIdx, 0, maxIdx, timeStart)

    with open(allResultsDir + "commands.txt", 'w') as f:
        for idx, command in enumerate(commands):
            print(f"{command}", file=f)

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

    time.sleep(5)
    if (args["prepCmdsOnly"] == 1):
        return
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

def plotIm(fig, ax, x, Vb, B, Vt, cbar_label):
    im = ax.imshow(x, extent=(Vb[0], Vb[-1], B[0], B[-1]),
                   aspect=(Vb[-1] - Vb[0]) / (B[-1] - B[0]),
                #    origin="lower", interpolation="bilinear", cmap="viridis")
                   origin="lower", interpolation="nearest", cmap=args["cmap"])

    ax.set_title(f"Vt = {Vt} V")
    ax.set_xlabel("Vb [V]")
    ax.set_ylabel("B [T]")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label)

def differenciate(image, x, y = None):
    if y is None:
        xAxis = 1
        dX = np.diff(x)
        dImage = np.diff(image, 1, axis = xAxis)
        dImagedX = dImage/dX
        return dImagedX
    else:
        yAxis = 0
        dY = np.diff(y)
        dImage = np.diff(image, 1, axis = yAxis)
        dImagedY = (dImage.transpose()/dY).transpose()
        return dImagedY

def plotVgtVgb(Vgt: np.ndarray,
               Vgb: np.ndarray,
               Vb: np.ndarray,
               B: np.ndarray,
               Vt: float) -> None:
    print("Plotting Voltages")

    fig, ax = plt.subplots(2, 2, figsize=(18, 16))

    Vgt = Vgt * au2V
    Vgb = Vgb * au2V

    plotIm(fig, ax[0,0], Vgt, Vb, B, Vt, "Vgt [V]")
    dVgtdVb = differenciate(Vgt, Vb)
    plotIm(fig, ax[0,1], dVgtdVb, Vb, B, Vt, r"$\frac{dVgt}{dVb}$ [V/V]")

    plotIm(fig, ax[1,0], Vgb, Vb, B, Vt, "Vgb [V]")
    dVgbdVb = differenciate(Vgb, Vb)
    plotIm(fig, ax[1,1], dVgbdVb, Vb, B, Vt, r"$\frac{dVgb}{dVb}$ [V/V]")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}VgtVgb.pdf")

def plotE0tE0b(E0t: np.ndarray,
               E0b: np.ndarray,
               Vb: np.ndarray,
               B: np.ndarray,
               Vt: float) -> None:
    print("Plotting Energies")

    fig, ax = plt.subplots(2, 2, figsize=(18, 16))

    E0t = filter(filter(E0t)) * au2eV
    E0b = filter(filter(E0b)) * au2eV

    plotIm(fig, ax[0,0], E0t, Vb, B, Vt, "E0t [eV]")
    dE0tdB = differenciate(E0t, Vb, B)
    plotIm(fig, ax[0,1], dE0tdB, Vb, B, Vt, r"$\frac{dE0t}{dB}$ [V/T]")

    plotIm(fig, ax[1,0], E0b, Vb, B, Vt, "E0b [V]")
    dE0bdB = differenciate(E0b, Vb, B)
    plotIm(fig, ax[1,1], dE0bdB, Vb, B, Vt, r"$\frac{dE0b}{dB}$ [V/T]")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}E0tE0b.pdf")

def plotDensities(nt: np.ndarray,
             nb: np.ndarray,
             Vb: np.ndarray,
             B: np.ndarray,
             Vt: float) -> None:
    print("Plotting Densities")

    fig, ax = plt.subplots(2, 2, figsize=(18, 16))



    # nt = nt * au2inv_cmSq
    nt = nt * au2inv_mSq
    # nb = nb * au2inv_cmSq
    nb = nb * au2inv_mSq

    maxN = np.max([np.max(nt), np.max(nb)])
    exp = int(np.log10(maxN)) - 1

    nt = nt / (10 ** exp)
    nb = nb / (10 ** exp)

    print(maxN, exp)

    plotIm(fig, ax[0,0], nt, Vb, B, Vt, r"nt [$10^{" + str(exp) + r"}$ $\frac{1}{\text{m}^2}$]")
    # plotIm(fig, ax[0,0], nt, Vb, B, Vt, r"nt [$10^{" + str(exp) + r"}$ $\frac{1}{\text{cm}^2}$]")
    dnt_dVb = differenciate(nt, Vb)
    plotIm(fig, ax[0,1], dnt_dVb, Vb, B, Vt, r"$\frac{dnt}{dVb}$")

    plotIm(fig, ax[1,0], nb, Vb, B, Vt, r"nt [$10^{" + str(exp) + r"}$ $\frac{1}{\text{m}^2}$]")
    # plotIm(fig, ax[1,0], nb, Vb, B, Vt, r"nt [$10^{" + str(exp) + r"}$ $\frac{1}{\text{cm}^2}$]")
    dnb_dVb = differenciate(nb, Vb)
    plotIm(fig, ax[1,1], dnb_dVb, Vb, B, Vt, r"$\frac{dnb}{dVb}$")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}densities.pdf")



def plotCrossSection(ax, image, Vb, B, Vt, y_label, frac = 0.5):
    Ycoord = int(image.shape[0] * frac)
    image_middle = image[Ycoord]

    # Get colors from the same normalization as the 2D plot
    norm = Normalize(vmin=image.min(), vmax=image.max())
    cmap = plt.get_cmap(args["cmap"])

    # Create line segments for LineCollection
    points = np.array([Vb, image_middle]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Plot with LineCollection
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(image_middle)
    lc.set_linewidth(2)
    ax.add_collection(lc)
    ax.autoscale()

    ax.set_title(f"B = {round(B[Ycoord], 2)} [T], Vt={Vt} [V]")
    ax.set_xlabel("Vb [V]")
    ax.grid()
    ax.set_ylabel(y_label)

def filter(array: np.ndarray) -> np.ndarray:
    # return array
    if (args["filter"] == 1):
        avg = np.mean(array)
        std = np.std(array)
        print (f"avg: {avg}, std: {std}")

        array[array > (avg + 2.5 * std)] = avg + std
        array[array < (avg - 2.5 * std)] = avg - std
    return array

def cutT(T_2D: np.ndarray,
         Vb:   np.ndarray,
         B:    np.ndarray):
    # percents
    min_Vb   = np.min(Vb)
    max_Vb   = np.max(Vb)
    range_Vb = max_Vb - min_Vb

    min_B   = np.min(B)
    max_B   = np.max(B)
    range_B = max_B - min_B

    assert_mess(min_Vb <= args["VbMin"], f"{min_Vb} !<= { args["VbMin"]}")
    assert_mess(max_Vb >= args["VbMax"], f"{max_Vb} !>= { args["VbMax"]}")
    assert_mess(min_B  <= args["BMin"] , f"{min_B } !<= {args["BMin"]  }")
    assert_mess(max_B  >= args["BMax"] , f"{max_B } !>= {args["BMax"]  }")

    lp_Vb = (args["VbMin"] - min_Vb) / range_Vb
    hp_Vb = (args["VbMax"] - min_Vb) / range_Vb
    lp_B  = (args["BMin"]  - min_B)  / range_B
    hp_B  = (args["BMax"]  - min_B)  / range_B

    # indexes
    li_Vb = int(T_2D.shape[1] * lp_Vb)
    hi_Vb = int(T_2D.shape[1] * hp_Vb)
    li_B  = int(T_2D.shape[0] * lp_B )
    hi_B  = int(T_2D.shape[0] * hp_B )

    T_2D = T_2D[li_B:hi_B, li_Vb:hi_Vb]
    Vb = Vb[li_Vb:hi_Vb]
    B = B[li_B:hi_B]

    return T_2D, Vb, B

def plotConductance(T_2D: np.ndarray,
                    Vb:   np.ndarray,
                    B:    np.ndarray,
                    Vt:   float) -> None:
    print("Plotting Conductance")

    G = T2Gau(T_2D)
    fig, ax = plt.subplots(2, 3, figsize=(28, 11), height_ratios=[3.5,1])

    plotIm(fig, ax[0,0], G, Vb, B, Vt, r"$G$ [$\frac{e^2}{h}$]")
    plotCrossSection(ax[1,0], G, Vb, B, Vt, r"$G$ [$\frac{e^2}{h}$]")

    dGdVb = differenciate(G, Vb)
    plotIm(fig, ax[0,1], dGdVb, Vb, B, Vt, r"$\frac{dG}{dVb}$")
    plotCrossSection(ax[1,1], dGdVb, Vb[:-1], B, Vt, r"$\frac{dG}{dVb}$")

    dGdB = differenciate(G, Vb, B)
    plotIm(fig, ax[0,2], dGdB, Vb, B, Vt, r"$\frac{dG}{dB}$")
    plotCrossSection(ax[1,2], dGdB, Vb, B[:-1], Vt, r"$\frac{dG}{dB}$")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}Conductance.pdf")

def plotResistance(T_2D: np.ndarray,
                   Vb:   np.ndarray,
                   B:    np.ndarray,
                   Vt:   float) -> None:
    print("Plotting Resistance")

    R = T2RSi(T_2D) / 1000 # get in kiloOhms
    log10R = np.log10(R)
    fig, ax = plt.subplots(2, 3, figsize=(28, 11), height_ratios=[3.5,1])

    plotIm(fig, ax[0,0], log10R, Vb, B, Vt, r"$log_{10}(R)$ [$log_{10}(k\Omega)$]")
    plotCrossSection(ax[1,0], log10R, Vb, B, Vt, r"$log_{10}(R)$ [$log_{10}(k\Omega)$]")

    dGdVb = differenciate(log10R, Vb)
    plotIm(fig, ax[0,1], dGdVb, Vb, B, Vt, r"$log_{10}(\frac{dR}{dVb})$")
    plotCrossSection(ax[1,1], dGdVb, Vb[:-1], B, Vt, r"$log_{10}(\frac{dR}{dVb})$")

    dGdB = differenciate(log10R, Vb, B)
    plotIm(fig, ax[0,2], dGdB, Vb, B, Vt, r"$log_{10}(\frac{dR}{dB})$")
    plotCrossSection(ax[1,2], dGdB, Vb, B[:-1], Vt, r"$log_{10}(\frac{dR}{dB})$")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}Resistance.pdf")

def T2Gau(T: np.ndarray) -> np.ndarray:
    return (2*e*e/h)*T

def T2GSi(T: np.ndarray) -> np.ndarray:
    return (2*eSi*eSi/hSi)*T

def T2RSi(T: np.ndarray) -> np.ndarray:

    T[T < 2] = 2
    R = 1 / T2GSi(T)
    return R

def plotdGdV(T_2D: np.ndarray,
             Vb:   np.ndarray,
             B:    np.ndarray,
             Vt:   float) -> None:
    print("Plotting dGdV")
    fig, ax = plt.subplots(figsize=(14, 9))

    dV = np.diff(Vb)
    BAxis = 0
    VAxis = 1
    dG = np.diff(T2Gau(T_2D), 1, axis = VAxis)
    dGdV = dG/dV

    plotIm(fig, ax, dGdV, Vb, B, Vt, r"$\frac{dG}{dV}$")

    ax.set_title(f"Vt={Vt}")
    ax.set_ylabel("B [T]")
    ax.set_xlabel("Vb [V]")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}dGdV.pdf")

def plotdGdB(T_2D, Vb, B, plotForVt):
    print("Plotting dGdB")

    fig, ax = plt.subplots(figsize=(14, 9))

    dV = np.diff(Vb)
    dB = np.diff(B)
    BAxis = 0
    VAxis = 1
    dG = np.diff(T2Gau(T_2D), 1, axis = BAxis)
    dGdB = np.divide(dG.transpose(), dB).transpose()

    plotIm(fig, ax, dGdB, Vb, B, Vt, r"$\frac{dG}{dB}$")

    ax.set_title(f"Vt={Vt}")
    ax.set_ylabel("B [T]")
    ax.set_xlabel("Vb [V]")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}dGdB.pdf")


def plotOnsites(E0t: np.ndarray,
                E0b: np.ndarray,
                Vgt: np.ndarray,
                Vgb: np.ndarray,
                Vb: np.ndarray,
                B: np.ndarray,
                plotForVt: float):
    print("Plotting Onsites")

    fig, ax = plt.subplots(2, 3, figsize=(26, 16))

    E0t = E0t * au2eV
    E0b = E0b * au2eV
    Vgt = Vgt * au2V
    Vgb = Vgb * au2V

    onsiteT = - E0t - Vgt
    onsiteB = - E0b - Vgb

    plotIm(fig, ax[0,0], onsiteT, Vb, B, Vt, "-E0t - Vgt [V]")
    dOnsiteTdVb = differenciate(onsiteT, Vb)
    plotIm(fig, ax[0,1], dOnsiteTdVb, Vb, B, Vt, r"$\frac{d(-E0t - Vgt)}{dVb}$ [V]")
    dOnsiteTdB = differenciate(onsiteT, Vb, B)
    plotIm(fig, ax[0,2], dOnsiteTdB, Vb, B, Vt, r"$\frac{d(-E0t - Vgt)}{dB}$ [V]")

    plotIm(fig, ax[1,0], onsiteB, Vb, B, Vt, "-E0b - Vgb [V]")
    dOnsiteBdVb = differenciate(onsiteB, Vb)
    plotIm(fig, ax[1,1], dOnsiteBdVb, Vb, B, Vt, r"$\frac{d(-E0b - Vgb)}{dVb}$ [V/V]")
    dOnsiteBdB = differenciate(onsiteB, Vb, B)
    plotIm(fig, ax[1,2], dOnsiteBdB, Vb, B, Vt, r"$\frac{d(-E0b - Vgb)}{dB}$ [V/V]")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}onsites.pdf")

def saveProcessed(T_2D: np.ndarray,
                  Vgt_2D: np.ndarray,
                  Vgb_2D: np.ndarray,
                  E0t_2D: np.ndarray,
                  E0b_2D: np.ndarray,
                  nb_2D: np.ndarray,
                  nt_2D: np.ndarray,
                  Vb: np.ndarray,
                  B: np.ndarray):
    np.savetxt(args["allResultsDir"] + "Vb.csv", Vb, delimiter=',')
    np.savetxt(args["allResultsDir"] + "B.csv", B, delimiter=',')
    np.savetxt(args["allResultsDir"] + "T.csv", T_2D, delimiter=',')
    np.savetxt(args["allResultsDir"] + "Vgt.csv", Vgt_2D, delimiter=',')
    np.savetxt(args["allResultsDir"] + "Vgb.csv", Vgb_2D, delimiter=',')
    np.savetxt(args["allResultsDir"] + "E0t.csv", E0t_2D, delimiter=',')
    np.savetxt(args["allResultsDir"] + "E0b.csv", E0b_2D, delimiter=',')
    np.savetxt(args["allResultsDir"] + "nb.csv", nb_2D, delimiter=',')
    np.savetxt(args["allResultsDir"] + "nt.csv", nt_2D, delimiter=',')

def processFiles(plotForVt: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    print("Processing Files")

    _, dirs = getFiles(args["allResultsDir"] + "dirs/")

    l = len(dirs)

    B_list = [] * l
    Vb_list = [] * l
    T_list = [] * l
    Vgt_list = [] * l
    Vgb_list = [] * l
    E0t_list = [] * l
    E0b_list = [] * l
    nb_list = [] * l
    nt_list = [] * l

    timeStart = time.time()

    for idx, dir in enumerate(dirs):
        B, Vb, Vt = getParamsFromDir(dir)
        if Vt != plotForVt:
            continue

        if (idx % 10 == 0) or (idx+1 == l):
            progressBar(idx+1, 1, l, timeStart, "\r")

        try:
            # Ef(au),T_total,E0t(au),E0b(au),Vgt(au),Vgb(au),nt(au),nb(au)
            data = read_csv(os.path.join(dir, "single_T.dat"), ',', header = [0])
            T   = data[1]
            E0t = data[2]
            E0b = data[3]
            Vgt = data[4]
            Vgb = data[5]
            nt = data[6]
            nb = data[7]
        except:
            print(f"folder with problem: {dir}")
            T   = 0
            E0t = 0
            E0b = 0
            Vgt = 0
            Vgb = 0
            nt = 0
            nb = 0

        if T == 0:
            print(f"T = 0 in ")

        B_list.append(B)
        Vb_list.append(Vb)
        T_list.append(T)
        E0t_list.append(E0t)
        E0b_list.append(E0b)
        Vgt_list.append(Vgt)
        Vgb_list.append(Vgb)
        nt_list.append(nt)
        nb_list.append(nb)

    if len(B_list) == 0:
        print(f"Nothing to plot for Vt = {plotForVt}")
        return np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)

    print(len(B_list))
    print(len(dirs))

    # Get unique sorted values for grid
    B_unique = np.sort(np.unique(B_list))
    Vb_unique = np.sort(np.unique(Vb_list))

    # Create 2D array (Vb rows, B columns)
    T_2D   = np.zeros((len(B_unique), len(Vb_unique)))
    Vgt_2D = np.zeros((len(B_unique), len(Vb_unique)))
    Vgb_2D = np.zeros((len(B_unique), len(Vb_unique)))
    E0t_2D = np.zeros((len(B_unique), len(Vb_unique)))
    E0b_2D = np.zeros((len(B_unique), len(Vb_unique)))
    nt_2D  = np.zeros((len(B_unique), len(Vb_unique)))
    nb_2D  = np.zeros((len(B_unique), len(Vb_unique)))

    # Fill the 2D array
    for B, Vb, T, Vgt, Vgb, E0t, E0b, nt, nb in zip(B_list, Vb_list, T_list, Vgt_list, Vgb_list, E0t_list, E0b_list, nt_list, nb_list):
        i = np.where(B_unique == B)[0][0]
        j = np.where(Vb_unique == Vb)[0][0]
        T_2D[i, j]   = T
        Vgt_2D[i, j] = Vgt
        Vgb_2D[i, j] = Vgb
        E0t_2D[i, j] = E0t
        E0b_2D[i, j] = E0b
        nt_2D[i, j]  = nt
        nb_2D[i, j]  = nb

    saveProcessed(T_2D, Vgt_2D, Vgb_2D, E0t_2D, E0b_2D, nt_2D, nb_2D, Vb_unique, B_unique)

    return T_2D, Vgt_2D, Vgb_2D, E0t_2D, E0b_2D, nt_2D, nb_2D, Vb_unique, B_unique

def readFiles(plotForVt: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    print("Reading files")
    Vb = np.loadtxt(args["allResultsDir"] + "Vb.csv", delimiter=',')
    B = np.loadtxt(args["allResultsDir"] + "B.csv", delimiter=',')
    T = np.loadtxt(args["allResultsDir"] + "T.csv", delimiter=',')
    Vgt = np.loadtxt(args["allResultsDir"] + "Vgt.csv", delimiter=',')
    Vgb = np.loadtxt(args["allResultsDir"] + "Vgb.csv", delimiter=',')
    E0t = np.loadtxt(args["allResultsDir"] + "E0t.csv", delimiter=',')
    E0b = np.loadtxt(args["allResultsDir"] + "E0b.csv", delimiter=',')
    nt = np.loadtxt(args["allResultsDir"] + "nt.csv", delimiter=',')
    nb = np.loadtxt(args["allResultsDir"] + "nb.csv", delimiter=',')

    return T, Vgt, Vgb, E0t, E0b, nt, nb, Vb, B

def plotAll(plotForVt : float) -> None:
    if (not os.path.exists(args["allResultsDir"] + "T.csv") or \
        not os.path.exists(args["allResultsDir"] + "Vgt.csv") or \
        not os.path.exists(args["allResultsDir"] + "Vgb.csv") or
        args["processFiles"]):
        T_2D, Vgt, Vgb, E0t, E0b, nt, nb, Vb, B = processFiles(plotForVt)
    else:
        T_2D, Vgt, Vgb, E0t, E0b, nt, nb, Vb, B = readFiles(plotForVt)

    if (args["filter"] == 1):
        T_2D[T_2D > 125] = 125
        T_2D = filter(T_2D)
    if len(T_2D) == 0: return

    if (args["cut"] == 1):
        T_2D, _, _ = cutT(T_2D, Vb, B)
        # T_2D = smoothZeros(T_2D)
        Vgt,  _, _ = cutT(Vgt,  Vb, B)
        Vgb,  _, _ = cutT(Vgb,  Vb, B)
        E0t,  _, _ = cutT(E0t,  Vb, B)
        E0b,  _, _ = cutT(E0b,  Vb, B)
        nt,   _, _ = cutT(nt,   Vb, B)
        nb,  Vb, B = cutT(nb,  Vb, B)

    print("Plotting")
    plotConductance(T_2D, Vb, B, plotForVt)
    plotResistance( T_2D, Vb, B, plotForVt)
    plotVgtVgb(Vgt, Vgb, Vb, B, plotForVt)
    plotE0tE0b(E0t, E0b, Vb, B, plotForVt)
    plotDensities(nt, nb, Vb, B, plotForVt)
    plotdGdV(T_2D, Vb, B, plotForVt)
    plotdGdB(T_2D, Vb, B, plotForVt)
    E0t_unique = np.unique(E0t)
    E0b_unique = np.unique(E0b)
    if (len(E0t_unique) > 1) or (len(E0b_unique) > 1):
        plotOnsites(E0t, E0b, Vgt, Vgb, Vb, B, plotForVt)

################################################################################

if __name__ == "__main__":

    if (len(sys.argv) == 2):
        if (sys.argv[1] == "help" or sys.argv[1] == "--help" or sys.argv[1] == "-help"):
            printUsage(args)
            exit(0)

    parseArgs(args)
    printArgs(args)

    if (args["runSim"] == 1) or (args["prepCmdsOnly"] == 1):
        runSim()

    if (args["plotAll"] == 1):
        VtTab = createTab(args["VtMin"], args["VtMax"], args["numVt"])
        for Vt in VtTab:
            plotAll(Vt)
        if (args["show"] == 1):
            plt.show()
