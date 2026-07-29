import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import AutoMinorLocator

import numpy as np
from utils import *
from args import *

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

    # E0t = filter(filter(E0t)) * au2eV
    # E0b = filter(filter(E0b)) * au2eV

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

def plotResistance(R: np.ndarray,
                   Vb:   np.ndarray,
                   B:    np.ndarray,
                   Vt:   float) -> None:
    print("Plotting Resistance")

    # R = R / 1000 # get in kiloOhms
    # log10R = np.abs(R)
    rUnit = r"\frac{\text{h}}{\text{e}}"

    fig, ax = plt.subplots(2, 3, figsize=(28, 11), height_ratios=[3.5,1])

    # R = np.maximum(np.minimum(R, 0.1), 0)
    plotIm(fig, ax[0,0], R, Vb, B, Vt, rf"$R$ [${rUnit}$]")
    plotCrossSection(ax[1,0], R, Vb, B, Vt, rf"$R$ [${rUnit}$]")

    dGdVb = differenciate(R, Vb)
    plotIm(fig, ax[0,1], dGdVb, Vb, B, Vt, r"$\frac{dR}{dVb}$")
    plotCrossSection(ax[1,1], dGdVb, Vb[:-1], B, Vt, r"$\frac{dR}{dVb}$")

    dGdB = differenciate(R, Vb, B)
    plotIm(fig, ax[0,2], dGdB, Vb, B, Vt, r"$\frac{dR}{dB}$")
    plotCrossSection(ax[1,2], dGdB, Vb, B[:-1], Vt, r"$\frac{dR}{dB}$")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}Resistance.pdf")


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

    plotIm(fig, ax, dGdB, Vb, B, plotForVt, r"$\frac{dG}{dB}$")

    ax.set_title(f"Vt={plotForVt}")
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

    plotIm(fig, ax[0,0], onsiteT, Vb, B, plotForVt, "-E0t - Vgt [V]")
    dOnsiteTdVb = differenciate(onsiteT, Vb)
    plotIm(fig, ax[0,1], dOnsiteTdVb, Vb, B, plotForVt, r"$\frac{d(-E0t - Vgt)}{dVb}$ [V]")
    dOnsiteTdB = differenciate(onsiteT, Vb, B)
    plotIm(fig, ax[0,2], dOnsiteTdB, Vb, B, plotForVt, r"$\frac{d(-E0t - Vgt)}{dB}$ [V]")

    plotIm(fig, ax[1,0], onsiteB, Vb, B, plotForVt, "-E0b - Vgb [V]")
    dOnsiteBdVb = differenciate(onsiteB, Vb)
    plotIm(fig, ax[1,1], dOnsiteBdVb, Vb, B, plotForVt, r"$\frac{d(-E0b - Vgb)}{dVb}$ [V/V]")
    dOnsiteBdB = differenciate(onsiteB, Vb, B)
    plotIm(fig, ax[1,2], dOnsiteBdB, Vb, B, plotForVt, r"$\frac{d(-E0b - Vgb)}{dB}$ [V/V]")

    fig.tight_layout()
    fig.savefig(f"{args["allResultsDir"]}onsites.pdf")


def plotIm(fig, ax, x, Vb, B, Vt, cbar_label):
    im = ax.imshow(x, extent=(Vb[0], Vb[-1], B[0], B[-1]),
                   aspect=(Vb[-1] - Vb[0]) / (B[-1] - B[0]),
                #    origin="lower", interpolation="bilinear", cmap="viridis")
                   origin="lower", interpolation="nearest", cmap=args["cmap"])

    ax.set_title(f"Vt = {Vt} V")
    ax.set_xlabel("Vb [V]")
    ax.set_ylabel("B [T]")

    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    if args["grid"]:
        ax.grid()
        # ax.grid(which="minor", color="white", linewidth=0.3, alpha=0.5)
        ax.grid(which="minor")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label)
