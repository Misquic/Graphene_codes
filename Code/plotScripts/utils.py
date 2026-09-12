import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
import time

T2au = 4.254382E-6
au2T = 1/T2au
eV2au = 0.03674932587122423
au2eV = 1/eV2au
nm2au = 1.0/0.0529
au2nm = 1/nm2au
e = 1.
h = 1.
V2au = 0.03674932587122423/e
au2V = 1/V2au
cm2au = 1e-2 * 1e9 * nm2au
inv_cmSq2au = 1. / cm2au / cm2au
au2inv_cmSq = cm2au * cm2au

m2au = 1e9 * nm2au
inv_mSq2au = 1. / m2au / m2au
au2inv_mSq = m2au * m2au

eSi = 1.602176634e-19 # [C]
hSi = 6.62607015e-34  # [Js]

def T2Gau(T: np.ndarray) -> np.ndarray:
    return (2*e*e/h)*T


def T2GSi(T: np.ndarray) -> np.ndarray:
    return (2*eSi*eSi/hSi)*T

def R2Si(R: np.ndarray) -> np.ndarray:
    return R*eSi/hSi

def read_csv(path: str, delimiter = ' ', header = None) -> np.ndarray:
    data = pd.read_csv(path, delimiter = delimiter, header = header)
    data = np.array(data)
    if(data.shape[0] == 1):
        data = data.flatten()
    return data


def getFiles(dir, ext = None) -> tuple[list[str], list[str]]:
    absPath = os.path.abspath(dir)
    files = []
    dirs = []

    fileOrFolderList = os.listdir(absPath)
    # print(f"files in folder {dir}: {fileOrFolderList}")
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


def printUsage(args: dict) -> None:
    if not hasattr(printUsage, "isUsagePrinted"):
        printUsage.isUsagePrinted = False

        print(f"Usage <arg>[default value]: {sys.argv[0]} ", end="")
        for key, val in args.items():
            print(f"<{key}>[{val}] ", end = "")
        print()
        printUsage.isUsagePrinted = True
    else:
        assert_mess(False, "usage already printed")


def parseArgs(args: dict) -> None:
    printUsage(args)

    # assert_mess(len(sys.argv[1:]) <= len(args), "Too much arguments")
    keys = [k for k in args.keys()]

    # parse positional args
    numParsed = 1
    for index, arg in enumerate(sys.argv[1:]):
        if (arg[0] == '-' and arg[1] not in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']) or (len(arg) > 1 and (arg[1] == '-')):
            break
        args[keys[index]] = float(arg)
        numParsed += 1

    # parse named args
    for index, arg in enumerate(sys.argv[numParsed:]):
        print(arg)
        if ((arg[0] == '-') or (arg[1] == '-')):
            replIndex = 2
            if (arg[1] != '-'): replIndex = 1
            if ("=" in arg):
                argName = arg.replace('-', "", replIndex)
                argValue = argName.split('=')[1]
                argName = argName.split('=')[0]
            elif (index != len(sys.argv) - 1):
                argName = arg.replace('-', "", replIndex)
                argValue = sys.argv[index + 1]
            else:
                print("Something wrong with arguments")
                exit(1)

            if argName in args.keys():
                argType = type(args[argName])
            else:
                argType = float
                print(f"{argName} not in arguments")

            argValue = argType(argValue)
            args[argName] = argValue

    if (args["allResultsDir"][len(args["allResultsDir"])-1] != '/'):
        args["allResultsDir"] += '/'


def printArgs(args: dict) -> None:
    for key, val in args.items():
        print(f"{key} = {val}")
    print()


def progressBar(current, min, max, timeStart, end = '\n') -> None:
    range = max - min + 1
    dist = current - min

    timeNow = time.time()
    timeTaken = timeNow - timeStart
    done = dist / range
    toDo = 1 - done
    timeToFinish = (timeTaken * toDo)/np.max([done, 0.000001])
    minutes = timeToFinish // 60
    print(f"Progress: {current}/{range} = {str(round(dist/range*100, 1)).rjust(4)}%, "
          f"ETA: {minutes} min {round(timeToFinish - minutes * 60, 1)} s   ", end = end, flush = True)
    if (current == max):
        print('')

def getParamsFromDir(dir: str) -> tuple[float, float, float]:
    if (dir[-1] == '/'):
        dir = dir.split('/')[-2]
    baseDirName = os.path.basename(dir)
    split = baseDirName.split('_')
    B = float(split[1])
    Vb = float(split[3])
    Vt = float(split[5])

    return B, Vb, Vt