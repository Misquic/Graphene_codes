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
eSi = 1.602176634e-19 # [C]
hSi = 6.62607015e-34  # [Js]

def read_csv(path: str, delimiter = ' ', header = None) -> np.ndarray:
    data = pd.read_csv(path, delimiter = delimiter, header = header)
    data = np.array(data)
    if(data.shape[0] == 1):
        data = data.flatten()
    return data

def getFiles(dir) -> tuple[list[str], list[str]]:
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

    return files, dirs

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
    for index, arg in enumerate(sys.argv[1:]):
        if (arg[0] == '-') and (arg[1] == '-'):
            break
        args[keys[index]] = float(arg)

    # parse named args
    for index, arg in enumerate(sys.argv):
        if ((arg[0] == '-') and (arg[1] == '-')):
            if ("=" in arg):
                argName = arg.replace('-', "", 2)
                argValue = argName.split('=')[1]
                argName = argName.split('=')[0]
            elif (index != len(sys.argv) - 1):
                argName = arg.replace('-', "", 2)
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

def progressBar(current, min, max, timeStart) -> None:
    range = max - min + 1
    dist = current - min

    timeNow = time.time()
    timeTaken = timeNow - timeStart
    done = dist / range
    toDo = 1 - done
    timeToFinish = (timeTaken * toDo)/np.max([done, 0.000001])
    minutes = timeToFinish // 60
    print(f"Progress: {current}/{range} = {str(round(dist/range*100, 1)).rjust(4)}%, "
          f"ETA: {minutes} min {round(timeToFinish - minutes * 60, 1)} s   ", end = "\r")
    if (current == max):
        print('')
