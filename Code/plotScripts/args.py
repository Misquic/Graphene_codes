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
    "cmap"         : "inferno_r",                                # cmap for plots
    "cut"          : 0,                                        # cut plots to Bmin Bmax Vbmin Vbmax
    "saveCurrents" : 0,                                        # saves and plots currents
    "leadInfo"     : "\"currentLeadFrom currentLeadTo voltageLeadHigh voltageLeadLow\"", \
                                                               # leadInfo for resistances postprocessing
    # "resistance"   : 0,
    "grid"         : 0,                                        # if 1 then plots imShow with grid
    "plotCurrents" : 0,                                        # if 1 then plots currents for points described by division with numB and numVb
    "seed"         : 12345,                                    # seed for simulations
    "Rmax"         : 0.4,                                      # maximum R when filtering is applied
    "averageFiles" : 0,                                        # if 1 then averages files first from subdirectories
    "saveBands"    : 0,
    "plotResults"  : 0,
}
