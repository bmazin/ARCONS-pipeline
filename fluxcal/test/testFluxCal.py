from fluxcal.fluxCal import FluxCal
from util.ObsFile import ObsFile
from util.FileName import FileName
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from matplotlib.backends.backend_pdf import PdfPages

"""
Created 2/7/2013 by Seth Meeker
Routine for testing the generation of FluxCal files.
usage: > python testFluxCal.py paramfilename
example parameter files are included in the fluxcal/test for all Pal2012 standards
outputs the fluxCalSoln h5 file and a .pdf of the debug plots. Can be set to plots=False to turn this off
"""

def main():

    params = []
    paramfile = sys.argv[1]
    f = open(paramfile,'r')
    for line in f:
        params.append(line)
    f.close()
    datadir = params[0].split('=')[1].strip()
    flatdir = params[1].split('=')[1].strip()
    wvldir = params[2].split('=')[1].strip()
    fluxfile = params[3].split('=')[1].strip()
    skyfile = params[4].split('=')[1].strip()
    flatfile = params[5].split('=')[1].strip()
    wvlfile = params[6].split('=')[1].strip()
    fluxdir = params[7].split('=')[1].strip()
    fluxoutfile = params[8].split('=')[1].strip()
    objectName = params[9].split('=')[1].strip()

    fluxFileName = os.path.join(datadir, fluxfile)
    skyFileName = os.path.join(datadir, skyfile)
    wvlCalFileName = os.path.join(wvldir, wvlfile)
    flatCalFileName = os.path.join(flatdir, flatfile)
    fluxCalFileName = os.path.join(fluxdir, fluxoutfile)

    print objectName

    fc = FluxCal(fluxFileName, skyFileName, wvlCalFileName,flatCalFileName, objectName, fluxCalFileName, plots=True)

if __name__ == '__main__':
    main()

