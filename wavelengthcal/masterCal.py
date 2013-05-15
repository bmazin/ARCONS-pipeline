'''
Author: Danica Marsden                          Date: May 1, 2013

For each night during the Palomar observing run, come up with a best
wavecal estimate.

'''

import numpy as np
import datetime as dt
import tables
import os.path
import matplotlib.pylab as mpl
import pyfits
import pickle
from utils import *
from tables import *
from hotPixels import *
from FileName import *
from ObsFile import *


class WaveCalSoln(IsDescription):
    roach = tables.UInt16Col()                 # ROACH board number
    pixelnum = tables.UInt16Col()              # pixel number on the roach
    pixelrow = tables.UInt16Col()              # physical x location - from beam map
    pixelcol = tables.UInt16Col()              # physical y location 
    polyfit = tables.Float64Col(3)             # polynomial to convert from phase amplitude to wavelength,
                                               #    double precision
    sigma = tables.Float64Col()                # 1 sigma (Gaussian width) in eV, for blue peak
    solnrange = tables.Float32Col(2)           # start and stop wavelengths for the fit in Angstroms
    wave_flag = tables.UInt16Col()             # flag to indicate if pixel is good (0), dead (1) or failed
                                               #    during wave cal fitting (2)


n_rows = 46
n_cols = 44
non_alloc_pix = 250
run = 'PAL2012'

## At present, picking which cal files to use based on the array plots for the wavecal solutions from that night.

## For Dec. 7, 2012.
calTimestamps = ['20121208-023522', '20121208-024054', '20121208-031248', '20121208-031719', '20121208-043620', '20121208-063741', '20121208-070505', '20121208-101613', '20121208-112506', '20121208-132558']
sunsetDates = ['20121207','20121207','20121207','20121207','20121207','20121207','20121207','20121207','20121207','20121207']

## For Dec. 8, 2012.
#calTimestamps = ['20121209-133419', '20121209-021609', '20121209-042051', '20121209-043937', '20121209-060704', '20121209-062404', '20121209-111730', '20121209-115205', '20121209-131132', '20121211-020441']
#sunsetDates = ['20121208','20121208','20121208','20121208','20121208','20121208','20121208','20121208','20121208','20121208']

## For Dec. 10, 2012.
#calTimestamps = ['20121211-020441', '20121211-044853', '20121211-052230', '20121211-072632', '20121211-074031', '20121211-090613', '20121211-115429', '20121211-133056', '20121211-135844']
#sunsetDates = ['20121210','20121210','20121210','20121210','20121210','20121210','20121210','20121210','20121210']

## For Dec. 11, 2012.
#calTimestamps = ['20121212-023031', '20121212-032455', '20121212-063518', '20121212-065247', '20121212-073251', '20121212-091828', '20121212-102334', '20121212-105050', '20121212-111847']
#sunsetDates = ['20121211','20121211','20121211','20121211','20121211','20121211','20121211','20121211','20121211']




wvlCalFiles = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
    
roachTable = np.zeros([46, 44])
pixelTable = np.zeros([46, 44])
wvlCalTable = np.zeros([46, 44, 3])
wvlErrorTable = np.zeros([46, 44])
wvlFlagTable = np.zeros([46, 44])
wvlRangeTable = np.zeros([46, 44, 2])
flagplot = np.zeros([46, 44])

counter = 0
n_pix = 0
for file1 in wvlCalFiles:
        
    print file1
    
    ## To plot the beammap:
    #obsfile = '/ScienceData/PAL2012/'+file1.split('/')[3]+'/'+file1.split('/')[4].split('sol')[0]+file1.split('/')[4].split('sol')[1]
    #print obsfile
    #data = tables.openFile(obsfile, mode='r')
    #beammap = data.root.beammap.beamimage.read()
    #beamtable = np.zeros([46, 44])  # get obs files, beam maps, plot them, check...
    #for i in range(n_rows):
    #    for j in range(n_cols):
    #        pstring = beammap[i][j]
    #        pixelNode = pstring.split('t')[0]
    #        roach = int(pixelNode.split('/')[1][1:])
    #        pixel = int(pixelNode.split('/')[2][1:])
    #        if (roach==0) & (pixel==non_alloc_pix):
    #            beamtable[i][j] = 0
    #        else:
    #            beamtable[i][j] = 1

    #plotArray(beamtable, showMe=False, plotFileName='beammap'+str(counter)+'.png')
        
    wvlCalFile = tables.openFile(file1, mode='r')
    wvlCalData = wvlCalFile.root.wavecal.calsoln

        
    for calPixel in wvlCalData:
        #print calPixel['pixelrow'], calPixel['pixelcol']
        if counter == 0:
            roachTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['roach']
            pixelTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['pixelnum']
            wvlFlagTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['wave_flag']
            wvlErrorTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['sigma']
            wvlCalTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['polyfit']
            wvlRangeTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['solnrange']
        else:
            if ((calPixel['wave_flag']==0) and (wvlFlagTable[calPixel['pixelrow']][calPixel['pixelcol']] != 0)):
                print 'improving pixel ', calPixel['pixelrow'], calPixel['pixelcol']
                wvlFlagTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['wave_flag']
                wvlErrorTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['sigma']
                wvlCalTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['polyfit']
                wvlRangeTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['solnrange']
                
    counter += 1

    
## write to file...
outfile = '/Scratch/waveCalSolnFiles/'+sunsetDates[0]+'/wvlCalFiles_combined_Dec'+sunsetDates[0][6:8]+'2012.h5'
h5out = tables.openFile(outfile, mode='w')
root = h5out.root
calgroup = h5out.createGroup(root, 'wavecal', 'Table of calibration parameters for each pixel')
caltable = h5out.createTable(calgroup, 'calsoln', WaveCalSoln, title='Wavelength Cal Table')
row = caltable.row

for i in range(46):
    for j in range(44):

        row['pixelrow'] = i
        row['pixelcol'] = j
        
        row['roach'] = roachTable[i][j]
        row['pixelnum'] = pixelTable[i][j]
                
        row['wave_flag'] = wvlFlagTable[i][j]
        row['sigma'] = wvlErrorTable[i][j]
        row['solnrange'] = wvlRangeTable[i][j]
        row['polyfit'] = wvlCalTable[i][j]

        row.append()

        if wvlFlagTable[i][j] == 0:
            flagplot[i][j] = 1
            n_pix += 1

caltable.flush()
h5out.close()

plotArray(flagplot, showMe=False, plotFileName='/Scratch/waveCalSolnFiles/'+sunsetDates[0]+'/figs/flags_combined_Dec'+sunsetDates[0][6:8]+'2012.png', plotTitle=str(n_pix)+' pixels')
plotArray(wvlErrorTable, showMe=False, cbar=True, plotFileName='/Scratch/waveCalSolnFiles/'+sunsetDates[0]+'/figs/combined_Dec'+sunsetDates[0][6:8]+'2012_arrayPlot.png',
                  plotTitle='Energy Resolution at 400nm')
