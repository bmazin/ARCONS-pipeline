#!/bin/python

'''
Author: Paul Szypryt		Date: October 29, 2012

Loads up an h5 observation file and uses the PyGuide package to calculate the centroid of an object,
for a given frame.  Requires an initial x and y pixel position guess for the centroid.  Uses this guess
if PyGuide algorithm fails to find the centroid.
'''


import numpy as np
from tables import *
import sys
import ephem
import PyGuide as pg
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util.ObsFile import ObsFile 

import os
from PyQt4.QtGui import *
import hotpix.hotPixels as hp


# Class to allow clicking of a pixel in plt.matshow and storing the xy position of the click in an array.
class MouseMonitor():
    def __init__(self):
        pass

    def on_click(self,event):
	self.xyguess = [event.xdata,event.ydata]
   
    def connect(self):
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

# Choose obs, wavecal, and flatcal files
app = QApplication(sys.argv)

startingObsDirectory = '/ScienceData/PAL2012/'
#obsFn =str(QFileDialog.getOpenFileName(None, 'Choose Observation File',startingObsDirectory, filter=str("H5 (*.h5)")))
obsFn = '/ScienceData/PAL2012/20121208/obs_20121209-120530.h5'

index1 = obsFn.find('_')
index2 = obsFn.find('-')
utcDate = int(obsFn[index1+1:index2])
sunsetDate = utcDate-1
startingWfnDirectory = '/Scratch/waveCalSolnFiles/' + str(sunsetDate)
#wfn = str(QFileDialog.getOpenFileName(None, 'Choose Wavelength Calibration File',startingWfnDirectory, filter=str("H5 (*.h5)")))
wfn = '/Scratch/waveCalSolnFiles/20121208/calsol_20121209-131132.h5'

startingFfnDirectory = '/Scratch/flatCalSolnFiles/'
#ffn = str(QFileDialog.getOpenFileName(None, 'Choose Flat Calibration File',startingFfnDirectory, filter=str("H5 (*.h5)")))
ffn = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

# Create ObsFile instance
ob = ObsFile(obsFn)

# Load wavelength and flat cal solutions
ob.loadWvlCalFile(wfn)
ob.loadFlatCalFile(ffn)
ob.setWvlCutoffs(3000,5000)

# Load/generate hot pixel mask file
hotPixFn = '/Scratch/timeMasks/timeMask' + obsFn[index1:]
if not os.path.exists(hotPixFn):
    hp.findHotPixels(obsFn,hotPixFn)
    print "Flux file pixel mask saved to %s"%(hotPixFn)
ob.loadHotPixCalFile(hotPixFn,switchOnMask=True)
print "Hot pixel mask loaded %s"%(hotPixFn)

# Get exptime from header.  Also choose guess time.  This will be the time over which a guess will be valid
# for the centroid position.  Actual centroid calculated using this guess.  For example, choosing a guess
# time of 300 will use the same guess for the entire observation file.
exptime = ob.getFromHeader('exptime')
guessTime = 300

# Pick an integration time over which to centroid.
integrationTime=10

# Get array size information from obs file
gridHeight = ob.nCol
gridWidth = ob.nRow

# Create saturated pixel mask.  Leave as zero since hot pixel masking already does this for us.
saturatedMask = np.zeros((gridWidth,gridHeight))

# Generate dead pixel mask, invert obsFile deadMask format to put it into PyGuide format
deadMask = ob.getDeadPixels()
deadMask = -1*deadMask + 1

# Specify CCDInfo (bias,readNoise,ccdGain,satLevel)
ccd = pg.CCDInfo(0,0.00001,1,2500)

# Create output file to save centroid data
outFn = '/home/pszypryt/Scratch/centroid_test/centroid_list.txt'
f = open(outFn,'w')

for iFrame in range(exptime):
    # Integrate over the guess time.  Click a pixel in the plot corresponding to the xy center guess.  This will be used to centroid for the duration of guessTime.
    if iFrame%guessTime == 0:
	# Use obsFile to get guessTime image.
        imageInformation = ob.getPixelCountImage(firstSec=iFrame, integrationTime= guessTime, weighted=True,fluxWeighted=False, getRawCount=False,scaleByEffInt=False)
        image=imageInformation['image']
        map = MouseMonitor()
        map.fig = plt.figure()
        map.ax = map.fig.add_subplot(111)
        map.ax.set_title('Object 1')
        map.ax.matshow(image,cmap = plt.cm.gray, origin = 'lower')
        map.connect()
        plt.show()
	try:
    	    xyguess = map.xyguess
        except AttributeError:
	    pass
        print 'Guess = ' + str(xyguess)
    # Centroid an image that has been integrated over integrationTime.
    if iFrame%integrationTime == 0:
	# Use obsFile to get integrationTime image.
        imageInformation = ob.getPixelCountImage(firstSec=iFrame, integrationTime= integrationTime, weighted=True,fluxWeighted=False, getRawCount=False,scaleByEffInt=False)
        image=imageInformation['image']        
	# Use PyGuide centroiding algorithm.
        pyguide_output = pg.centroid(image,deadMask,saturatedMask,xyguess,3,ccd,0,False)
	# Use PyGuide centroid positions, if algorithm failed, use xy guess center positions instead
        try:
            xycenter = [float(pyguide_output.xyCtr[0]),float(pyguide_output.xyCtr[1])]
            print 'Calculated [x,y] center = ' + str((xycenter)) + ' for frame ' + str(iFrame) +'.'
        except TypeError:
            print 'Cannot centroid frame' + str(iFrame) + ', using guess instead'
            xycenter = xyguess
        # Write data to file
        f=open(outFn,'a')
        f.write(str(iFrame) + '\t' + str(xycenter[0]) + '\t' + str(xycenter[1]) + '\n')
        f.close()

