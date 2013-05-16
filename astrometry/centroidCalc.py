#!/bin/python

'''
Author: Paul Szypryt		Date: October 29, 2012

Loads up an h5 observation file and uses the PyGuide package to calculate the centroid of an object,
for a given frame.  Requires an initial x and y pixel position guess for the centroid.  Uses this guess
if PyGuide algorithm fails to find the centroid.

History:
March 25, 2013 - Now calculates the hour angle using the right ascension of the target object as a 1st
    order approximation.  To get more exact hour angle, would need right ascension of the center of 
    rotation position, but this is unknown and non-constant.
'''


import numpy as np
import tables
import sys
import ephem
import PyGuide as pg
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util.ObsFile import ObsFile 
import os
from PyQt4.QtGui import *
import hotpix.hotPixels as hp
from tables import *
from util.FileName import FileName


# Class to allow clicking of a pixel in plt.matshow and storing the xy position of the click in an array.
class MouseMonitor():
    def __init__(self):
        pass

    def on_click(self,event):
        self.xyguess = [event.xdata,event.ydata]
   
    def connect(self):
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

# Some useful conversion functions
def arcsec_to_radians(total_arcsec):
    total_degrees = total_arcsec/3600.0
    total_radians = total_degrees*d2r
    return total_radians

def radians_to_arcsec(total_radians):
    total_degrees = total_radians*r2d
    total_arcsec = total_degrees*3600.0
    return total_arcsec

def saveTable(centroidListFileName,timeList,xPositionList,yPositionList,hourAngleList,flagList):

    if os.path.isabs(centroidListFileName) == True:
        fullCentroidListFileName = centroidListFileName
    else:
        scratchDir = os.getenv('INTERM_PATH')
        centroidDir = os.path.join(scratchDir,'centroidListFiles')
        fullCentroidListFileName = os.path.join(centroidDir,centroidListFileName)
    
    try:
        centroidListFile = tables.openFile(fullCentroidListFileName,mode='w')
    except:
        print 'Error: Couldn\'t create centroid list file, ',fullCentroidListFileName
        return
    print 'wrote to', centroidListFileName

    centroidgroup = centroidListFile.createGroup(centroidListFile.root,'centroidlist','Table of times, x positions, y positions, hour angles, and flags')
    caltable = tables.Array(centroidgroup,'times',object=timeList,title='Times at which centroids were calculated')
    caltable = tables.Array(centroidgroup,'xPositions',object=xPositionList,title='X centroid positions')
    caltable = tables.Array(centroidgroup,'yPositions',object=yPositionList,title='Y centroid positions')
    caltable = tables.Array(centroidgroup,'hourAngles',object=hourAngleList,title='Hour angles at specified times')
    caltable = tables.Array(centroidgroup,'flags',object=flagList,title='Flags whether or not guess had to be used, 1 for guess used')
    centroidListFile.flush()
    centroidListFile.close()

#seq0 = ['120530', '121033','121536', '122039', '122542', '123045', '123548', '124051', '124554', '125057', '125601', '130103', '130606']

# Obs file info
run = 'PAL2012'
sunsetDate='20121208'
utcDate='20121209'

centroidTimestamp = '20121209-120530'
calTimestamp = '20121209-131132'

# Specify input parameters.
centroid_RA = '09:26:38.7'
centroid_DEC = '36:24:02.4'
HA_offset = 16.0
guessTime = 300
integrationTime=30
secondMaxCountsForDisplay = 500

# Some useful conversions
d2r = np.pi/180.0
r2d = 180.0/np.pi

# Choose obs, wavecal, and flatcal files
app = QApplication(sys.argv)

obsFn = FileName(run=run,date=sunsetDate,tstamp=centroidTimestamp).obs()
wfn = FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln()
ffn = FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln()

ffn = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

# Create ObsFile instance
ob = ObsFile(obsFn)

# Load wavelength and flat cal solutions
ob.loadWvlCalFile(wfn)
ob.loadFlatCalFile(ffn)
ob.setWvlCutoffs(3000,8000)

# Load/generate hot pixel mask file
index1 = obsFn.find('_')
index2 = obsFn.find('-')
hotPixFn = '/Scratch/timeMasks/timeMask' + obsFn[index1:]
if not os.path.exists(hotPixFn):
    hp.findHotPixels(obsFn,hotPixFn)
    print "Flux file pixel mask saved to %s"%(hotPixFn)
ob.loadHotPixCalFile(hotPixFn,switchOnMask=False)
print "Hot pixel mask loaded %s"%(hotPixFn)

# Get exptime and LST from header.  
exptime = ob.getFromHeader('exptime')
original_lst = ob.getFromHeader('lst')
print 'Original LST from telescope:', original_lst

# Initial RA and LST
centroid_RA_radians = ephem.hours(centroid_RA).real
centroid_RA_arcsec = radians_to_arcsec(centroid_RA_radians)

centroid_DEC_radians = ephem.degrees(centroid_DEC).real
centroid_DEC_arcsec = radians_to_arcsec(centroid_DEC_radians)

original_lst_radians = ephem.hours(original_lst).real
original_lst_seconds = radians_to_arcsec(original_lst_radians)/15.0

# Get array size information from obs file
gridHeight = ob.nCol
gridWidth = ob.nRow

# Create saturated pixel mask to apply to PyGuide algorithm.
print 'Creating saturation mask...'
nFrames = int(exptime/integrationTime)
saturatedMask = np.zeros(((nFrames,gridWidth,gridHeight)))
hotPixInfo = hp.readHotPixels(hotPixFn)
intervalsMatrix = hotPixInfo['intervals']
for t in range(nFrames):
    for x in range(gridHeight):
        for y in range(gridWidth):
            if intervalsMatrix[y][x] == []:
                pass
            else:
                saturatedMask[t][y][x]=1

# Generate dead pixel mask, invert obsFile deadMask format to put it into PyGuide format
print 'Creating dead mask...'
deadMask = ob.getDeadPixels()
deadMask = -1*deadMask + 1

# Specify CCDInfo (bias,readNoise,ccdGain,satLevel)
ccd = pg.CCDInfo(0,0.00001,1,2500)

# Create output file to save centroid data
centroidListFileName=FileName(run=run,date=sunsetDate,tstamp=centroidTimestamp).centroidList()
print centroidListFileName

norm = mpl.colors.Normalize(vmin=0,vmax=secondMaxCountsForDisplay*guessTime)

timeList=[]
xPositionList=[]
yPositionList=[]
hourAngleList=[]
flagList=[]

flag=0
print 'Retrieving images...'
for iFrame in range(exptime):
    # Integrate over the guess time.  Click a pixel in the plot corresponding to the xy center guess.  This will be used to centroid for the duration of guessTime.
    if iFrame%guessTime == 0:
    # Use obsFile to get guessTime image.
        imageInformation = ob.getPixelCountImage(firstSec=iFrame, integrationTime= guessTime, weighted=True,fluxWeighted=False, getRawCount=False,scaleByEffInt=False)
        image=imageInformation['image']
        map = MouseMonitor()
        map.fig = plt.figure()
        map.ax = map.fig.add_subplot(111)
        map.ax.set_title('Centroid Guess')
        map.ax.matshow(image,cmap = plt.cm.gray, origin = 'lower',norm=norm)
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
        satMask=saturatedMask[int(iFrame/integrationTime)]
        imageInformation = ob.getPixelCountImage(firstSec=iFrame, integrationTime= integrationTime, weighted=True,fluxWeighted=False, getRawCount=False,scaleByEffInt=False)
        image=imageInformation['image']        
        # Use PyGuide centroiding algorithm.
        pyguide_output = pg.centroid(image,deadMask,satMask,xyguess,3,ccd,0,False)
        # Use PyGuide centroid positions, if algorithm failed, use xy guess center positions instead
        try:
            xycenter = [float(pyguide_output.xyCtr[0]),float(pyguide_output.xyCtr[1])]
            print 'Calculated [x,y] center = ' + str((xycenter)) + ' for frame ' + str(iFrame) +'.'
            flag = 0
        except TypeError:
            print 'Cannot centroid frame' + str(iFrame) + ', using guess instead'
            xycenter = xyguess
            flag = 1
        # Begin RA/DEC mapping
        current_lst_seconds = original_lst_seconds + iFrame
        current_lst_radians = arcsec_to_radians(current_lst_seconds*15.0)
        HA_variable = current_lst_radians - centroid_RA_radians
        HA_static = HA_offset*d2r
        HA_current = HA_variable + HA_static
        HA_degrees = HA_current * r2d
        # Make lists to save to h5 file
        timeList.append(iFrame)
        xPositionList.append(xycenter[0])
        yPositionList.append(xycenter[1])
        hourAngleList.append(HA_degrees)
        flagList.append(flag)
# Save to h5 table
saveTable(centroidListFileName=centroidListFileName,timeList=timeList,xPositionList=xPositionList,yPositionList=yPositionList,hourAngleList=hourAngleList,flagList=flagList)

