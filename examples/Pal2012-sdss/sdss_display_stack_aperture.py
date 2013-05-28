#!/bin/python

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util.readDict import readDict
from util import utils
import tables
import matplotlib.pyplot as plt
import hotpix.hotPixels as hp
import os
from time import time

def circleAperture(startpx,startpy,radius=3):
    offset = np.floor(radius+0.5)
    xBinsPerPixel = yBinsPerPixel = 50
    totalBinsPerPixel = xBinsPerPixel * yBinsPerPixel
    xPixels = yPixels = 1 + 2*np.floor(0.5 + radius)
    xBins = yBins = int(xPixels*xBinsPerPixel)
    pixelsInSquare = int(xPixels * yPixels)
    squareVals = np.linspace(-xPixels/2.,xPixels/2.,xBins)
    squareVals = squareVals[:-1]

    xvals = []
    yvals = []
    for x in range(xBins-1):
        for y in range(yBins-1):
            if squareVals[x]**2 + squareVals[y]**2 <= radius**2:
                xvals.append(squareVals[x])
                yvals.append(squareVals[y])

    counts = np.zeros((int(yPixels),int(xPixels)))
    for x in range(int(xPixels)):
        for y in range(int(yPixels)):
            xPos = x - xPixels/2.0
            yPos = y - yPixels/2.0
            for i in range(len(xvals)):
                if (xvals[i] < xPos + 1.0) and (xvals[i] >= xPos) and (yvals[i] < yPos + 1.0) and (yvals[i] >= yPos):
                    counts[y][x] += 1
    pixelWeights = counts / totalBinsPerPixel
#    print pixelWeights
    pixelMask = np.zeros((46,44))
    for x in range(int(xPixels)):
        for y in range(int(yPixels)):
            if (startpy + y - offset >=0) and (startpx + x - offset >= 0):
                pixelMask[startpy + y - offset][startpx + x - offset] = pixelWeights[y][x]
    return pixelMask
#pixMask=circleAperture(0,0,3)
#print 'pixMask',pixMask
#plt.matshow(pixMask,origin = 'lower')
#plt.show()

param = readDict()
param.read_from_file('0926params.dict')

run = param['run']
flatSunsetLocalDate = param['flatSunsetLocalDate']
flatTimestamp = param['flatTimestamp']

wvlCalSunsetLocalDate = param['wvlCalSunsetLocalDate']
wvlCalTimestamp = param['wvlCalTimestamp']

expTime = param['expTime']
integrationTime = param['integrationTime']

radius1 = param['apertureRadius1']
radius2 = param['apertureRadius2']
ObjPosFile = param['ObjPosFile']

utcDates = param['utcDates']
sunsetDates = param['sunsetDates']

calTimestamps = param['calTimestamps']
seqs = param['seqs']

wvlLowerCutoff = param['wvlLowerCutoff']
wvlUpperCutoff = param['wvlUpperCutoff']
npzFileName = param['npzFileName']
gifFileName = param['gifFileName']

timestampLists = [[utcDate+'-'+str(ts) for ts in seq] for utcDate,seq in zip(utcDates,seqs)]

#wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]

wvlCalFilenames = ['/Scratch/waveCalSolnFiles/%s/calsol_%s.h5'%(sunsetDate,sunsetDate) for sunsetDate in sunsetDates]

#flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]

flatCalFilenames = ['/Scratch/flatCalSolnFiles/%s/flatsol_%s.h5'%(sunsetDate,sunsetDate) for sunsetDate in sunsetDates]

#use for Dec8
flatCalFilenames[0] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'
flatCalFilenames[1] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

frames = []
showframes = []
times = []
titles = []
LightCurve = []
count= 0

NumFiles = []
for seq in seqs:
    NumFiles.append(len(seq))
NumFiles = sum(NumFiles)
print (NumFiles)*expTime/integrationTime,'frames to make'

backgroundPerPixel=[]
for iSeq in range(len(seqs)):
    timestampList = timestampLists[iSeq]
    print timestampList
    wfn = wvlCalFilenames[iSeq]
    ffn = flatCalFilenames[iSeq]
    sunsetDate = sunsetDates[iSeq]

    for position in ObjPosFile:
        if iSeq+1 == position[0]:
            print 'finding aperture and sky masks'
            guessX = position[1]
            guessY = position[2]
            apertureMask = circleAperture(guessX,guessY,radius1)
            bigMask = circleAperture(guessX,guessY,radius2)
            skyMask = bigMask - apertureMask
#            print skyMask[0:20][5:25]
       
    for i,ts in enumerate(timestampList):
        print 'loading',ts
        obsFn = FileName(run=run,date=sunsetDate,tstamp=ts).obs()
        ob = ObsFile(obsFn)
        ob.loadTimeAdjustmentFile(FileName(run=run).timeAdjustments())
	index1 = obsFn.find('_')
	hotPixFn = '/Scratch/timeMasks/timeMask' + obsFn[index1:]
        if not os.path.exists(hotPixFn):
            hp.findHotPixels(obsFn,hotPixFn)
            print "Flux file pixel mask saved to %s"%(hotPixFn)
        ob.loadHotPixCalFile(hotPixFn,switchOnMask=True)
        ob.loadWvlCalFile(wfn)
        ob.loadFlatCalFile(ffn)
        ob.setWvlCutoffs(wvlLowerCutoff,wvlUpperCutoff)

	bad_solution_mask=np.zeros((46,44))
	bad_count=0;
	for y in range(46):
	    for x in range(44):
		if (ob.wvlRangeTable[y][x][1] < wvlUpperCutoff) or (ob.wvlRangeTable[y][x][0] > wvlLowerCutoff):
		    bad_solution_mask[y][x] = 1

        unix = ob.getFromHeader('unixtime')
        startJD = unix/86400.+2440587.5
        nSecInFile = ob.getFromHeader('exptime')
        deadMask = ob.getDeadPixels()

        y_values, x_values = np.where(np.logical_and(bad_solution_mask == 0, np.logical_and(apertureMask != 0, deadMask == 1)))
        y_sky, x_sky = np.where(np.logical_and(bad_solution_mask == 0, np.logical_and(skyMask != 0, deadMask == 1)))
	#tic = time()
        for sec in np.arange(0,nSecInFile,integrationTime):
            jd = startJD + sec/(24.*3600.) + integrationTime/2./(24.*3600.)#add seconds offset to julian date, move jd to center of bin
            print count,jd
            count+=1
            times.append(jd)
            titles.append('%.6f'%jd)
            frameData = ob.getAperturePixelCountImage(firstSec=sec, integrationTime=integrationTime, y_values=y_values, x_values=x_values, y_sky=y_sky, x_sky=x_sky, apertureMask=apertureMask, skyMask=skyMask, weighted=True)
#            frameData = ob.getPixelCountImage(firstSec=sec,integrationTime=integrationTime,weighted=True)
            backgroundPerPixel.append(frameData['SkyCountSubtractedPerPixel'])
            lightcurve = frameData['lightcurve']
	    frame = frameData['image']         
            showFrame = np.array(frame)
            showframes.append(showFrame)
            frame[deadMask == 0] = np.nan
            frame[bad_solution_mask == 1] = np.nan
            frame[apertureMask == 0] = np.nan
            frames.append(frame)
            LightCurve.append(lightcurve)
          
cube = np.dstack(frames)
times = np.array(times)
print 'average sky count', np.average(backgroundPerPixel)
print 'standard deviation of sky count', np.std(backgroundPerPixel)
np.savez(npzFileName,stack=cube,jd=times,lightcurve=LightCurve,AveSkyCount=np.average(backgroundPerPixel),StdSkyCout=np.std(backgroundPerPixel))
print 'saved'
utils.makeMovie(showframes,cbar=True,frameTitles=titles,outName=gifFileName)

