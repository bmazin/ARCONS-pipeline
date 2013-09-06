#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program opens a series of observations of the crab pulsar.  It calculates the period, and calcluates the phase for every photon in the series of observations, It then plots a light curve, and a histogram of counts per period
'''
from util.ObsFile import ObsFile
from util.FileName import FileName
from util.popup import PopUp
import matplotlib.pyplot as plt
import numpy as np
import datetime
import tables
import ephem
import matplotlib
import matplotlib.cm as cm
import os
import hotpix.hotPixels as hp

def circ(startpx,startpy,radius=5):
    r = radius
    length = 2*r
    height = length
    allx = xrange(startpx-int(np.ceil(length/2.0)),startpx+int(np.floor(length/2.0))+1)
    ally = xrange(startpy-int(np.ceil(height/2.0)),startpy+int(np.floor(height/2.0))+1)
    pixx = []
    pixy = []
    for x in allx:
        for y in ally:
            if (np.abs(x-startpx))**2+(np.abs(y-startpy))**2 <= (r)**2:
                if y==31 and x==29:
                    y = 3
                    x = 17
                elif y==3 and x==17:
                    y = 31
                    x = 29
                pixx.append(x)
                pixy.append(y)
    return pixx,pixy


def halfTorus(startpx,startpy,radiusInner=5,radiusOuter=10):
    length = 2*radiusOuter
    height = length
    allx = xrange(startpx-int(np.ceil(length/2.0)),startpx+int(np.floor(length/2.0))+1)
    ally = xrange(startpy-int(np.ceil(height/2.0)),startpy+int(np.floor(height/2.0))+1)
    pixx = []
    pixy = []
    for x in allx:
        for y in ally:
            dist2 = (np.abs(x-startpx))**2+(np.abs(y-startpy))**2 
            if dist2 <= (radiusOuter)**2 and dist2 > (radiusInner)**2 and x >= startpx:
                if y==31 and x==29:
                    y = 3
                    x = 17
                elif y==3 and x==17:
                    y = 31
                    x = 29
                pixx.append(x)
                pixy.append(y)
    return pixx,pixy
def main():

    obsSequence1="""
    033323
    033825
    034327
    034830
    """

    obsSequence2="""
    """

    obsSequence3="""
    055930
    060432
    060934
    061436
    061938
    062440
    062942
    """

#    obsSequence1="""
#    035332
#    035834
#    040336
#    040838
#    041341
#    041843
#    042346
#    042848
#    043351
#    043853
#    044355
#    044857
#    045359
#    045902
#    """
#
#    obsSequence2="""
#    050404
#    050907
#    051409
#    051912
#    052414
#    052917
#    053419
#    053922
#    """
#
#    obsSequence3="""
#    054926
#    055428
#    """

    run = 'PAL2012'
    obsSequences = [obsSequence1,obsSequence2,obsSequence3]
    wvlCals = ['063518','063518','063518']
    flatCals = ['20121211','20121211','20121211']
    fluxCalDates = ['20121206','20121206','20121206']
    fluxCals = ['20121207-072055','20121207-072055','20121207-072055']

    #Row coordinate of center of crab pulsar for each obsSequence
    centersRow = [29,29,10]
    #Col coordinate of center of crab pulsar for each obsSequence
    centersCol = [29,30,14]


    obsUtcDate = '20121212'
    obsUtcDates = ['20121212','20121212','20121212']

    obsFileNames = []
    obsFileNameTimestamps = []
    wvlFileNames = []
    flatFileNames = []
    fluxFileNames = []
    timeMaskFileNames = []
    plFileNames = []
    skyFileNames = []
    
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        obsSequence = obsSequence.strip().split()
        obsFileNameTimestamps.append(obsSequence)
        obsUtcDate = obsUtcDates[iSeq]
        sunsetDate = str(int(obsUtcDate)-1)
        obsSequence = [obsUtcDates[iSeq]+'-'+ts for ts in obsSequence]
        obsFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).obs() for ts in obsSequence])
        timeMaskFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).timeMask() for ts in obsSequence])
        wvlCalTstamp = obsUtcDate+'-'+wvlCals[iSeq]
        wvlFileNames.append(FileName(run=run,date=sunsetDate,tstamp=wvlCalTstamp).calSoln())
        fluxFileNames.append(FileName(run=run,date=fluxCalDates[iSeq],tstamp=fluxCals[iSeq]).fluxSoln())
        flatFileNames.append(FileName(run=run,date=flatCals[iSeq],tstamp='').flatSoln())
        plFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).crabList() for ts in obsSequence])
        skyFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).crabSkyList() for ts in obsSequence])

    for iSeq,obsSequence in enumerate(obsSequences):
        obsSequence = obsSequence.strip().split()
        print obsSequence
        for iOb,obs in enumerate(obsSequence):
            timeMaskFileName = timeMaskFileNames[iSeq][iOb]
            if not os.path.exists(timeMaskFileName):
                print 'Running hotpix for ',obs
                hp.findHotPixels(obsFileNames[iSeq][iOb],timeMaskFileName)
                print "Flux file pixel mask saved to %s"%(timeMaskFileName)

    apertureRadius = 5
    obLists = [[ObsFile(fn)for fn in seq ] for seq in obsFileNames]
    tstampFormat = '%H:%M:%S'
    #print 'fileName','headerUnix','headerUTC','logUnix','packetReceivedUnixTime'
    for iSeq,obList in enumerate(obLists):
        for iOb,ob in enumerate(obList):
            print ob.fileName
            centerRow = centersRow[iSeq]
            centerCol = centersCol[iSeq]
            circCol,circRow = circ(centerCol,centerRow,radius=apertureRadius)
            skyCol,skyRow = halfTorus(centerCol,centerRow,radiusInner=apertureRadius,radiusOuter=2*apertureRadius)
            ob.loadTimeAdjustmentFile(FileName(run='PAL2012').timeAdjustments())
            ob.loadWvlCalFile(wvlFileNames[iSeq])
            ob.loadFlatCalFile(flatFileNames[iSeq])
            ob.loadFluxCalFile(fluxFileNames[iSeq])
            timeMaskFileName = timeMaskFileNames[iSeq][iOb]
            ob.loadHotPixCalFile(timeMaskFileName)
            ob.setWvlCutoffs(4000,11000)
            ob.writePhotonList(rowList=skyRow,colList=skyCol,filename=skyFileNames[iSeq][iOb])
            print 'wrote sky list'
            ob.writePhotonList(rowList=circRow,colList=circCol,filename=plFileNames[iSeq][iOb])
            print 'wrote crab list'


if __name__=='__main__':
    main()

