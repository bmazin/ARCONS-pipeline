#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program opens a series of observations of the crab pulsar.  It calculates the period, and calcluates the phase for every photon in the series of observations, It then plots a light curve, and a histogram of counts per period
'''
from util.ObsFile import ObsFile
from util.FileName import FileName
from util.popup import PopUp
from util import utils
import matplotlib.pyplot as plt
import numpy as np
import datetime
import tables
import ephem
import matplotlib
import matplotlib.cm as cm
import os
import hotpix.hotPixels as hp

rad=5
def circ(startpx,startpy,radius=rad):
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
                pixx.append(x)
                pixy.append(y)
    return pixx,pixy


def main():

    obsSequence0="""
    051516
    052520
    """
    obsSequence1="""
    033323
    041843
    045902
    """

    obsSequence2="""
    050404
    054424
    """

    obsSequence3="""
    054926
    062942
    """

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

    for iSeq,obsSequence in enumerate(obsSequences):
        obsSequence = obsSequence.strip().split()
        print obsSequence
        for iOb,obs in enumerate(obsSequence):
            timeMaskFileName = timeMaskFileNames[iSeq][iOb]
            if not os.path.exists(timeMaskFileName):
                print 'Running hotpix for ',obs
                hp.findHotPixels(obsFileNames[iSeq][iOb],timeMaskFileName)
                print "Flux file pixel mask saved to %s"%(timeMaskFileName)

    apertureRadius = 4
    obLists = [[ObsFile(fn)for fn in seq ] for seq in obsFileNames]
    tstampFormat = '%H:%M:%S'
    #print 'fileName','headerUnix','headerUTC','logUnix','packetReceivedUnixTime'
    for iSeq,obList in enumerate(obLists):
        for iOb,ob in enumerate(obList):
            print ob.fileName
            centerRow = centersRow[iSeq]
            centerCol = centersCol[iSeq]
            circCol,circRow = circ(centerCol,centerRow)
            ob.loadTimeAdjustmentFile(FileName(run='PAL2012').timeAdjustments())
            ob.loadWvlCalFile(wvlFileNames[iSeq])
            ob.loadFlatCalFile(flatFileNames[iSeq])
            ob.loadFluxCalFile(fluxFileNames[iSeq])
            timeMaskFileName = timeMaskFileNames[iSeq][iOb]
            ob.loadHotPixCalFile(timeMaskFileName)
            ob.setWvlCutoffs(None,None)

    for iSeq,obList in enumerate(obLists):
        for iOb,ob in enumerate(obList):
            print ob.fileName
            
            centerRow = centersRow[iSeq]
            centerCol = centersCol[iSeq]
            circCol,circRow = circ(centerCol,centerRow)
            imgDict = ob.getPixelCountImage()
            img = imgDict['image']
            utils.plotArray(img,showMe=False)
            aperture = plt.Circle((centerCol,centerRow),rad,fill=False,color='g')
            aperture2 = plt.Circle((centerCol,centerRow),2*rad,fill=False,color='g')
            plt.gca().add_patch(aperture)
            plt.gca().add_patch(aperture2)
            plt.show()
            

    
if __name__=='__main__':
    main()

