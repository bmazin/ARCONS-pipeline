#!/bin/python
'''
-------------
Adapted from code by Matt to write photon lists out
for the Palomar Crab data - see notes below. May 21 2013, JvE.
-------------

Author: Matt Strader        Date: March 6,2013

This program opens a series of observations of the crab pulsar.  It calculates
the period, and calcluates the phase for every photon in the series of
observations, It then plots a light curve, and a histogram of counts per period
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
import matplotlib.pylab as mpl
import os
import hotpix.hotPixels as hp
import astrometry.CentroidCalc as cc

def main():

    #Try remapping the pixels - set to None to not do this....
    pixRemapFileName = FileName(run='PAL2012').pixRemap()

    obsSequence0="""
    051516
    052018
    052520
    """
    obsSequence1="""
    033323
    033825
    034327
    034830
    035332
    035834
    040336
    040838
    041341
    041843
    042346
    042848
    043351
    043853
    044355
    044857
    045359
    045902
    """

    obsSequence2="""
    050404
    050907
    051409
    051912
    052414
    052917
    053419
    053922
    054424
    """

    obsSequence3="""
    054926
    055428
    055930
    060432
    060934
    061436
    061938
    062440
    062942
    """

    run = 'PAL2012'
    obsSequences = [obsSequence0,obsSequence1,obsSequence2,obsSequence3]
    
    #TEMPORARY FUDGE TO JUST LOOK AT A COUPLE OF IMAGES...
    obsSequences = ['9999999',
                    '''
                    033323
                    ''',
                    '9999999',
                    '9999999']
    
    wvlCals = ['051341','063518','063518','063518']
    flatCals = ['20121211','20121211','20121211','20121211']
    fluxCalDates = ['20121206','20121206','20121206','20121206']
    fluxCals = ['20121207-072055','20121207-072055','20121207-072055','20121207-072055']

    #Row coordinate of center of crab pulsar for each obsSequence
    centersRow = [9,30,29,10]
    #Col coordinate of center of crab pulsar for each obsSequence
    centersCol = [13,30,30,14]

    centerRA = '05:34:31.93830'
    centerDec = '+22:00:52.1758'

    obsUtcDate = '20121212'
    obsUtcDates = ['20121206','20121212','20121212','20121212']

    obsFileNames = []
    obsFileNameTimestamps = []
    wvlFileNames = []
    flatFileNames = []
    fluxFileNames = []
    timeMaskFileNames = []
    centroidFileNames = []
    
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        obsSequence = obsSequence.strip().split()
        obsFileNameTimestamps.append(obsSequence)
        obsUtcDate = obsUtcDates[iSeq]
        sunsetDate = str(int(obsUtcDate)-1)
        obsSequence = [obsUtcDates[iSeq]+'-'+ts for ts in obsSequence]
        obsFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).obs() for ts in obsSequence])
        timeMaskFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).timeMask() for ts in obsSequence])
        centroidFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).centroidList() for ts in obsSequence])

        wvlCalTstamp = obsUtcDate+'-'+wvlCals[iSeq]
        wvlFileNames.append(FileName(run=run,date=sunsetDate,tstamp=wvlCalTstamp).calSoln())
        fluxFileNames.append(FileName(run=run,date=fluxCalDates[iSeq],tstamp=fluxCals[iSeq]).fluxSoln())
        flatFileNames.append(FileName(run=run,date=flatCals[iSeq],tstamp='').flatSoln())
        centroidFileNames.append(FileName(run=run,date=sunsetDate,tstamp=ts).centroidList())


    for iSeq,obsSequence in enumerate(obsSequences):
        obsSequence = obsSequence.strip().split()
        print obsSequence
        for iOb,obs in enumerate(obsSequence):
            timeMaskFileName = timeMaskFileNames[iSeq][iOb]
            if os.path.exists(obsFileNames[iSeq][iOb]) and not os.path.exists(timeMaskFileName):
                print 'Running hotpix for ',obs
                hp.findHotPixels(obsFileNames[iSeq][iOb],timeMaskFileName)
                print "Flux file pixel mask saved to %s"%(timeMaskFileName)
            else:
                print 'Skipping hot pixel mask creation for non-existent file '+obsFileNames[iSeq][iOb]


    apertureRadius = 4
    obLists = [[ObsFile(fn)for fn in seq if os.path.exists(fn)] for seq in obsFileNames]
    tstampFormat = '%H:%M:%S'
    #print 'fileName','headerUnix','headerUTC','logUnix','packetReceivedUnixTime'
   
    print '---------Getting centroids-----------' 
    for iSeq,obList in enumerate(obLists):
        for iOb,ob in enumerate(obList):
            timeAdjFileName = FileName(run='PAL2012').timeAdjustments()
            wvlCalFileName = wvlFileNames[iSeq]
            flatCalFileName = flatFileNames[iSeq]
            fluxCalFileName = fluxFileNames[iSeq]
            timeMaskFileName = FileName(obsFile=ob).timeMask() #timeMaskFileNames[iSeq][iOb]
            centroidFileName = FileName(obsFile=ob).centroidList()
            if not os.path.exists(centroidFileName):
                print 'Loading calibration files:'
                print [os.path.basename(x) for x in [timeAdjFileName,wvlCalFileName,flatCalFileName, \
                        fluxCalFileName,timeMaskFileName, centroidFileName]]
                ob.loadTimeAdjustmentFile(timeAdjFileName)
                ob.loadWvlCalFile(wvlCalFileName)
                ob.loadFlatCalFile(flatCalFileName)
                ob.loadFluxCalFile(fluxCalFileName)
                ob.loadHotPixCalFile(timeMaskFileName)
                ob.setWvlCutoffs(None,None)
                print 'Running CentroidCalc for ',ob.fileName
                cc.centroidCalc(ob, centerRA, centerDec, outputFileName=centroidFileName,
                                guessTime=300, integrationTime=30, secondMaxCountsForDisplay=500,
                                xyapprox=[centersCol[iSeq],centersRow[iSeq]])
                print "Centroiding calculations saved to %s"%(centroidFileName)
            else:
                print "File already exists - skipping: "+centroidFileName
                
    print '---------Writing photon lists----------'
    for iSeq,obList in enumerate(obLists):
        for ob in obList:
#            ob.loadTimeAdjustmentFile(FileName(run='PAL2012').timeAdjustments())
#            ob.loadWvlCalFile(wvlFileNames[iSeq])
#            ob.loadFlatCalFile(flatFileNames[iSeq])
#            ob.loadFluxCalFile(fluxFileNames[iSeq])
#            timeMaskFileName = timeMaskFileNames[iSeq][iOb]
#            ob.loadHotPixCalFile(timeMaskFileName)
#            ob.setWvlCutoffs(None,None)
#                        timeAdjFileName = FileName(run='PAL2012').timeAdjustments()
            photListFileName = FileName(obsFile=ob).photonList()
            if os.path.exists(photListFileName):
                print 'Phot. list file already exists - skipping: ',photListFileName
            else:
                wvlCalFileName = wvlFileNames[iSeq]
                flatCalFileName = flatFileNames[iSeq]
                fluxCalFileName = fluxFileNames[iSeq]
                timeMaskFileName = FileName(obsFile=ob).timeMask() #timeMaskFileNames[iSeq][iOb]
                centroidFileName = FileName(obsFile=ob).centroidList()
                print 'Loading calibration files:'
                print [os.path.basename(x) for x in [timeAdjFileName,wvlCalFileName,flatCalFileName, \
                      fluxCalFileName,timeMaskFileName, centroidFileName]]
                ob.loadTimeAdjustmentFile(timeAdjFileName)
                ob.loadWvlCalFile(wvlCalFileName)
                ob.loadFlatCalFile(flatCalFileName)
                ob.loadFluxCalFile(fluxCalFileName)
                ob.loadHotPixCalFile(timeMaskFileName)
                ob.loadCentroidListFile(centroidFileName)
                ob.setWvlCutoffs(None,None)
                
                #Show the image in detector space
                #ob.displaySec()
                #Mark any pixels that were bad at any point:
                #badPix = hp.getHotPixels(ob.hotPixTimeMask)
                #x = np.arange(ob.nCol)
                #y = np.arange(ob.nRow)
                #xx, yy = np.meshgrid(x, y)
                #if np.sum(badPix) > 0: mpl.scatter(xx[badPix], yy[badPix], c='y')
                
                print 'Writing: '+FileName(obsFile=ob).photonList()
                ob.writePhotonList(pixRemapFileName=pixRemapFileName)

    print '-----------All done--------------'

if __name__=='__main__':
    main()

