#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program opens a series of observations of the crab pulsar.  It calculates the period, and calcluates the phase for every photon in the series of observations, It then plots a light curve, and a histogram of counts per period
'''
from flatcal.flatCal import FlatCal
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
import time
import linecache

def circ(startpx,startpy,radius=3):
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

def calculateFreq(mjd,secOffset=0):
    #goldstone ephemeris
    F0=29.7414493465249770#Hz
    deltaF0=0.0000000055983574
    F1=-3.669005118205e-10#df0/dt Hz/s
    F2=-3.085573120457e-20#d^2f0/dt^2 Hz/s^2
    pEpoch = 54781.604891 #Modified Julian Date corresponding to F0
    #pEpoch = pEpoch#convert mjd to jd
    secsPerDay = (24.*3600)
    pEpoch *= secsPerDay #in seconds
    startTime = mjd*secsPerDay#in seconds
    dt = startTime-pEpoch+secOffset#seconds since pepoch

    #freq = F0+F1*dt+F2/2.*dt**2
    freq = F0+F1*dt
    return freq

def main():
    overwriteOutFile = False
    outFile = 'mjdTimesSeq2day38.npz'


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

    needPlusOne = """
    034327
    040838
    041843
    042848
    050404
    062942
    """

    needPlusOne = needPlusOne.strip().split()
    obsSequences = [obsSequence0,obsSequence1,obsSequence2,obsSequence3]
    wvlCals = ['051341','063518','063518','063518']

    #Row coordinate of center of crab pulsar for each obsSequence
    centersRow = [9,29,29,10] 
    #Col coordinate of center of crab pulsar for each obsSequence
    centersCol = [13,28,30,14] 


    obsFileNames = []
    obsUtcDate = '20121212'
    obsUtcDates = ['20121206','20121212','20121212','20121212']
    obsFileNameTimestamps = []
    wvlFileNames = []
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        obsSequence = obsSequence.strip().split()
        obsFileNameTimestamps.append(obsSequence)
        obsUtcDate = obsUtcDates[iSeq]
        sunsetDate = str(int(obsUtcDate)-1)
        obsSequence = [obsUtcDates[iSeq]+'-'+ts for ts in obsSequence]
        obsFileNames.append([FileName(run='PAL2012',date=sunsetDate,tstamp=ts).obs() for ts in obsSequence])
        wvlCalTstamp = obsUtcDate+'-'+wvlCals[iSeq]
        wvlFileNames.append(FileName(run='PAL2012',date=sunsetDate,tstamp=wvlCalTstamp).calSoln())

    obLists = [[ObsFile(fn)for fn in seq ] for seq in obsFileNames]
    tstampFormat = '%H:%M:%S'
    #print 'fileName','headerUnix','headerUTC','logUnix','packetReceivedUnixTime'
    for iSeq,obList in enumerate(obLists):
        for iOb,ob in enumerate(obList):
            ob.loadWvlCalFile(wvlFileNames[iSeq])
            ob.setWvlCutoffs(3000,8000)
            obsUtcDate = obsUtcDates[iSeq]
            sunsetDate = str(int(obsUtcDate)-1)
            pmLogFileName = FileName(run='PAL2012',date=sunsetDate,tstamp=obsUtcDate+'-'+obsFileNameTimestamps[iSeq][iOb]).packetMasterLog()
            try:
                
                #f = open('/ScienceData/PacketMasterLogs/obs_%s-%s.log'%(obsUtcDate,obsSequences[iSeq].strip().split()[iOb]),'r')
                f = open(pmLogFileName,'r')
            except:
                print 'missing Packet Master Log ',pmLogFileName
                continue
            lastTstampLines = np.zeros(8)
            firstTstampLines = np.zeros(8)
            for line in f:
                if line.split(' ')[0] == 'bundle':
                    tstampLine = line
                    break
            for line in f:
                if line.split(' ')[0] == 'bundle':
                    at = float(line.split('at')[1].split()[0])
                    lastTstampLine = at
                    iRoach = int(line.split('roach')[1].split('took')[0].strip())
                    lastTstampLines[iRoach] = at
                    if firstTstampLines[iRoach] == 0:
                        firstTstampLines[iRoach] = at
                
#                for line in f:
#                    if line.split(' ')[0] == 'bundle':
#                        lastTstampLine = line
            packetReceivedUnixTimestamp = float((tstampLine.split('took')[1].split('total')[0].strip()))
#                lastPacketReceivedUnixTimestamp = float(lastTstampLine.split('at')[1].split()[0])
            firstPacketDelay = packetReceivedUnixTimestamp-int(ob.getFromHeader('unixtime'))
            #print os.path.split(obsFileNames[iSeq][iOb])[1],time.strftime(tstampFormat,time.gmtime(int(ob.getFromHeader('unixtime')))),time.strftime(tstampFormat,time.gmtime(packetReceivedUnixTimestamp)),packetReceivedUnixTimestamp,lastPacketReceivedUnixTimestamp
            #print os.path.split(obsFileNames[iSeq][iOb])[1],time.strftime(tstampFormat,time.gmtime(int(ob.getFromHeader('unixtime')))),'%.3f'%firstPacketDelay,'%.2f'%(lastPacketReceivedUnixTimestamp),'%.2f'%(lastPacketReceivedUnixTimestamp+firstPacketDelay)
#            print os.path.split(obsFileNames[iSeq][iOb])[1],
#            for i in range(8):
#                print '%.2f'%firstTstampLines[i],
#            print ''
            
            print os.path.split(obsFileNames[iSeq][iOb])[1],
            for i in range(8):
                print '%d'%((lastTstampLines[i]+firstPacketDelay-300)//1.),
            print ''

    #obsDate = obList[0].getFromHeader('jd')

    unixEpochJD = 2440587.5
    epochMJD = 2400000.5
    secsPerDay = (24.*3600)

    #Time of beginning of observation in mjd, (use beginning of second sequence, second day)
    obsDate0 = obLists[0][0].getFromHeader('unixtime')/secsPerDay+unixEpochJD
    obsDate = obsDate0

    integrationTime = -1
    firstSec = 0
    
    F0=29.7414493465249770#Hz
    F1=-3.669005118205e-10#df0/dt Hz/s
    F2=-3.085573120457e-20#d^2f0/dt^2 Hz/s^2

    adjustF0=1#-2e-6
    #adjustF1 = -200 #for seq2
    #adjustF1 = 30 #for all seqs (adjustF0=1.)
    adjustF1 = 1
    #adjustF2 = 1e7
    adjustF2 = 0
    freq0 = adjustF0*calculateFreq(obsDate0-epochMJD,firstSec)*secsPerDay
    freqDot = adjustF1*F1*secsPerDay**2#convert to days
    freqDotDot = adjustF2*F2* secsPerDay**3

    nPhaseBins = 200
    #np.set_printoptions(threshold=np.nan)
    if os.path.exists(outFile) and overwriteOutFile == False:
        npz = np.load(outFile)
        mjdTimes = npz['mjdTimes']
        mjdTimesList = npz['mjdTimesList']
    else:
        
        jdTimes = np.array([],dtype=np.float64)
        mjdTimesList = []
        for iSeq,obList in enumerate(obLists):
            centerRow = centersRow[iSeq]
            centerCol = centersCol[iSeq]
            circCol,circRow = circ(centerCol,centerRow,radius=4)
            for iOb,ob in enumerate(obList):
                print iOb,'of',len(obList)
                #Time of beginning of observation in mjd
                obsDate = ob.getFromHeader('unixtime')/secsPerDay+unixEpochJD
                jdTimesOb = np.array([],dtype=np.float64)
                firstSec = 0
                for i in range(len(circCol)):
                    iRow = circRow[i]
                    iCol = circCol[i]
                    pixelWvlList = ob.getPixelWvlList(iRow,iCol,firstSec,integrationTime,excludeBad=True)
                    timestamps = np.array(pixelWvlList['timestamps'],dtype=np.float64)

                    jdTimestamps = obsDate+timestamps /secsPerDay
                    jdTimesOb = np.append(jdTimesOb,jdTimestamps)
                if obsFileNameTimestamps[iSeq][iOb] in needPlusOne:
                    jdTimesOb += 1./secsPerDay
                mjdTimesList.append(jdTimesOb-epochMJD)
                jdTimes = np.append(jdTimes,jdTimesOb)
        mjdTimes = jdTimes - epochMJD #convert to modified jd
        np.savez(outFile,mjdTimes=mjdTimes,mjdTimesList=mjdTimesList)
    mjdTime0 = obsDate0-epochMJD

    needPlusOneMask = []
    for iSeq,obList in enumerate(obLists):
        for iOb,ob in enumerate(obList):
            if obsFileNameTimestamps[iSeq][iOb] in needPlusOne:
                needPlusOneMask.append(1)
            else:
                needPlusOneMask.append(0)

    phaseOffset = 0
    print freq0
    print mjdTime0
    totalPhases = freq0*(mjdTimes-mjdTime0)+freqDot*(mjdTimes-mjdTime0)**2/2.+freqDotDot*(mjdTimes-mjdTime0)**3/6.+phaseOffset

    totalPhasesList = []
    histList = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    obsFileNameTimestamps = [ts for seq in obsFileNameTimestamps for ts in seq]
    for iTimes,times in enumerate(mjdTimesList):
        phases = freq0*(times-mjdTime0)+freqDot*(times-mjdTime0)**2/2.+freqDotDot*(times-mjdTime0)**3/6.+phaseOffset
        phases = phases % 1.0 
        totalPhasesList.append(phases)
        histPhases,phaseBinEdges = np.histogram(phases,bins=nPhaseBins)
        histList.append(histPhases)
        ax.plot(phaseBinEdges[0:-1],histPhases,c=cm.Paired((iTimes+1.)/len(mjdTimesList)),label=obsFileNameTimestamps[iTimes])

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
              fancybox=True, shadow=True, ncol=8)
    plt.show()
        
    #in case phaseOffset pushed a phase past 1. roll it back 
    periodIndices=np.array(totalPhases,dtype=np.int)
    phases = totalPhases % 1.0 
    #determine which period each mjd timestamp falls in
    #make folded phase hist
    phases = phases[0:len(phases)]
    histPhases,phaseBinEdges = np.histogram(phases,bins=nPhaseBins)

    plt.plot(phaseBinEdges[0:-1],histPhases)
    plt.show()

    firstPeriod = periodIndices.min()
    lastPeriod = periodIndices.max()
    nPeriods = lastPeriod-firstPeriod+1
    periodBins = np.arange(firstPeriod,lastPeriod+2)

    #Count how many timestamps/photons occur in each period
    singlePeriodCounts,periodBinEdges = np.histogram(periodIndices,bins=periodBins)
    plt.plot(periodBinEdges[:-1],singlePeriodCounts)
    plt.show()

    #make a histogram of counts in each period
    histSinglePeriodCounts,countBinEdges = np.histogram(singlePeriodCounts,bins=50)
    plt.plot(countBinEdges[:-1],histSinglePeriodCounts)
    plt.show()

if __name__=='__main__':
    main()

