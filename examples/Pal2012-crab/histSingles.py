#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program opens a series of observations of the crab pulsar.  It calculates the period, and calcluates the phase for every photon in the series of observations, It then plots a light curve, and a histogram of counts per period
'''
from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
from util.FileName import FileName
from util.popup import PopUp
from hotpix import hotPixelsMatt as hotPixels
import matplotlib.pyplot as plt
import numpy as np
import datetime
import tables
import ephem
import matplotlib
import os

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

def calculatePeriod(mjd,secOffset=0):
    #goldstone ephemeris
    F0=29.7414493465249770#Hz
    deltaF0=0.0000000055983574
    F1=-3.669005118205e-10#df0/dt Hz/s
    F2=-3.085573120457e-20#d^2f0/dt^2 Hz/s^2
    pEpoch = 54781.604891 #Modified Julian Date corresponding to F0
    #pEpoch = pEpoch#convert mjd to jd
    pEpoch *= 24*3600 #in seconds
    startTime = mjd*24*3600#in seconds
    dt = startTime-pEpoch+secOffset#seconds since pepoch

    #freq = F0+F1*dt+F2/2*dt**2
    freq = F0+F1*dt
    period = 1.0/freq
    return period

def calculateFreq(mjd,secOffset=0):
    #goldstone ephemeris
    F0=29.7414493465249770#Hz
    deltaF0=0.0000000055983574
    F1=-3.669005118205e-10#df0/dt Hz/s
    F2=-3.085573120457e-20#d^2f0/dt^2 Hz/s^2
    pEpoch = 54781.604891 #Modified Julian Date corresponding to F0
    #pEpoch = pEpoch#convert mjd to jd
    pEpoch *= 24*3600 #in seconds
    startTime = mjd*24*3600#in seconds
    dt = startTime-pEpoch+secOffset#seconds since pepoch

    #freq = F0+F1*dt+F2/2*dt**2
    freq = F0+F1*dt
    period = 1.0/freq
    return freq

def main():
    overwriteOutFile = False
    outFile = 'mjdTimesSeqStableObs.npy'

    obsSequence="""
    060934
    """
#    obsSequence="""
#    055428
#    055930
#    060432
#    060934
#    061436
#    061938
#    062440
#    062942
#    """
    obsSequence = obsSequence.strip().split()
    obsUtcDate = '20121212'
    obsSequence = [obsUtcDate+'-'+ts for ts in obsSequence]

    obsFileNames = [FileName(run='PAL2012',date='20121211',tstamp=ts).obs() for ts in obsSequence]
    wvlFileName = FileName(run='PAL2012',date='20121211',tstamp='20121212-063518').calSoln()

    obList = [ObsFile(fn) for fn in obsFileNames]
    for iOb,ob in enumerate(obList):
        ob.loadWvlCal(wvlFileName)
        ob.setWvlCutoffs(4000,6000)

    #obsDate = obList[0].getFromHeader('jd')

    unixEpochJD = 2440587.5
    epochMJD = 2400000.5
    unixSecsInDay = 86400.
    headerUnixtime = obList[0].getFromHeader('unixtime')
    #Time of beginning of observation in mjd
    obsDate = headerUnixtime/unixSecsInDay+unixEpochJD



    #period=0.03367660643405371
    #period=0.03367664238573182
    #period=0.03367662440988317

    iRow = 10
    iCol = 14
    integrationTime = -1
    circCol,circRow = circ(iCol,iRow,radius=4)
    firstSec = 0
    
    F0=29.7414493465249770#Hz
    F1=-3.669005118205e-10#df0/dt Hz/s

    pdot = -F1/F0**2#s/s
    print 'pdot',pdot
    period0 = calculatePeriod(obsDate-epochMJD,firstSec)
    freq0 = calculateFreq(obsDate-epochMJD,firstSec)*24.*3600
    freqDot = F1*(24.*3600)**2#convert to days
    print obsDate,'period=',period0,'s'

    nPhaseBins = 200
    #np.set_printoptions(threshold=np.nan)
    if os.path.exists(outFile) and overwriteOutFile == False:
        mjdTimes = np.load(outFile)
    else:
        
        jdTimes = np.array([],dtype=np.float64)
        times = np.array([])
        for iOb,ob in enumerate(obList):
            print iOb,'of',len(obList)
            ob.setWvlCutoffs(4000,6000)
            #obsDate = ob.getFromHeader('jd')
            #Time of beginning of observation in mjd
            obsDate = ob.getFromHeader('unixtime')/unixSecsInDay+unixEpochJD
            for i in range(len(circCol)):
                iRow = circRow[i]
                iCol = circCol[i]
                timestamps,peaks,baselines = ob.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
                timestamps = np.array(timestamps,dtype=np.float64)
                jdTimestamps = obsDate+timestamps /(24.*3600.)
                jdTimes = np.append(jdTimes,jdTimestamps)
                times = np.append(times,timestamps)

        mjdTimes = jdTimes - epochMJD #convert to modified jd

        np.save(outFile,mjdTimes)
    mjdTime0 = obsDate-epochMJD
    #periodDays = period/(24.*3600.)
    periods = calculatePeriod(mjdTimes,0)
    avgPeriods = (periods+period0)/2.
    diffPeriods = periods-period0

    periodDays = periods/(24.*3600)
    periodDays0 = period0/(24.*3600)
    avgPeriodDays = avgPeriods/(24.*3600)
    #periodDays[0:len(periodDays)//2] = periodDays[0]
    #periodDays[len(periodDays)//2:] = periodDays[-1]

    phaseOffset = 0
    #phases = (mjdTimes % avgPeriodDays)/periodDays+phaseOffset
    #phases = np.log(1+pdot/period0*(mjdTimes-(obsDate-epochMJD))*(24.*3600))/pdot
    #phases = np.log(1+pdot/periodDays0*(mjdTimes-(obsDate-epochMJD)))/pdot
    #phases = mjdTimes/periodDays0-pdot*mjdTimes**2/(2*periodDays0**2)
    #totalPhases = freq0*(mjdTimes-mjdTime0)+1.65e2*freqDot*(mjdTimes-mjdTime0)**2/2.+phaseOffset
    totalPhases = freq0*(mjdTimes-mjdTime0)+freqDot*(mjdTimes-mjdTime0)**2/2.+phaseOffset
    phases0 = ((mjdTimes-mjdTime0) % avgPeriodDays[0])/periodDays+phaseOffset
    #in case phaseOffset pushed a phase past 1. roll it back 
    periodIndices=np.array(totalPhases,dtype=np.int)
    phases = totalPhases % 1.0 
    #determine which period each mjd timestamp falls in
    #make folded phase hist
    phases = phases[0:len(phases)]
    histPhases,phaseBinEdges = np.histogram(phases,bins=nPhaseBins)
    histPhases0,phaseBinEdges = np.histogram(phases0,bins=nPhaseBins)
       

    plt.plot(phaseBinEdges[0:-1],histPhases)
    plt.plot(phaseBinEdges[0:-1],histPhases0+1,'r')
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
