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

def calculatePeriod(jd,secOffset):
    #goldstone ephemeris
    F0=29.7414493465249770#Hz
    deltaF0=0.0000000055983574
    F1=-3.669005118205e-10#df0/dt Hz/s
    F2=-3.085573120457e-20#d^2f0/dt^2 Hz/s^2
    pEpoch = 54781.604891 #Modified Julian Date corresponding to F0
    pEpoch = pEpoch+2400000.5#convert mjd to jd
    pEpoch *= 24*3600 #in seconds
    startTime = jd*24*3600#in seconds
    dt = startTime-pEpoch+secOffset#seconds since pepoch

    #freq = F0+F1*dt+F2/2*dt**2
    freq = F0+F1*dt
    period = 1.0/freq
    return period

def main():
    outFile = 'jdTimesSeq2.npy'
    obsSequence="""
    055428
    055930
    060432
    060934
    061436
    061938
    062440
    062942
    """
    obsSequence = obsSequence.strip().split()
    obsUtcDate = '20121212'
    obsSequence = [obsUtcDate+'-'+ts for ts in obsSequence]

    obsFileNames = [FileName(run='PAL2012',date='20121211',tstamp=ts).obs() for ts in obsSequence]

    obList = [ObsFile(fn) for fn in obsFileNames]
    obsDate = obList[0].getFromHeader('jd')



    #period=0.03367660643405371
    #period=0.03367664238573182
    #period=0.03367662440988317

    iRow = 10
    iCol = 14
    integrationTime = -1
    circCol,circRow = circ(iCol,iRow,radius=4)
    firstSec = 0
    
    period = calculatePeriod(obsDate,firstSec)
    print 'period=',period,'s'

    nPhaseBins = 200
    #np.set_printoptions(threshold=np.nan)
    
    jdTimes = np.array([],dtype=np.float64)
    times = np.array([])
    for iOb,ob in enumerate(obList):
        print iOb,'of',len(obList)
        obsDate = ob.getFromHeader('jd')
        for i in range(len(circCol)):
            iRow = circRow[i]
            iCol = circCol[i]
            timestamps,peaks,baselines = ob.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
            timestamps = np.array(timestamps,dtype=np.float64)
            jdTimestamps = obsDate+timestamps /(24.*3600.)
            jdTimes = np.append(jdTimes,jdTimestamps)
            times = np.append(times,timestamps)

    jdTimes -= 2400000.5 #convert to modified jd

    np.save(outFile,jdTimes)
    periodDays = period/(24.*3600.)
    phaseOffset = .2
    phases = (jdTimes % periodDays)/periodDays+phaseOffset
    #in case phaseOffset pushed a phase past 1. roll it back 
    phases = phases % 1.0 
    #determine which period each jd timestamp falls in
    periodIndices = (jdTimes+phaseOffset*periodDays)//periodDays
    #make folded phase hist
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
