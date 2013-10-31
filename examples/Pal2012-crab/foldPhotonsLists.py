#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program opens a series of observations of the crab pulsar.  It calculates the period, and calcluates the phase for every photon in the series of observations, It then plots a light curve, and a histogram of counts per period
'''
from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
from photlist import PhotList
from FileName import FileName
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
import subprocess
import numexpr


def main():

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

    verbose = True
    run = 'PAL2012'
    path = '/Scratch/dataProcessing/crabData2/'
    obsSequences = [obsSequence1,obsSequence2,obsSequence3]
    #obsSequences = [obsSequence1]
    wvlCals = ['063518','063518','063518']
    flatCals = ['20121211','20121211','20121211']
    fluxCalDates = ['20121206','20121206','20121206']
    fluxCals = ['20121207-072055','20121207-072055','20121207-072055']

    #Row coordinate of center of crab pulsar for each obsSequence
    centersRow = [29,29,10]
    #Col coordinate of center of crab pulsar for each obsSequence
    centersCol = [28,30,14]


    obsUtcDates = ['20121212','20121212','20121212']


    obsFileNames = []
    obsFileNameTimestamps = []
    wvlFileNames = []
    flatFileNames = []
    fluxFileNames = []
    timeMaskFileNames = []
    plFileNames = []
    plFileNames2 = []

    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        obsSequence = obsSequence.strip().split()
        obsFileNameTimestamps.append(obsSequence)
        obsUtcDate = obsUtcDates[iSeq]
        sunsetDate = str(int(obsUtcDate)-1)
        obsSequence = [obsUtcDates[iSeq]+'-'+ts for ts in obsSequence]
        plFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).crabList() for ts in obsSequence])

    np.set_printoptions(precision=11,threshold=np.nan)
    path = '/Scratch/dataProcessing/crabData2/'


    nBins = 33000 #1 bin = ~1 us

    histStart = 0.
    histEnd = 1.
    pulsarPeriod = 33e-3 #s, approximately
    plList = [PhotList(fn) for seq in plFileNames for fn in seq]
    plMins = []
    plMaxs = []
#    for iPL,pl in enumerate(plList):
#        minPulseNumber,maxPulseNumber = pl.file.root.pulseNumberRange.read()
#        plMins.append(minPulseNumber)
#        plMaxs.append(maxPulseNumber)

    for iPL,pl in enumerate(plList):
        minPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][0]]
        maxPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][-1]]
        print pl.fileName,minPulseNumber,maxPulseNumber
        plMins.append(minPulseNumber)
        plMaxs.append(maxPulseNumber)
    plMins = np.array(plMins)
    plMaxs = np.array(plMaxs)
        
    validWvlCutoff = 11000 #angstroms
    #wvlBandEdges = np.array([4000,5500,7000,8500,10000,11000])
    #wvlBandEdges = np.arange(4000,11001,1400)
    #wvlBandEdges = np.array([4000,6333,8667,11000])
    wvlBandEdges = np.array([4000,11000])
    
    nWvlBands = len(wvlBandEdges)-1
    #outFile = tables.openFile(outFilePath,mode='w')
    avgProfiles = np.zeros((nWvlBands,nBins),dtype=np.double)
    for iPL,pl in enumerate(plList):
        print ''
        print iPL
        photons = pl.photTable.readWhere('(waveUpperLimit > {}) & (wavelength > {})'.format(validWvlCutoff,wvlBandEdges[0]))
        print len(photons)
        for iWvlBand,wvlBandEdge in enumerate(wvlBandEdges[1:]):
            print iWvlBand
            photonsInBand = photons[photons['wavelength']<wvlBandEdge]
            phases = photonsInBand['phase']
            profile,phaseBinEdges = np.histogram(phases,bins=nBins,range=(histStart,histEnd))
            avgProfiles[iWvlBand]+=profile
            
            #remove wavelengths from the band for the next iteration
            photons = photons[photons['wavelength']>=wvlBandEdge]

    np.savez(path+'fullProfile.npz',avgProfiles=avgProfiles,phaseBinEdges=phaseBinEdges,wvlBandEdges=wvlBandEdges)    
    print 'done saving'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for iWvlBand,wvlBand in enumerate(wvlBandEdges[1:]):
        color=cm.jet((iWvlBand+1.)/len(wvlBandEdges))
        ax.plot(phaseBinEdges[0:-1],avgProfiles[iWvlBand],color=color)

    plt.show()

if __name__ == '__main__':
    main()
