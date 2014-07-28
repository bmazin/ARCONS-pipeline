#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program opens a series of photon lists for a pulsar that have been timed by tempo2.  It uses the values from tempo2 to make a pulse profile in various wavelength bands.
'''
from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
from hackedPipeline.photlist import PhotList
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
import subprocess
import numexpr

def plotHist(ax,histBinEdges,hist,**kwargs):
    ax.plot(histBinEdges,np.append(hist,hist[-1]),drawstyle='steps-post',**kwargs)

def main():

    obsSequence0="""
    112731
    113234
    113737
    114240
    114743
    115246
    """

    obsSequence1="""
    115749
    120252
    120755
    121258
    121801
    122304
    """

    run = 'PAL2012'
    obsSequences = [obsSequence0]
    wvlCals = ['112506']
    flatCals = ['20121211']
    fluxCalDates = ['20121206']
    fluxCals = ['20121207-072055']
    bApplyWeight = False
    if bApplyWeight:
        weightLabel = 'weighted'
    else:
        weightLabel = 'unweighted'

    #Row coordinate of center of crab pulsar for each obsSequence
    centersRow = [30]
    #Col coordinate of center of crab pulsar for each obsSequence
    centersCol = [30]


    obsUtcDates = ['20121208']
    path = '/Scratch/dataProcessing/psr0656/'

    obsFileNames = []
    obsFileNameTimestamps = []
    wvlFileNames = []
    flatFileNames = []
    fluxFileNames = []
    timeMaskFileNames = []
    plFileNames = []

    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        obsSequence = obsSequence.strip().split()
        obsFileNameTimestamps.append(obsSequence)
        obsUtcDate = obsUtcDates[iSeq]
        sunsetDate = str(int(obsUtcDate)-1)
        obsSequence = [obsUtcDates[iSeq]+'-'+ts for ts in obsSequence]
        plFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).timedPhotonList() for ts in obsSequence])

    np.set_printoptions(precision=11,threshold=np.nan)

    nBins = 10

    histStart = 0.
    histEnd = 1.
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
    #wvlBandEdges = np.array([4000,7000,9500])
    wvlBandEdges = np.array([4000,11000])
    
    nWvlBands = len(wvlBandEdges)-1
    #outFile = tables.openFile(outFilePath,mode='w')
    avgProfiles = np.zeros((nWvlBands,nBins),dtype=np.double)
    nPhotons = np.zeros(nWvlBands)
    for iPL,pl in enumerate(plList):
        print iPL,
        photons = pl.photTable.read()
        print len(photons),'photons in all wvl'
        for iWvlBand,wvlBandEdge in enumerate(wvlBandEdges[1:]):
            photonsInBand = photons[photons['wavelength']<wvlBandEdge]
            phases = photonsInBand['phase']
            if bApplyWeight:
                totalWeights = photonsInBand['flatWeight']*photonsInBand['fluxWeight']
            else:
                totalWeights = None
            profile,phaseBinEdges = np.histogram(phases,bins=nBins,range=(histStart,histEnd),weights=totalWeights)
            nPhotons[iWvlBand]+=len(photonsInBand)#np.sum(profile)
            avgProfiles[iWvlBand]+=profile
                
            #remove wavelengths from the band for the next iteration
            photons = photons[photons['wavelength']>=wvlBandEdge]

    
    print 'nPhotons by band:',nPhotons,', per bin:',nPhotons/nBins
    print 'expected std of band',1./np.sqrt(nPhotons/nBins)
    avgProfiles = avgProfiles / np.mean(avgProfiles,axis=1)[:,np.newaxis]
    print 'std of band',np.std(avgProfiles,axis=1)
    np.savez(path+'fullProfile{}bins_{}bands_{}.npz'.format(nBins,nWvlBands,weightLabel),avgProfiles=avgProfiles,phaseBinEdges=phaseBinEdges,wvlBandEdges=wvlBandEdges)    
    print 'done saving'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    phaseBinCenters = phaseBinEdges[0:-1]+np.diff(phaseBinEdges)/2.
    doublePhaseBinCenters = np.append(phaseBinCenters,phaseBinCenters+1)
    doublePhaseBinEdges = np.append(phaseBinEdges[0:-1],phaseBinEdges+1)
    
    for iWvlBand,wvlBand in enumerate(wvlBandEdges[1:]):
        doubleProfile = np.append(avgProfiles[iWvlBand],avgProfiles[iWvlBand])
        color=cm.jet((iWvlBand+1.)/len(wvlBandEdges))
        plotHist(ax=ax,histBinEdges=doublePhaseBinEdges,hist=doubleProfile,color=color)
        #ax.plot(phaseBinEdges[0:-1],avgProfiles[iWvlBand],color=color)
        #ax.plot(doublePhaseBinCenters,doubleProfile,color=color)

    plt.show()

if __name__ == '__main__':
    main()
