#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program opens a series of observations of the crab pulsar.  It calculates the period, and calcluates the phase for every photon in the series of observations, It then plots a light curve, and a histogram of counts per period
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

    plList = [PhotList(fn) for seq in plFileNames for fn in seq]
    plMins = []
    plMaxs = []
#    for iPL,pl in enumerate(plList):
#        minPulseNumber,maxPulseNumber = pl.file.root.pulseNumberRange.read()
#        plMins.append(minPulseNumber)
#        plMaxs.append(maxPulseNumber)

    msecInDay = 3600*24*1000
    sampleRate = 20. #Hz
    dt = 1/sampleRate
    dtPerDay = (3600*24.)/dt
    
    wvlUpperLimit = 8000
    timesList = []
    for iPL,pl in enumerate(plList):
        photons = pl.photTable.read()
        photonsInBand = photons[photons['wavelength']<wvlUpperLimit]
        bjd = photonsInBand['bjd']
        sampledTimes = np.array(numexpr.evaluate('bjd*dtPerDay'),dtype=np.int)
        timesList.append(sampledTimes)

    timesList = np.concatenate(timesList)
    timesList -= np.min(timesList)
    nTimes = np.max(timesList)
    print len(timesList),nTimes,nTimes*dt,'elapsed'
    timestream = np.bincount(timesList)
    print timestream[0],timestream[-1]
    fftSize = 1000
    nFfts = nTimes//fftSize
    stop = fftSize*nFfts

    spec = np.zeros(fftSize)
    
    for iSample in xrange(0,stop,fftSize):
        spec += np.abs(np.fft.fft(timestream[iSample:iSample+fftSize]))

    freqs = np.fft.fftfreq(fftSize)[1:fftSize/2]/dt
    periods = 1./freqs
    spec = spec[1:fftSize/2]

    fig,ax = plt.subplots()
    ax.plot(freqs,spec)
    plt.show()
if __name__ == '__main__':
    main()
