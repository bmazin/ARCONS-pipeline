#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program opens pulsar photon lists that have been pushed through tempo2, so that the pulsar period used can be found.
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

    secInDay = 3600*24
    for iPL,pl in enumerate(plList):
        minPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][0]]
        maxPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][-1]]
        photon0 = pl.photTable[0]
        print photon0['bjd'],photon0['totalPhase']
        photon1 = pl.photTable[1]
        print photon1['bjd'],photon1['totalPhase']
        photon2 = pl.photTable[2]
        print photon2['bjd'],photon2['totalPhase']
        period0=(photon1['bjd']-photon0['bjd'])/(photon1['totalPhase']-photon0['totalPhase'])*secInDay
        print period0,'s'
        period1=(photon2['bjd']-photon0['bjd'])/(photon2['totalPhase']-photon0['totalPhase'])*secInDay
        print period1,'s'
        period2=(photon2['bjd']-photon1['bjd'])/(photon2['totalPhase']-photon1['totalPhase'])*secInDay
        print period2,'s'

        freq0 = 1/period0
        print freq0,'Hz'


if __name__ == '__main__':
    main()
