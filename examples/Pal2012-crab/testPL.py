#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

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
    """

    verbose = True
    run = 'PAL2012'
    path = '/Scratch/dataProcessing/crabData/'
    nIdxToCheck = 81
    nBins = 250
    outFilePath = path+'indPulseProfiles_3sigma_%d_%dphaseBins_sky.h5'%(nIdxToCheck,nBins)

    obsSequences = [obsSequence1,obsSequence2,obsSequence3]
    #obsSequences = [obsSequence1]
    wvlCals = ['063518','063518','063518']
    flatCals = ['20121211','20121211','20121211']
    fluxCalDates = ['20121206','20121206','20121206']
    fluxCals = ['20121207-072055','20121207-072055','20121207-072055']

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
        plFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).crabList() for ts in obsSequence])
        skyFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).crabSkyList() for ts in obsSequence])

    np.set_printoptions(precision=11,threshold=np.nan)


    labels=np.array('iBTOA   pulseNumber BTOA    Noise_Offset    Noise_RMS   Max_P2  Mean_P2 Index_P2    TestMax_P2  TestMean_P2 TestIndex_P2'.split())

    table = np.loadtxt(path+'giantPulseList_P2_3sigma_indices.txt',skiprows=1,usecols=range(len(labels)))
    peakDetectionMask = table[:,np.argmax('Max_P2'==labels)]!=0
    table = table[peakDetectionMask]#cut out test_P2 (false) detections

    dimMask = np.ones(len(table))
    #dimMask[13292:22401]=0
    dimMask = dimMask==1

    radioStrength = table[:,5]
    radioStrengthCutoff = 0.175
    radioCutoffMask = radioStrength >= radioStrengthCutoff
    strongMask = np.logical_and(radioCutoffMask,dimMask)
    table = table[strongMask]

    giantDict = dict()
    for iLabel,label in enumerate(labels):
        giantDict[label] = table[:,np.argmax(label==labels)]

    peakToCheck = 2
    noiseStrengthLabel = 'TestMax_P{}'.format(peakToCheck)
    strengthProbLabel = 'Prob_P{}'.format(peakToCheck)


    #count number of detections in the test range and the P2 peak range
    nNoiseDetections = np.sum(giantDict[noiseStrengthLabel]!=0)
    print nNoiseDetections,'noise detections'
    nPeakDetections = np.sum(giantDict['Max_P2']!=0)
    print nPeakDetections,'peak detections'
    
    giantPulseNumbers = giantDict['pulseNumber']
    radioMax = giantDict['Max_P2']
    radioMean = giantDict['Mean_P2']
    radioDetectedIndices = giantDict['Index_P2']


    nGRP = len(giantPulseNumbers)
    if verbose:
        print 'Number of GRPs',nGRP
    startIdxs = np.zeros((nGRP,nIdxToCheck))
    endIdxs = np.zeros((nGRP,nIdxToCheck))
    idxOffsets = np.array(np.linspace(-(nIdxToCheck//2),nIdxToCheck//2,nIdxToCheck),dtype=np.int)
    nIdxOffsets = len(idxOffsets)

    histStart = 0.
    histEnd = 1.
    pulsarPeriod = 33e-3 #s, approximately
    plList = [PhotList(fn) for seq in plFileNames for fn in seq]
    skyList = [PhotList(fn) for seq in skyFileNames for fn in seq]
    plMins = []
    plMaxs = []
    for iPL,pl in enumerate(plList):
        print pl.fileName,
        try:
            pns = pl.photTable.read(0)
            if pl.photTable.colindexed['pulseNumber'] == True:
                minPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][0]]
                maxPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][-1]]
                print minPulseNumber,maxPulseNumber
                plMins.append(minPulseNumber)
                plMaxs.append(maxPulseNumber)
            else:
                print 'unindexed'
        except:
            print 'corrupted',pl.fileName
    for iPL,pl in enumerate(skyList):
        print pl.fileName,
        try:
            pns = pl.photTable.read(0)
            if pl.photTable.colindexed['pulseNumber'] == True:
                minPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][0]]
                maxPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][-1]]
                print minPulseNumber,maxPulseNumber
                plMins.append(minPulseNumber)
                plMaxs.append(maxPulseNumber)
            else:
                print 'unindexed'
        except:
            print 'corrupted',pl.fileName
    plMins = np.array(plMins)
    plMaxs = np.array(plMaxs)
        


if __name__ == '__main__':
    main()
