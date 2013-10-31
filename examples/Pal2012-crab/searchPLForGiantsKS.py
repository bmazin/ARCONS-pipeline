#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

'''
from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
from photlist import PhotList
from FileName import FileName
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


    pulseLabel = 1 #1 for interpulse, 2 for main pulse
    verbose = True
    run = 'PAL2012'
    path = '/Scratch/dataProcessing/crabData2/'
    nIdxToCheck = 81
    nBins = 250
    outFilePath = path+'indPulseProfiles_3sigma_P{}_KS.h5'.format(pulseLabel)
    wvlBinEdges = ObsFile.makeWvlBins(wvlStart=4000,wvlStop=11000)
    wvlBinWidths = np.diff(wvlBinEdges)
    nWvlBins = len(wvlBinWidths)
    peakIdx=167#the index of the phaseBin at the main pulse peak

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


    labels=np.array('iBTOA   pulseNumber BTOA    Noise_Offset    Noise_RMS   Max  Mean Index    TestMax  TestMean TestIndex'.split())

    table = np.loadtxt(path+'giantPulseList_P{}_3sigma_indices.txt'.format(pulseLabel),skiprows=1,usecols=range(len(labels)))
    peakDetectionMask = table[:,np.argmax('Max'==labels)]!=0
    table = table[peakDetectionMask]#cut out test_P2 (false) detections

    dimMask = np.ones(len(table))
    #dimMask[13292:22401]=0
    dimMask = dimMask==1

    radioStrength = table[:,5]
    allGiantPulseNumbers = table[:,1]
    if pulseLabel == 2:
        radioStrengthCutoff = 0.175
    elif pulseLabel == 1:
        radioStrengthCutoff = 0.01

    radioCutoffMask = radioStrength >= radioStrengthCutoff
    strongMask = np.logical_and(radioCutoffMask,dimMask)
    table = table[strongMask]

    giantDict = dict()
    for iLabel,label in enumerate(labels):
        giantDict[label] = table[:,np.argmax(label==labels)]

    noiseStrengthLabel = 'TestMax'


    #count number of detections in the test range and the P2 peak range
    nNoiseDetections = np.sum(giantDict[noiseStrengthLabel]!=0)
    print nNoiseDetections,'noise detections'
    nPeakDetections = np.sum(giantDict['Max']!=0)
    print nPeakDetections,'peak detections'
    
    giantPulseNumbers = giantDict['pulseNumber']
    radioMax = giantDict['Max']
    radioMean = giantDict['Mean']
    radioDetectedIndices = giantDict['Index']

    nGRP = len(giantPulseNumbers)
    if verbose:
        print 'Number of GRPs',nGRP
    counts = np.zeros((nGRP,nIdxToCheck))
    skyCounts = np.zeros((nGRP,nIdxToCheck))

    idxOffsets = np.array(np.linspace(-(nIdxToCheck//2),nIdxToCheck//2,nIdxToCheck),dtype=np.int)
    nIdxOffsets = len(idxOffsets)

    histStart = 0.
    histEnd = 1.
    pulsarPeriod = 33e-3 #s, approximately
    indProfiles = np.zeros((nGRP,nIdxOffsets,nBins))
    skyIndProfiles = np.zeros((nGRP,nIdxOffsets,nBins))
    fullSpectra = np.zeros((nIdxOffsets,nWvlBins),dtype=np.double)
    skyFullSpectra = np.zeros((nIdxOffsets,nWvlBins),dtype=np.double)
    peakSpectra = np.zeros((nIdxOffsets,nWvlBins),dtype=np.double)
    skyPeakSpectra = np.zeros((nIdxOffsets,nWvlBins),dtype=np.double)
    plList = [PhotList(fn) for seq in plFileNames for fn in seq]
    skyList = [PhotList(fn) for seq in skyFileNames for fn in seq]
    plMins = []
    plMaxs = []
    for iPL,pl in enumerate(plList):
        minPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][0]]
        maxPulseNumber = pl.photTable.cols.pulseNumber[pl.photTable.colindexes['pulseNumber'][-1]]
        print pl.fileName,minPulseNumber,maxPulseNumber
        plMins.append(minPulseNumber)
        plMaxs.append(maxPulseNumber)

    plMins = np.array(plMins)
    plMaxs = np.array(plMaxs)
        
    pulseNumberTable = np.array([gpn+idxOffsets for gpn in giantPulseNumbers])
    giantPulseNumberMask = np.zeros(np.shape(pulseNumberTable))
    validWvlCutoff = 11000 #angstroms
    if verbose:
        print 'filling giant pulse number mask'
    for iGiantPN,giantPN in enumerate(allGiantPulseNumbers):
        if verbose and iGiantPN % 500 == 0:
            print 'mask',iGiantPN,'of',len(allGiantPulseNumbers)
        giantMatches = (pulseNumberTable == giantPN)
        giantPulseNumberMask = np.logical_or(giantPulseNumberMask,giantMatches)

    outFile = tables.openFile(outFilePath,mode='w')

    nLivePixels = []
    for iPL,pl in enumerate(plList):
        pixels = np.unique(pl.photTable.cols.xyPix[:])
        print pixels
        nPixels = len(pixels)
        nLivePixels.append(nPixels)

    nLivePixels = np.array(nLivePixels)
    nSkyLivePixels = []
    for iPL,pl in enumerate(skyList):
        pixels = np.unique(pl.photTable.cols.xyPix[:])
        print pixels
        nPixels = len(pixels)
        nSkyLivePixels.append(nPixels)
    nSkyLivePixels = np.array(nSkyLivePixels)

    pulseNumberList = np.unique(pulseNumberTable.ravel())       
    nPulseNumbers = len(pulseNumberList)                        
    pulseNumberListMask = np.ones(nPulseNumbers)                

    floatAtom = tables.Float64Atom()
    uint8Atom = tables.UInt8Atom()

    nExpectedRows=1e7
    #phaseArrays = [outFile.createEArray(outFile.root,'phases{:03d}'.format(iIdxOffset),floatAtom,(0,),title='phases at pulse number {} w.r.t. GRP'.format(idxOffset),expectedrows=nExpectedRows) for iIdxOffset,idxOffset in enumerate(idxOffsets)]
    #skyPhaseArrays = [outFile.createEArray(outFile.root,'skyPhases{:03d}'.format(iIdxOffset),floatAtom,(0,),title='phases at pulse number {} w.r.t. GRP in sky'.format(idxOffset),expectedrows=nExpectedRows) for iIdxOffset,idxOffset in enumerate(idxOffsets)]
    #wavelengthArrays = [outFile.createEArray(outFile.root,'wavelengths{:03d}'.format(iIdxOffset),floatAtom,(0,),title='wavelengths at pulse number {} w.r.t. GRP'.format(idxOffset),expectedrows=nExpectedRows) for iIdxOffset,idxOffset in enumerate(idxOffsets)]
    #skyWavelengthArrays = [outFile.createEArray(outFile.root,'skyWavelengths{:03d}'.format(iIdxOffset),floatAtom,(0,),title='wavelengths at pulse number {} w.r.t. GRP in sky'.format(idxOffset),expectedrows=nExpectedRows) for iIdxOffset,idxOffset in enumerate(idxOffsets)]
    grpPeakWavelengthArrays = outFile.createVLArray(outFile.root,'grpWavelengths',floatAtom,'photon wavelength list for GRPs at main peak',tables.Filters(complevel=1))
    nongrpPeakWavelengthArrays = outFile.createVLArray(outFile.root,'nongrpWavelengths',floatAtom,'photon wavelength list for nonGRPs at main peak',tables.Filters(complevel=1))

    for iGiantPN,giantPN in enumerate(giantPulseNumbers):
        if verbose and iGiantPN % 100 == 0:
            print iGiantPN
        
        iPL = np.searchsorted(plMaxs,giantPN)

        if iPL >= len(plMins):
            if verbose:
                print iGiantPN,'GRP not found in optical'
            continue
        if plMins[iPL] > giantPN:
            if verbose:
                print iGiantPN,'GRP not found in optical'

            continue

        if plMins[iPL] >= giantPN+idxOffsets[0] or plMaxs[iPL] <= giantPN+idxOffsets[-1]:
            if verbose:
                print iGiantPN,'optical pulses surrounding GRP not found'
            continue

        pl = plList[iPL]
        skyPL = skyList[iPL]
        #grab all photons in the pulseNumber range covered by all idxOffsets for this GRP
        pulseSelectionIndices = pl.photTable.getWhereList('({} <= pulseNumber) & (pulseNumber <= {})'.format(giantPN+idxOffsets[0],giantPN+idxOffsets[-1]))
        skyPulseSelectionIndices = skyPL.photTable.getWhereList('({} <= pulseNumber) & (pulseNumber <= {})'.format(giantPN+idxOffsets[0],giantPN+idxOffsets[-1]))
        photonsInPulseSelection = pl.photTable.readCoordinates(pulseSelectionIndices)
        skyPhotonsInPulseSelection = skyPL.photTable.readCoordinates(skyPulseSelectionIndices)

        nPulsesInSelection = len(np.unique(photonsInPulseSelection['pulseNumber']))
        nSkyPulsesInSelection = len(np.unique(skyPhotonsInPulseSelection['pulseNumber']))
        if  nPulsesInSelection < nIdxOffsets or nSkyPulsesInSelection < nIdxOffsets:
            if verbose:
                print 'only ',nPulsesInSelection,' pulses for ',iGiantPN,giantPN
            continue
        startIdx = 0
        
        nPixels = 1.*nLivePixels[iPL]
        nSkyPixels = 1.*nSkyLivePixels[iPL]
        nongrpPeakWvls = np.array([])
        for iIdxOffset,idxOffset in enumerate(idxOffsets):
            pulseNumber = giantPN+idxOffset
            if giantPulseNumberMask[iGiantPN,iIdxOffset] == False or idxOffset==0:
                photons = photonsInPulseSelection[photonsInPulseSelection['pulseNumber']==pulseNumber]
                skyPhotons = skyPhotonsInPulseSelection[skyPhotonsInPulseSelection['pulseNumber']==pulseNumber]
                phases = photons['phase']
                wavelengths = photons['wavelength']
                skyPhases = skyPhotons['phase']
                skyWavelengths = skyPhotons['wavelength']
                count = 1.*len(photons)/nPixels
                skyCount = 1.*len(skyPhotons)/nSkyPixels
                counts[iGiantPN,iIdxOffset] = count
                skyCounts[iGiantPN,iIdxOffset] = skyCount

                profile,phaseBinEdges = np.histogram(phases,bins=nBins,range=(histStart,histEnd))
                skyProfile,phaseBinEdges = np.histogram(skyPhases,bins=nBins,range=(histStart,histEnd))
                
                indProfiles[iGiantPN,iIdxOffset] = profile/nPixels
                skyIndProfiles[iGiantPN,iIdxOffset] = skyProfile/nSkyPixels

                spectrum,_ = np.histogram(wavelengths,bins=wvlBinEdges)
                spectrum = 1.0*spectrum/(nPixels*wvlBinWidths)#convert to counts per pixel per angstrom
                fullSpectra[iIdxOffset] += spectrum

                skySpectrum,_ = np.histogram(skyWavelengths,bins=wvlBinEdges)
                skySpectrum = 1.0*skySpectrum/(nSkyPixels*wvlBinWidths)#convert to counts per pixel per angstrom
                skyFullSpectra[iIdxOffset] += skySpectrum

                phasesBinned = np.digitize(phases,phaseBinEdges)-1
                peakPhaseMask = np.logical_or(phasesBinned==(peakIdx-1),phasesBinned==peakIdx)
                peakPhaseMask = np.logical_or(peakPhaseMask,phasesBinned==(peakIdx+1))
                peakWvls = wavelengths[peakPhaseMask]
                if idxOffset == 0:
                    grpPeakWavelengthArrays.append(peakWvls)
                else:
                    nongrpPeakWvls = np.append(nongrpPeakWvls,peakWvls)

                peakSpectrum,_ = np.histogram(peakWvls,bins=wvlBinEdges)
                peakSpectrum = 1.0*peakSpectrum/(nPixels*wvlBinWidths)#convert to counts per pixel per angstrom
                peakSpectra[iIdxOffset] += peakSpectrum

                skyPhasesBinned = np.digitize(skyPhases,phaseBinEdges)-1
                skyPeakPhaseMask = np.logical_or(skyPhasesBinned==(peakIdx-1),skyPhasesBinned==peakIdx)
                skyPeakPhaseMask = np.logical_or(skyPeakPhaseMask,skyPhasesBinned==(peakIdx+1))
                skyPeakWvls = skyWavelengths[skyPeakPhaseMask]
                skyPeakSpectrum,_ = np.histogram(skyPeakWvls,bins=wvlBinEdges)
                skyPeakSpectrum = 1.0*skyPeakSpectrum/(nSkyPixels*wvlBinWidths)#convert to counts per pixel per angstrom
                skyPeakSpectra[iIdxOffset] += skyPeakSpectrum
#                phaseArrays[iIdxOffset].append(phases)
#                skyPhaseArrays[iIdxOffset].append(skyPhases)
#                wavelengthArrays[iIdxOffset].append(photons['wavelength'])
#                skyWavelengthArrays[iIdxOffset].append(skyPhotons['wavelength'])
        nongrpPeakWavelengthArrays.append(nongrpPeakWvls)


    if verbose:
        print 'done searching'
    outFile.createArray(outFile.root,'counts',counts)
    outFile.createArray(outFile.root,'skyCounts',skyCounts)
    outFile.createArray(outFile.root,'idxOffsets',idxOffsets)
    outFile.createArray(outFile.root,'radioMax',radioMax)
    outFile.createArray(outFile.root,'radioMean',radioMean)
    outFile.createArray(outFile.root,'indProfiles',indProfiles)
    outFile.createArray(outFile.root,'skyIndProfiles',skyIndProfiles)
    outFile.createArray(outFile.root,'phaseBinEdges',phaseBinEdges)
    outFile.createArray(outFile.root,'giantPulseNumbers',giantPulseNumbers)
    outFile.createArray(outFile.root,'giantPulseNumberMask',giantPulseNumberMask)
    outFile.createArray(outFile.root,'pulseNumberTable',pulseNumberTable)
    outFile.createArray(outFile.root,'radioIndices',radioDetectedIndices)
    outFile.createArray(outFile.root,'nPixels',nLivePixels)
    outFile.createArray(outFile.root,'nSkyPixels',nSkyLivePixels)
    outFile.createArray(outFile.root,'fullSpectra',fullSpectra)
    outFile.createArray(outFile.root,'skyFullSpectra',skyFullSpectra)
    outFile.createArray(outFile.root,'peakSpectra',peakSpectra)
    outFile.createArray(outFile.root,'skyPeakSpectra',skyPeakSpectra)
    outFile.createArray(outFile.root,'wvlBinEdges',wvlBinEdges)

    outFile.flush()
    outFile.close()
    if verbose:
        print 'done saving'


if __name__ == '__main__':
    main()
