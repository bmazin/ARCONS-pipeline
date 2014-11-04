#!/bin/python
'''
Author: Matt Strader        Date: March 20,2014

This program takes previously created photon lists (customized pulsar lists with empty columns for timing information) and feeds them through tempo2 to fill the jd,bjd,pulseNumber,totalPhase, and phase columns.  It then creates an index in the PLs so that they may be searched for pulseNumber efficiently. It also fills the waveOutOfRange and waveUpperLimit columns using the corresponding obs file's wavecal solution.

In the course of running, this program creates temporary text files filled with timestamps to feed to tempo2. Tempo2 makes it's own temporary output files for this program to read.
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


def main():

    obsSequence0="""
    112633
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
    obsUtcDates = ['20121208']
    wvlCals = ['112506']
    flatCals = ['20121211']
    fluxCalDates = ['20121206']
    fluxCals = ['20121207-072055']

    #Row coordinate of center of pulsar for each obsSequence
    centersRow = [30]
    #Col coordinate of center of pulsar for each obsSequence
    centersCol = [30]

    path = '/Scratch/dataProcessing/psr0656/'
    parFile = 'B0656+14_LAT.par'
    #parFile = 'B0531+21.par'

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
        obsFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).obs() for ts in obsSequence])
        plFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).timedPhotonList() for ts in obsSequence])

        timeMaskFileNames.append([FileName(run=run,date=sunsetDate,tstamp=ts).timeMask() for ts in obsSequence])
        wvlCalTstamp = obsUtcDate+'-'+wvlCals[iSeq]
        wvlFileNames.append(FileName(run=run,date=sunsetDate,tstamp=wvlCalTstamp).calSoln())
        fluxFileNames.append(FileName(run=run,date=fluxCalDates[iSeq],tstamp=fluxCals[iSeq]).fluxSoln())
        flatFileNames.append(FileName(run=run,date=flatCals[iSeq],tstamp='').flatSoln())

    apertureRadius = 4
    obLists = [[ObsFile(fn) for fn in seq ] for seq in obsFileNames]
    plLists = [[PhotList(fn) for fn in seq ] for seq in plFileNames]
    tstampFormat = '%H:%M:%S'
    #print 'fileName','headerUnix','headerUTC','logUnix','packetReceivedUnixTime'
    wvlRangeTable = np.zeros([46, 44, 2])
    wvlCalData = plLists[0][0].file.root.wavecal.calsoln
    for calPixel in wvlCalData:
        if calPixel['wave_flag'] == 0:
            wvlRangeTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['solnrange']

    #obsDate = obList[0].getFromHeader('jd')

    unixEpochJD = 2440587.5
    epochMJD = 2400000.5
    secsPerDay = (24.*3600)

    
    processList = []
    
    for iSeq,(obList,plList) in enumerate(zip(obLists,plLists)):
        for iOb,(ob,pl) in enumerate(zip(obList,plList)):
            print iOb,'of',len(obList),ob.fileName
            ob.loadTimeAdjustmentFile(FileName(run='PAL2012').timeAdjustments())
            pl.loadData()
            print 'loaded...',
            obsDate = ob.getFromHeader('unixtime')/secsPerDay+unixEpochJD
            timestamps = np.array(pl.data['arrivalTime'],dtype=np.float64)
            jdTimestamps = obsDate+timestamps/secsPerDay - epochMJD
            tempJDPath = path+'temp_timestamps_{:.3f}'.format(time.time())
            tempBJDPath = path+'temp_bary_timestamps_{:.3f}'.format(time.time())
            np.savetxt(tempJDPath,jdTimestamps)
            strCommand = 'tempo2 -gr arcons -f {} -a {} -o {}'.format(parFile,tempJDPath,tempBJDPath)
            print strCommand
            proc = subprocess.Popen(strCommand,shell=True,stdout=subprocess.PIPE)
            proc.wait()
            processList.append(proc)
            tempoOutput = np.loadtxt(tempBJDPath,usecols=[1,2])
            bjdTimestamps = tempoOutput[:,0]
            phases = tempoOutput[:,1]
            pl.photTable.cols.jd[:] = jdTimestamps
            pl.photTable.cols.bjd[:] = bjdTimestamps
            pl.photTable.cols.totalPhase[:] = phases
            print 'tempo done...',
            pl.photTable.cols.pulseNumber[:] = np.array(phases,dtype=np.uint32)
            pl.photTable.cols.phase[:] = numexpr.evaluate('phases%1.')
            print 'phases done...',
            try:
                pl.photTable.cols.pulseNumber.createCSIndex()
                print 'index made...',
            except:
                print 'WARNING:couldn\'t make index (probably exists already)...'
            os.remove(tempJDPath)
            os.remove(tempBJDPath)
            for photon in pl.photTable.iterrows():
                x = photon['xPix']
                y = photon['yPix']
                wave = photon['wavelength']
                photon['waveUpperLimit']=wvlRangeTable[y,x][1]
                photon.update()
            print 'wave checked.'
            pl.file.flush()
            del pl

if __name__ == '__main__':
    main()
