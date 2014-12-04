#!/bin/python
'''
Author: Matt Strader        Date: March 6,2013

This program takes previously created photon lists (customized crab aperture lists with empty columns for timing information) and feeds them through tempo2 to fill the jd,bjd,pulseNumber,totalPhase, and phase columns.  It then creates an index in the PLs so that they may be searched for pulseNumber efficiently. It also fills the waveOutOfRange and waveUpperLimit columns using the corresponding obs file's wavecal solution.

In the course of running, this program creates temporary text files filled with timestamps to feed to tempo2. Tempo2 makes it's own temporary output files for this program to read.
'''
from photlist import PhotList
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
import sys


if __name__ == '__main__':

    filenameList = np.recfromtxt('/home/mstrader/ARCONS-pipeline/params/crabphotonlist.txt')
    print filenameList
    path = '/Scratch/photonLists/20121211/'
    newFilenameList = ['pulsarP'+filename[1:] for filename in filenameList]
    fullPlPaths = [os.path.join(path,filename) for filename in filenameList]
    fullPulsarPlPaths = [os.path.join(path,filename) for filename in newFilenameList]
    plTimestamps = [filename[8:23] for filename in filenameList]
    parFile = 'B0531+21_LyneDM_TZRCorrect.par'

    #obsDate = obList[0].getFromHeader('jd')

    unixEpochJD = 2440587.5
    epochMJD = 2400000.5
    secsPerDay = (24.*3600)
    nRowsPerChunk = 1e6

#        self.firmwareDelay = self.timeAdjustFile.root.timeAdjust.firmwareDelay.read()[0]['firmwareDelay']
#        roachDelayTable = self.timeAdjustFile.root.timeAdjust.roachDelays
#        try:
#            self.roachDelays = roachDelayTable.readWhere('obsFileName == "%s"'%self.fileName)[0]['roachDelays']
#            self.timeAdjustFileName = os.path.abspath(timeAdjustFileName)
#            entry += np.max(self.roachDelays)
#            entry += self.firmwareDelay
#           plFile.root.header.header.cols.unixtime[0]


    for iPl,plPath in enumerate(fullPulsarPlPaths):
        print iPl,plPath

        pl = PhotList(plPath)
        plTimestamp = plTimestamps[iPl]
        print 'loaded...',

        unixtime = pl.file.root.header.header.cols.unixtime[0]
        firmwareDelay = pl.file.root.timeAdjust.timeAdjust.firmwareDelay.read()[0]['firmwareDelay']
        roachDelayTable = pl.file.root.timeAdjust.timeAdjust.roachDelays
        print 'obsFileName == "obs_{}.h5"'.format(plTimestamp)
        roachDelays = roachDelayTable.readWhere('obsFileName == "obs_{}.h5"'.format(plTimestamp))[0]['roachDelays']
        roachDelay = np.max(roachDelays)
        unixtime += roachDelay
        unixtime += firmwareDelay

        obsDate = unixtime/secsPerDay+unixEpochJD

        nPhotons = len(pl.photTable)
        nChunks = int(np.ceil(1.*nPhotons / nRowsPerChunk))
        print nPhotons,' photons in obs'
        for iChunk in xrange(nChunks):
            print 'chunk',iChunk,'of',nChunks
            iPhotonStart = iChunk*nRowsPerChunk
            iPhotonEnd = (iChunk+1)*nRowsPerChunk
            data = pl.photTable.read(iPhotonStart,iPhotonEnd)
            nPhotonsRead = len(data)
            timestamps = np.array(data['arrivalTime'],dtype=np.float64)
            jdTimestamps = obsDate+timestamps/secsPerDay - epochMJD
            tempJDPath = path+'temp_timestamps_{:.3f}'.format(time.time())
            tempBJDPath = path+'temp_bary_timestamps_{:.3f}'.format(time.time())
            np.savetxt(tempJDPath,jdTimestamps)
            strCommand = 'tempo2 -gr arcons -f {} -a {} -o {}'.format(parFile,tempJDPath,tempBJDPath)
            proc = subprocess.Popen(strCommand,shell=True,stdout=subprocess.PIPE)
            proc.wait()
            tempoOutput = np.loadtxt(tempBJDPath,usecols=[1,2])

            nPhotonsTimed = len(tempoOutput)
            if nPhotonsTimed != nPhotonsRead:
                print 'ERROR: number of photons sent to tempo2 does not match number of photons returned'
                exit(1)
            bjdTimestamps = tempoOutput[:,0]
            phases = tempoOutput[:,1]
            pl.photTable.cols.jd[iPhotonStart:iPhotonEnd] = jdTimestamps
            pl.photTable.cols.bjd[iPhotonStart:iPhotonEnd] = bjdTimestamps
            pl.photTable.cols.pulseNumber[iPhotonStart:iPhotonEnd] = np.array(phases,dtype=np.uint32)
            pl.photTable.cols.phase[iPhotonStart:iPhotonEnd] = numexpr.evaluate('phases%1.')
            os.remove(tempJDPath)
            os.remove(tempBJDPath)

        pl.file.flush()
        pl.file.close()
        del pl
        print 'closed.'
        sys.stdout.flush()

