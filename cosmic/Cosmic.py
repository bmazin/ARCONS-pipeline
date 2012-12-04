#!/bin/python
'''Author:  Chris Stoughton    Date:  November 28, 2012

Identify synchronous photons
'''

import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.utils import confirm
from headers import ArconsHeaders
from util import utils
from util.ObsFile import ObsFile

class Cosmic:

    def __init__(self, fileName, beginTime=0, endTime='exptime'):
        """
        Opens fileName in MKID_DATA_DIR, sets roachList
        """
        self.fileName = fileName
        self.file = ObsFile(fileName)
        self._setRoachList()
        self._setAllSecs()
        self.exptime = self.file.getFromHeader('exptime')
        if endTime =='exptime':
            self.endTime = self.exptime
        else:
            self.endTime = endTime
        if (self.endTime > self.exptime or endTime < 0):
            raise RuntimeError("bad endTime:  endTime=%d exptime=%d" % \
                                   (endTime,self.exptime))

        self.beginTime = beginTime
    def __del__(self):
        """
        Close any open files
        """
        try:
            self.file.close()
        except AttributeError:
            pass

    def _setRoachList(self):
        self.roachList = []
        for row in self.file.beamImage:
            for name in row:
                roachName = name.split("/")[1]
                if roachName not in self.roachList:
                    self.roachList.append(roachName)
        self.roachList.sort()

    def _setAllSecs(self):
       nRow = self.file.nRow
       nCol = self.file.nCol
       self.allSecs =  \
           dict( ((i,j),None) for i in range(nRow) for j in range(nCol))
       for iRow in np.arange(nRow):
           for iCol in np.arange(nCol):
               self.allSecs[iRow,iCol] = \
                   self.file.file.getNode(self.file.beamImage[iRow][iCol])

    def nPhoton(self, beginTime=0, endTime='expTime'):
        nPhoton = 0
        for iRow in range(self.file.nRow):
            for iCol in range(self.file.nCol):
                for iSec in range(self.beginTime, self.endTime):
                    sec = self.allSecs[iRow,iCol][iSec]
                    nPhoton += len(sec)
        return nPhoton

    def makeTimeHgs(self, nBinsPerSec=10, beginTime=0, endTime="exptime"):
        """
        Fill in the timeHgs variable
        This is a dictionary, indexed by the roach name, of the time histograms
        """
        self.timeHgs = {}
        self.nBinsPerSec = nBinsPerSec
        

        self.times = \
        np.arange(self.beginTime, self.endTime, 1.0/self.nBinsPerSec)

        for iSec in range(self.beginTime, self.endTime):
            print "Cosmic.makeTimeHgs:  iSec=%4d / %4d" % (iSec,self.endTime)
            hgsThisSec = {}
            for iRow in range(self.file.nRow):
                for iCol in range(self.file.nCol):
                    sec = self.allSecs[iRow,iCol][iSec]
                    if len(sec) > 0:
                        times = sec & self.file.timeMask
                        hg,edges = \
                        np.histogram(times,bins=nBinsPerSec, \
                                     range=(0,1.0/self.file.tickDuration))
                        roachName = \
                        self.file.beamImage[iRow][iCol].split("/")[1]
                        if not hgsThisSec.has_key(roachName):
                            hgsThisSec[roachName] = \
                            np.zeros(nBinsPerSec,dtype=np.int64)
                        hgsThisSec[roachName] += hg
            for roachName in hgsThisSec.keys():
                if not self.timeHgs.has_key(roachName):
                    self.timeHgs[roachName] = []
                self.timeHgs[roachName] += list(hgsThisSec[roachName])

    def plotTimeHgs(self):
        """
        Plot the time HGS in plt structure, with legend
        """
        plt.clf()
        keys = self.timeHgs.keys()
        keys.sort()
        for roachName in keys:
            hg = self.timeHgs[roachName]
            plt.plot(self.times, hg,label=roachName)
        plt.legend()
        plt.xlabel("time (sec)")
        dt = 1.0/self.nBinsPerSec
        plt.ylabel("photons/%.2f sec" % dt)
        plt.title("Cosmic timeHgs "+ self.fileName)
