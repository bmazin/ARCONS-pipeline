#!/bin/python
'''Author:  Chris Stoughton    Date:  November 28, 2012

Identify synchronous photons
'''

import sys,os
import tables
from tables.nodes import filenode
import numpy as np
import matplotlib.pyplot as plt
from util.utils import confirm
from headers import ArconsHeaders
from util import utils
from util.ObsFile import ObsFile
from util import meanclip
from interval import interval, inf, imath
from cosmic import TimeMask

class Cosmic:

    def __init__(self, fn, beginTime=0, endTime='exptime', \
                     nBinsPerSec=10, flashMergeTime=1.0):
        """
        Opens fileName in MKID_DATA_DIR, sets roachList
        """
        self.fn = fn
        self.fileName = fn.obs();
        self.file = ObsFile(self.fileName)
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
        self.timeHgs = "none"
        self.nBinsPerSec = nBinsPerSec
        self.flashMergeTime = flashMergeTime
        self.times = \
        np.arange(self.beginTime, self.endTime, 1.0/self.nBinsPerSec)

        # for measuring flashes, indexed by roach name
        self.rMean = {}   # mean from meanclip
        self.rSigma = {}  # sigma from meanclip
        self.rNSurvived = {} # number of survivors from meanclip
        self.rNormed = {} # (value-mean)/sigma
        self.flashInterval = {}
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
        """
        trivial example of counting the number of photons in a file
        """
        nPhoton = 0
        for iRow in range(self.file.nRow):
            for iCol in range(self.file.nCol):
                for iSec in range(self.beginTime, self.endTime):
                    sec = self.allSecs[iRow,iCol][iSec]
                    nPhoton += len(sec)
        return nPhoton

    def findFlashes(self, clipsig=3.0, maxiter=5,\
                        converge_num=0.05, verbose=0, flashsig=6):
        """
        find flashes by looking at the time histograms.  Calculate the
        mean,sigma using meanclip and the parameters clipsig, maxiter, converge_num
        A time bin has a flash if the normalized signal (measured-mean)/sigma
        is larger than flashsig.

        
        """
        # make the blank data structures
        if self.timeHgs == "none":
            self.makeTimeHgs()
            self.flashInterval["all"] = []
        # find the flashes in each roach 
        for roach in self.roachList:
            self.rMean[roach],self.rSigma[roach],self.rNSurvived[roach] = \
                meanclip.meanclip(\
                np.array(self.timeHgs[roach]), clipsig, maxiter, converge_num,\
                    verbose)
            self.rNormed[roach] = \
                (self.timeHgs[roach]-self.rMean[roach])/self.rSigma[roach]

            self.flashInterval[roach] = interval()
            prev = 0
            a = self.rNormed[roach]
            for i in range(len(a)):
                this = a[i] > flashsig
                if (this != prev):
                    if (this):
                        iBegin = i
                    else:
                        self.flashInterval[roach] = \
                        self.flashInterval[roach] | interval[iBegin,i]
                prev = this
            if (prev):
                self.flashInterval[roach] = \
                    self.flashInterval[roach] | interval[iBegin,i]

        # union of all flashes
        self.flashInterval["all"] = interval()
        for roach in self.roachList:
            self.flashInterval["all"] = \
                self.flashInterval["all"] | self.flashInterval[roach]

        # look for gaps smaller than self.flashMergeTime and plug them
        dMax = self.nBinsPerSec*self.flashMergeTime
        extrema = self.flashInterval["all"].extrema
        for i in range(len(extrema)/2-1):
            i0 = extrema[2*i+1][0]
            i1 = extrema[2*(i+1)][0]
            if (i1-i0) <= dMax:
                t0 = self.beginTime + float(i0)/self.nBinsPerSec
                t1 = self.beginTime + float(i1)/self.nBinsPerSec
                self.flashInterval["all"] = \
                    self.flashInterval["all"] | interval[i0,i1]

        # convert to ticks since the beginning of the data file
        rlAll = list(self.roachList)
        rlAll.append("all")
        ticksPerSecond = int(1.0/self.file.tickDuration)
        print "ticksPerSecond=",ticksPerSecond
        offset = self.beginTime*ticksPerSecond
        scale = 1.0/(self.file.tickDuration*self.nBinsPerSec)
        print "offset=",offset,"  scale=",scale
        for roach in rlAll:
            self.flashInterval[roach] = offset+scale*self.flashInterval[roach]

    def writeFlashesToHdf5(self,overwrite=1):
        """
        write intervals with flashes to the timeMask file
        """
        # get the output file name, and make the directory if you need to
        tmFileName = self.fn.timeMask()
        (tmDir,name) = os.path.split(tmFileName)
        if not os.path.exists(tmDir):
            os.makedirs(tmDir)

        # write parameters used to find flashes
        h5f = tables.openFile(tmFileName, 'w')
        fnode = filenode.newNode(h5f, where='/', name='timeMaskHdr')
        fnode.attrs.beginTime      = self.beginTime
        fnode.attrs.endTime        = self.endTime
        fnode.attrs.nBinsPerSec    = self.nBinsPerSec
        fnode.attrs.flashMergeTime = self.flashMergeTime
        fnode.close();

        # write the times where flashes are located
        tbl =  h5f.createTable('/','timeMask',TimeMask.TimeMask,"Time Mask")
        rlAll = list(self.roachList)
        rlAll.append("all")
        for roach in rlAll:
            extrema = self.flashInterval[roach].extrema
            for i in range(len(extrema)/2):
                row = tbl.row
                row['tBegin'] = int(extrema[2*i][0])
                row['tEnd'] = int(extrema[2*i+1][0])
                if (roach == "all"):
                    reason = "Merged Flash"
                else:
                    reason = "Flash in %s" % roach
                row['reason'] = TimeMask.timeMaskReason[reason]
                row.append()
                tbl.flush()
        tbl.close()
        h5f.close()

    def makeTimeHgs(self):
        """
        Fill in the timeHgs variable
        This is a dictionary, indexed by the roach name, of the time histograms
        """
        self.timeHgs = {}
        for iSec in range(self.beginTime, self.endTime):
            print "Cosmic.makeTimeHgs:  iSec=%4d / %4d" % (iSec,self.endTime)
            hgsThisSec = {}
            for iRow in range(self.file.nRow):
                for iCol in range(self.file.nCol):
                    sec = self.allSecs[iRow,iCol][iSec]
                    if len(sec) > 0:
                        times = sec & self.file.timeMask
                        hg,edges = \
                        np.histogram(times,bins=self.nBinsPerSec, \
                                     range=(0,1.0/self.file.tickDuration))
                        roachName = \
                        self.file.beamImage[iRow][iCol].split("/")[1]
                        if not hgsThisSec.has_key(roachName):
                            hgsThisSec[roachName] = \
                            np.zeros(self.nBinsPerSec,dtype=np.int64)
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
        plt.figure(1)
        keys = self.timeHgs.keys()
        keys.sort()

        plt.subplot(211)
        for roachName in keys:
            hg = self.timeHgs[roachName]
            plt.plot(self.times, hg,label=roachName)
        plt.legend()
        dt = 1.0/self.nBinsPerSec
        plt.ylabel("photons/%.2f sec" % dt)
        plt.title("Cosmic timeHgs "+ self.fileName)

        plt.subplot(212)
        for roachName in keys:
            plt.plot(self.times, \
                         self.rNormed[roachName],label=roachName)
        plt.xlabel("time (sec)")
        plt.ylim(-23,30)
        dt = 1.0/self.nBinsPerSec
        plt.ylabel("normalized photons/%.2f sec" % dt)

        y = -5
        x0 = self.beginTime + 0.1*(self.endTime-self.beginTime)
        xmax = plt.xlim()[1]
        rlAll = list(self.roachList)
        rlAll.append("all")
        for roach in rlAll:
            print "plot for roach=",roach
            plt.plot([x0,xmax],[y,y], linestyle=":", color="gray")
            plt.text(x0, y, roach, fontsize=8, va="center")
            extrema = self.flashInterval[roach].extrema
            for i in range(len(extrema)/2):
                t0 = (extrema[2*i][0]   - 0.5)*self.file.tickDuration
                t1 = (extrema[2*i+1][0] - 0.5)*self.file.tickDuration
                plt.plot([t0,t1],[y,y],'r', linewidth=4)
            y -= 2

    def findCosmics(self):
        print "begin findCosmics flashInterval=",self.flashInterval['all']
        tickDur = self.file.tickDuration
        interFullSec = interval[0, (1.0/tickDur)-1]
        for iSec in range(self.beginTime, self.endTime):
            inter = (self.flashInterval['all'] & \
                interval[iSec/tickDur, ((iSec+1)/tickDur)-1]) - iSec/tickDur
            print "findCosmics:  iSec=",iSec,"  inter=",inter
            hg = self.getHgForOneSec(iSec, inter)

    def getHgForOneSec(self, iSec, inter):
        timeHgValues = np.zeros(self.file.ticksPerSec)
        for iRow in range(self.file.nRow):
            for iCol in range(self.file.nCol):
                sec = self.allSecs[iRow,iCol][iSec]
                timestamps,parabolaPeaks,baselines = \
                    ObsFile.parsePhotonPackets(self.file, sec, inter,\
                                                   doParabolaFitPeaks=False,\
                                                   doBaselines = False)
        
                hg = np.histogram(timestamps,self.file.ticksPerSec,\
                                     range=(0,self.file.ticksPerSec))
                timeHgValues += hg[0]
        return timeHgValues
