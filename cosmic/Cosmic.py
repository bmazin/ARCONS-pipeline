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
from util import FileName
import inspect
from interval import interval, inf, imath
from headers import TimeMask
from cosmic import tsBinner
from scipy.optimize import curve_fit
from scipy.stats import expon
import time
import pickle
class Cosmic:

    def __init__(self, fn, beginTime=0, endTime='exptime', \
                     nBinsPerSec=10, flashMergeTime=1.0):
        """
        Opens fileName in MKID_DATA_DIR, sets roachList
        endTime is exclusive
        """
        self.fn = fn
        self.fileName = fn.obs();
        self.file = ObsFile(self.fileName)
        # apply Matt's time fix
        timeAdjustments = self.fn.timeAdjustments()
        self.file.loadTimeAdjustmentFile(timeAdjustments)
        # apply Julian's time masks
        timeMaskFile = self.fn.timeMask();
        if os.path.exists(timeMaskFile):
            self.file.loadHotPixCalFile(timeMaskFile,switchOnMask=True)
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
        #print "now in Cosmic.__del__ for ",self.fileName
        try:
            del self.file
        except:
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
        cmFileName = self.fn.cosmicMask()
        (cmDir,name) = os.path.split(cmFileName)
        if not os.path.exists(cmDir):
            os.makedirs(cmDir)

        # write parameters used to find flashes
        h5f = tables.openFile(cmFileName, 'w')
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
                        times = sec & self.file.timestampMask
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

    def findCosmics(self, stride=10, threshold=100, 
                    populationMax=2000, nSigma=5):
        """
        Find cosmics ray suspects.  Histogram the number of photons
        recorded at each timeStamp.  When the number of photons in a
        group of stride timeStamps is greater than threshold in second
        iSec, add (iSec,timeStamp) to cosmicTimeLists.  Also keep
        track of the hisogram of the number of photons per stride
        timeStamps.

        return a dictionary of 'populationHg', 'cosmicTimeLists',
        'binContents', 'timeHgValues', 'interval'  and 'frameSum'
      
        populationHg is a histogram of the number of entries per time bin.
        Used to create the poisson histogram in CosmicMaskExample.py
        
        cosmicTimeLists is a numpy array  of all the sequences that are
        suspects for cosmic rays
        
        binContents is a numpy array  of the total number of photons in each 
        bin
        
        timeHgValues is a histogram of the number of photons in each time
        interval
        
        frameSum is a two dimensional  numpy array of the number of photons
        detected by each pixel

        interval is the interval of data to be masked out

        """

        exptime = self.endTime-self.beginTime
        nBins = int(np.round(self.file.ticksPerSec*exptime+1))
        bins = np.arange(0, nBins, 1)
        timeHgValues = np.zeros(nBins, dtype=np.int64)
        frameSum = np.zeros((self.file.nRow,self.file.nCol))
        integrationTime = self.endTime - self.beginTime
        for iRow in range(self.file.nRow):
            #print "Cosmic.findCosmics:  iRow=",iRow
            for iCol in range(self.file.nCol):
                gtpl = self.file.getTimedPacketList(iRow,iCol,self.beginTime, 
                                                    integrationTime)
                timestamps = gtpl['timestamps']
                if timestamps.size > 0:
                    timestamps = \
                        (timestamps - self.beginTime)*self.file.ticksPerSec
                    # per Matt S. suggestion 2013-07-09
                    ts64 = np.round(timestamps).astype(np.uint64)
                    tsBinner.tsBinner(ts64, timeHgValues)
                    frameSum[iRow,iCol] += ts64.size
        pfthgv = Cosmic.populationFromTimeHgValues\
            (timeHgValues,populationMax,stride,threshold)
        #now build up all of the intervals in seconds
        i = interval()
        for cosmicTime in pfthgv['cosmicTimeList']:
            t0 = max(0,self.beginTime+(cosmicTime-50)/1.e6)
            t1 = min(self.endTime,self.beginTime+(cosmicTime+50)/1.e6)
            intTime = t1-t0
            tempTs = np.array([])
            for iRow in range(self.file.nRow):
                for iCol in range(self.file.nCol):
                    gtpl = self.file.getTimedPacketList(iRow, iCol, 
                                                        t0, intTime)
                    tempTs = np.append(tempTs, gtpl['timestamps'])
            mean = tempTs.mean()
            sigma = tempTs.std()
            left  = mean-nSigma*sigma
            #the cosmic events tend to have more peaks on the right side,
            #so more data is masked out in that direction
            right = mean+2*nSigma*sigma
            i = i | interval[left,right]
        retval = {}
        retval['timeHgValues'] = timeHgValues
        retval['populationHg'] = pfthgv['populationHg']
        retval['cosmicTimeList'] = pfthgv['cosmicTimeList']
        retval['binContents'] = pfthgv['binContents']
        retval['frameSum'] = frameSum
        retval['interval'] = i
        return retval


    @staticmethod
    def populationFromTimeHgValues(timeHgValues,populationMax,stride,threshold):
        """
        Rebin the timgHgValues histogram by combining stride bins.  If
        stride > 1, then bin a second time after shifting by stride/2
        Create populationHg, a histogram of the number of photons in
        the large bins.  Also, create (and then sort) a list
        cosmicTimeList of the start of bins (in original time units)
        of overpopulated bins that have more than threshold number of
        photons.

        return a dictionary containing populationHg and cosmicTimeList
        """
        popRange = (-0.5,populationMax-0.5)
        if stride==1:
            populationHg = np.histogram(\
                timeHgValues, populationMax, range=popRange)
            cosmicTimeList = np.where(timeHgValues > threshold)[0]
            binContents = np.extract(timeHgValues > threshold, timeHgValues)
        else:
            # rebin the timeHgValues before counting the populations
            length = timeHgValues.size
            remainder = length%stride
            if remainder == 0:
                end = length
            else:
                end = -remainder

            timeHgValuesTrimmed = timeHgValues[0:end]

            timeHgValuesRebinned0 = np.reshape(
                timeHgValuesTrimmed, [length/stride, stride]).sum(axis=1)
            populationHg0 = np.histogram(
                timeHgValuesRebinned0, populationMax, range=popRange)
            cosmicTimeList0 = stride*np.where(\
                timeHgValuesRebinned0 > threshold)[0]
            binContents0 = np.extract(timeHgValuesRebinned0 > threshold,
                                      timeHgValuesRebinned0)

            timeHgValuesRebinned1 = np.reshape(
                timeHgValuesTrimmed[stride/2:-stride/2], 
                [(length-stride)/stride, stride]).sum(axis=1)
            populationHg1 = np.histogram(\
                timeHgValuesRebinned1, populationMax, range=popRange)
            cosmicTimeList1 = (stride/2)+stride*np.where(\
                timeHgValuesRebinned1 > threshold)[0]
            binContents1 = np.extract(timeHgValuesRebinned1 > threshold,
                                      timeHgValuesRebinned1)

            populationHg = (populationHg0[0]+populationHg1[0],\
                                populationHg0[1])
            cosmicTimeList = np.concatenate((cosmicTimeList0,cosmicTimeList1))
            binContents = np.concatenate((binContents0, binContents1))
            args = np.argsort(cosmicTimeList)
            cosmicTimeList = cosmicTimeList[args]
            binContents = binContents[args]
            cosmicTimeList.sort()

        retval = {}
        retval['populationHg'] = populationHg
        retval['cosmicTimeList'] = cosmicTimeList
        retval['binContents'] = binContents
        return retval

    def makeMovies(self,beginTick, endTick, backgroundFrame, accumulate=False):
        tick0 = np.uint64(beginTick)
        tick1 = np.uint64(endTick)
        for iRow in range(cosmic.file.nRow):
            for iCol in range(cosmic.file.nCol):
                gtpl = self.getTimedPacketList(iRow,iCol,sec0,1)
        timestamps = gtpl['timestamps']
        timestamps *= cosmic.file.ticksPerSec
        ts64 = timestamps.astype(np.uint64)
        for ts in ts64:
            tindex = ts-t0
            try:
                listOfPixelsToMark[tindex].append((iRow,iCol))
            except IndexError:
                pass
            for tick in range(t0,t1):
                frames.append(frameSum)
                title = makeTitle(tick,t0,t1)
                titles.append(title)

                mfn0 = "m-%s-%s-%s-%s-%010d-%010d-i.gif"%(run,sundownDate,obsDate,seq,t0,t1)
                utils.makeMovie(frames, titles, outName=mfn0, delay=0.1, colormap=mpl.cm.gray,
                                listOfPixelsToMark=listOfPixelsToMark,
                                pixelMarkColor='red')

        for i in range(len(listOfPixelsToMark)-1):
            listOfPixelsToMark[i+1].extend(listOfPixelsToMark[i])

        mfn1 = "m-%s-%s-%s-%s-%010d-%010d-a.gif"%(run,sundownDate,obsDate,seq,t0,t1)
        utils.makeMovie(frames, titles, outName=mfn1, delay=0.1, colormap=mpl.cm.gray,
                        listOfPixelsToMark=listOfPixelsToMark,
                        pixelMarkColor='green')

    def fitDecayTime(self,t0Sec,lengthSec=200,plotFileName='none'):
        print "hello from fitDecayTime"
        timedPacketList = self.file.getTimedPacketList(
            iRow, iCol, sec0, lengthSec)

    
    def fitExpon(self, t0, t1):
        """
        Fit an exponential to all photons from time t0 to time t1
        t0 and t1 are in ticks, 1e6 ticks per second
        return a dictionary of:  timeStamps,fitParams,chi2

        """
         
        firstSec = int(t0/1e6)  # in seconds
        integrationTime = 1+int((t1-t0)/1e6) # in seconds
        nBins = integrationTime*1e6 # number of microseconds; one bin per microsecond
        timeHgValues = np.zeros(nBins, dtype=np.int64)
        print "firstSec=",firstSec," integrationTime=",integrationTime
        for iRow in range(self.file.nRow):
            for iCol in range(self.file.nCol):
                timedPacketList = self.file.getTimedPacketList(
                    iRow, iCol, firstSec=firstSec, 
                    integrationTime=integrationTime)
                timeStamps = timedPacketList['timestamps']
                if (len(timeStamps) > 0):
                    # covert the time values to microseconds, and
                    # make it the type np.uint64
                    # round per Matt S. suggestion 2013-07-09
                    #ts64 = (timeStamps).astype(np.uint64)
                    ts64round = np.round(timeStamps).astype(np.uint64)
                    tsBinner.tsBinner(ts64round, timeHgValues)
        
                    temp = 1e6*(timeStamps-firstSec)
                    for i in range(len(timeStamps)):
                        ts64 = ((timeStamps-firstSec)*1e6).astype(np.uint64)
                    # add these timestamps to the histogram timeHgValues
        remain0 = int(t0%1e6)
        remain1 = int(t1%1e6)
        timeHgValues = timeHgValues[remain0:remain1]
        x = np.arange(len(timeHgValues))
        y = timeHgValues
        
        xArray = np.arange(0, dtype=np.int64)
        yArray = np.arange(0, dtype=np.int64)

        for i in range(len(x)):
            if y[i] > 2:
                xArray = np.append(xArray,i)
                yArray = np.append(yArray,y[i])
        ySigma = np.sqrt(yArray)


        mean = (x*y).sum()/float(y.sum())
        bExponGuess = 1/mean
        aExponGuess = bExponGuess*timeHgValues.sum()
        cExponGuess = 0
        dExponGuess = 0
        pExponGuess = [aExponGuess, bExponGuess, cExponGuess, dExponGuess]

        bGaussGuess = mean
        avgx2 = (x*x*y).sum()/float(y.sum())        
        cGaussGuess = np.sqrt(avgx2-bGaussGuess*bGaussGuess)
        
        aGaussGuess = (timeHgValues.sum()/(cGaussGuess*np.sqrt(2*np.pi)))
        pGaussGuess = [aGaussGuess, bGaussGuess, cGaussGuess]
        
        xLimit = [bGaussGuess-4*cGaussGuess, bGaussGuess+4*cGaussGuess]
        retval = {'timeHgValues':timeHgValues, 'pExponFit':pExponGuess, 
                  'pGaussFit':pGaussGuess, 'xLimit':xLimit, 
                  'cGaussGuess':cGaussGuess, 'timeStamps':timeStamps}

        return retval

    @staticmethod
    def intervalTime(intervals):
        """
        return the time (in seconds) masked by the intervals
        """
        time = 0
        for interval in intervals:
            time += interval[1]-interval[0]
        return time

    @staticmethod
    def funcExpon(x, a, b, c, d):
        retval = a*np.exp(-b*(x-d)) + c
        retval[x < d] = 0
        return retval

    @staticmethod
    def funcGauss(x, a, b, c):
        return a*np.exp(-(x-b)**2/(2.*c**2))
       
