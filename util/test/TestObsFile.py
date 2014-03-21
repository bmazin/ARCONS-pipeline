import os
import unittest
import numpy as np
from interval import interval, inf, imath
import matplotlib.pyplot as plt
import inspect
from util import ObsFile, FileName
from cosmic.Cosmic import Cosmic
import time
import cProfile
import re
import StringIO, pstats
#import interval
#from util.FileName import FileName
#from util.ObsFile import ObsFile

class TestObsFile(unittest.TestCase):
    def testNegativeTimeIssue(self):
        """
        This particular sequence has a photon close to time=0.  Make sure
        it returns a positive value after loadTimeAdjustmentFile fixes
        everything else.
        """
        run = 'PAL2012'
        sundownDate = '20121211'
        obsDate = '20121212'
        seq = '112709'
        fileName = FileName.FileName(run, sundownDate, obsDate+"-"+seq)
        timeAdjustments = fileName.timeAdjustments()
        obsFile = ObsFile.ObsFile(fileName.obs())
        obsFile.loadTimeAdjustmentFile(timeAdjustments)
        iRow = 9
        iCol = 22
        firstSec = 0
        integrationTime = -1
        gtpl = obsFile.getTimedPacketList(iRow, iCol, 
                                          firstSec, integrationTime)
        ts = gtpl['timestamps']
        self.assertTrue(ts[0] >= 0)

    def testGetTimedPacketList(self):
        """
        Sample call for ObsFile method getTimedPacketList

        Checks that some values are the same as they were on April 8
        """
        fn = FileName.FileName('LICK2012','20120919',  '20120920-092626')
        obsFile = ObsFile.ObsFile(fn.obs())
        exptime0 = obsFile.getFromHeader('exptime')
        self.assertEquals(exptime0, 300)
        timeAdjustments = fn.timeAdjustments()
        obsFile.loadTimeAdjustmentFile(timeAdjustments)
        exptime1 = obsFile.getFromHeader('exptime')
        self.assertEquals(exptime1, 298)

        iRow = 30
        iCol = 32
        
        tpl = obsFile.getTimedPacketList(iRow,iCol)
        timestamps = tpl['timestamps']

        # Chris S. found these values on April 8, 2012
        self.assertEquals(timestamps.size,145542)
        self.assertAlmostEquals(timestamps[0],0.0028879999999,places=10)
        self.assertAlmostEquals(timestamps[-1],297.9997779999999, places=10)

        timestampsUint64 = (timestamps*1e6).astype(np.uint64)
        self.assertEquals(timestampsUint64.size,145542)
        self.assertEquals(timestampsUint64[0],2888)
        self.assertEquals(timestampsUint64[-1],297999778)

        del obsFile

    def testMakeMask(self):
        """
        tests that makeMask gets the same results as the __contains__
        method of interval
        """
        timestamps = np.array([1.1, 2.2, 3.3, 4.4, 8.8])
        inter = interval([1.9,2.5],[4.3,5.3])
        mask = ObsFile.ObsFile.makeMask(timestamps,inter)
        for i in range(len(timestamps)):
            self.assertTrue(mask[i] == inter.__contains__(timestamps[i]))

    def testGetPackets(self):
        """
        test the getPackets method, which could replace getTimedPacketList.
        by masking more efficiently.  Report hthe number of packets
        returned and the time for the two methods
        """

        fn = FileName.FileName('LICK2012','20120919',  '20120920-092626')
        obsFile = ObsFile.ObsFile(fn.obs())
        exptime0 = obsFile.getFromHeader('exptime')
        self.assertEquals(exptime0, 300)
        timeAdjustments = fn.timeAdjustments()
        obsFile.loadTimeAdjustmentFile(timeAdjustments)
        exptime1 = obsFile.getFromHeader('exptime')
        self.assertEquals(exptime1, 298)


        exptime1 = 4
        iRow = 30
        iCol = 32

        fields = ['peakHeights']

        # mask out from 0.5 to 0.75 each second
        inter = interval()
        for sec in range(exptime1):
            inter = inter | interval([sec+0.5, sec+0.75])
        obsFile.cosmicMask = inter
        obsFile.switchOnCosmicTimeMask()

        print "now call getPackets"
        start = time.clock()
        tpl = obsFile.getPackets(iRow, iCol, 0, exptime1, fields=fields)
        elapsed = (time.clock() - start)
        print "getPackets elapsed=",elapsed

        print "now call getTimedPacketList"
        start = time.clock()
        gtpl = obsFile.getTimedPacketList(iRow, iCol, 0, exptime1)
        elapsed = (time.clock() - start)
        print "getTimedPackeList elapsed=",elapsed

        print "compare"
        print "len(tpl['timestamps'])=",len(tpl['timestamps'])
        print "len(gtpl['timestamps'])=",len(gtpl['timestamps'])
        print "now make plot"
        plt.clf()
        plt.plot(tpl['timestamps'], tpl['peakHeights'])
        plt.savefig(inspect.stack()[0][3]+".png")

    def testGetPacketsProfile(self):
        fn = FileName.FileName('LICK2012','20120919',  '20120920-092626')
        obsFile = ObsFile.ObsFile(fn.obs())
        exptime0 = obsFile.getFromHeader('exptime')
        self.assertEquals(exptime0, 300)
        timeAdjustments = fn.timeAdjustments()
        obsFile.loadTimeAdjustmentFile(timeAdjustments)
        exptime1 = obsFile.getFromHeader('exptime')
        self.assertEquals(exptime1, 298)

        exptime1 = 10
        iRow = 30
        iCol = 32

        fields = ['peakHeights']

        # mask out from 0.5 to 0.75 each second
        inter = interval()
        for sec in range(exptime1):
            inter = inter | interval([sec+0.1, sec+0.11])
            inter = inter | interval([sec+0.2, sec+0.21])
            inter = inter | interval([sec+0.3, sec+0.31])
            inter = inter | interval([sec+0.4, sec+0.41])
            inter = inter | interval([sec+0.5, sec+0.51])
            inter = inter | interval([sec+0.6, sec+0.61])
            inter = inter | interval([sec+0.7, sec+0.71])
            inter = inter | interval([sec+0.8, sec+0.81])
            inter = inter | interval([sec+0.9, sec+0.91])
        obsFile.cosmicMask = inter
        obsFile.switchOnCosmicTimeMask()


        print "now call getPackets"
        start = time.clock()
        pr = cProfile.Profile()
        pr.enable()
        start = time.clock()
        tpl = obsFile.getPackets(iRow, iCol, 0, exptime1, fields=fields)
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        #print s.getvalue()
        elapsed = (time.clock() - start)
        print "getPackets elapsed=",elapsed


    def testGetPacketsAndGetTimedPacketList(self):
        fn = FileName.FileName('LICK2012','20120919',  '20120920-092626')
        obsFile = ObsFile.ObsFile(fn.obs())
        exptime0 = obsFile.getFromHeader('exptime')
        self.assertEquals(exptime0, 300)
        timeAdjustments = fn.timeAdjustments()
        obsFile.loadTimeAdjustmentFile(timeAdjustments)
        exptime1 = obsFile.getFromHeader('exptime')
        self.assertEquals(exptime1, 298)

        exptime1 = 100
        iRow = 30
        iCol = 32

        fields = ['peakHeights', 'baselines']

        # mask out from 0.5 to 0.75 each second
        inter = interval()
        for sec in range(exptime1):
            inter = inter | interval([sec+0.1, sec+0.11])
            inter = inter | interval([sec+0.2, sec+0.21])
            inter = inter | interval([sec+0.3, sec+0.31])
            inter = inter | interval([sec+0.4, sec+0.41])
            inter = inter | interval([sec+0.5, sec+0.51])
            inter = inter | interval([sec+0.6, sec+0.61])
            inter = inter | interval([sec+0.7, sec+0.71])
            inter = inter | interval([sec+0.8, sec+0.81])
            inter = inter | interval([sec+0.9, sec+0.91])
        obsFile.cosmicMask = inter
        obsFile.switchOnCosmicTimeMask()

        print "now call getTimedPacketList"
        start = time.clock()
        gtpl = obsFile.getTimedPacketList(iRow, iCol, 0, exptime1)
        elapsed = (time.clock() - start)
        print "getTimedPackeList elapsed=",elapsed

        print "now call getPackets"
        start = time.clock()
        tpl = obsFile.getPackets(iRow, iCol, 0, exptime1, fields=fields)
        elapsed = (time.clock() - start)
        print "getPackets elapsed=",elapsed


        print "compare"
        print "tpl keys=",tpl.keys()
        print "gtpl keys=",gtpl.keys()
        print "len(tpl['timestamps'])=",len(tpl['timestamps'])
        print "len(gtpl['timestamps'])=",len(gtpl['timestamps'])
        self.assertTrue((tpl['timestamps'] == gtpl['timestamps']).all())
        self.assertTrue((tpl['baselines'] == gtpl['baselines']).all())
    def testGetFrame(self):
        """
        Demonstrate the getFrame method of ObsFile and check a specific result
        """
        fn = FileName.FileName('LICK2012','20120919',  '20120920-092626')
        obsFile = ObsFile.ObsFile(fn.obs())
        frame = obsFile.getFrame(0,3)
        shape = frame.shape
        self.assertEquals(obsFile.nRow, shape[0])
        self.assertEquals(obsFile.nCol, shape[1])
        self.assertEquals(frame.sum(), 331994)

    def testCalculateSlicesMiddle(self):
        "one interval in the middle of the set of timestamps"
        secs = np.arange(10)
        inter = interval([2,6])
        slices = ObsFile.calculateSlices(inter, secs)
        self.assertEquals(slices, ["0:2","7:10"])
        new = ObsFile.repackArray(secs, slices)
        self.assertTrue(npEquals(new,np.array([0,1,7,8,9])))

    def testCalculateSlicesBeginningAndEnd(self):
        "one interval at the beginning and one at the end"
        secs = np.arange(10)
        inter = interval([0,3]) | interval([8,9])
        slices = ObsFile.calculateSlices(inter, secs)
        new = ObsFile.repackArray(secs, slices)
        self.assertTrue(npEquals(new, np.array([4,5,6,7])))

    def testParsePhotonPacketsWithIntervalMask(self):
        """
        parse one second of data with and without masking a range of times.
        plot the result
        """
        fn = FileName.FileName('LICK2012','20120919',  '20120920-092626')
        fullFileName = fn.obs()
        obsFile = ObsFile.ObsFile(fn.obs())
        pixel = obsFile.getPixel(3,4)
        packet = pixel[123]
        ts,pfp,bl = obsFile.parsePhotonPackets(packet)
        # specific to the file, pixel, and second we parsed
        self.assertEquals(len(ts),66)
        xs = [200000,400000]
        inter = interval([xs[0],xs[1]])
        tsm,pfpm,blm = obsFile.parsePhotonPackets(packet, inter)
        # specific to the file, pixel, and second we parsed
        self.assertEquals(len(tsm),50)
        plt.plot(ts,pfp,label="all points");
        plt.plot(tsm,pfpm,'ro',label="not masked points");
        ys = plt.axis()
        plt.plot([xs[0],xs[0]],[ys[2],ys[3]],'green')
        plt.plot([xs[1],xs[1]],[ys[2],ys[3]],'green')
        plt.legend(loc="lower right")
        plt.xlabel("Time (ticks)")
        plt.ylabel("peak flux parabola fit")
        plt.title("mask times from %d to %d ticks" % (xs[0],xs[1]))
        plt.savefig(fn.makeName(inspect.stack()[0][3]+"_",""))

    
    def testNumpyPytablesConflict(self):
        '''
        Author: Julian van Eyken                Date: June 17 2013
        
        Test for Numpy/PyTables version conflict: it seems a bug shows
        up when using Numpy 1.7.1 with PyTables 2.3.1 and 2.4 (at least -
        and maybe other versions of PyTables too?). Most likely it
        appears to be related to a change to the way Numpy handles
        integer overflows in v. 1.7.1. It appears to cause occasional
        random failures when calling the pytables getNode() function.
        This test just tries to recreate the error as a check to make
        sure the conflict is not happening.
        
        Note that this is not 100% guaranteed to pick up the bug
        every time, since it's not deterministic.... But it seems
        to usually do the trick.
        '''
        
        fn = FileName.FileName(run='PAL2012',date='20121211',
                               tstamp='20121212-033323').obs()
        obsFile = ObsFile.ObsFile(fn)
        #Just try to get an empty image - that should call
        #getNode() on every pixel....
        
        print 'Getting a test image...' 
        try:
            dummyImage = obsFile.getPixelCountImage(firstSec=0,integrationTime=-1,
                                                getRawCount=True)
        except AttributeError:
            print '!!! Looks like the Numpy/PyTables bug probably showed up.'
            print '!!! Probably because you''re running Numpy 1.7.1 with PyTables...?'
            raise
        
        print 'Done. Seemed to work okay that time....'

    def testCosmicTimeMask(self):
        """

Test applying time masks to one file.
Note that this file has some hot pixels detected for some seconds,
for example, this one:
iRow= 44  iCol= 43
   interval= interval([0.0, 29.0])
   interval= interval([59.0, 61.0])
   interval= interval([62.0, 220.0])
   interval= interval([225.0, 283.0])
   interval= interval([286.0, 300.0])
"""

        fn = FileName.FileName(run='PAL2012',date='20121211',
                               tstamp='20121212-033323').obs()
        obsFile = ObsFile.ObsFile(fn)

        fileName = FileName.FileName(obsFile=fn)
        timeAdjustments = fileName.timeAdjustments()

        obsFile.loadTimeAdjustmentFile(timeAdjustments)

        iRow = 44
        iCol = 43
        firstSec = 58
        integrationTime = 5
        obsFile.switchOffHotPixTimeMask()
        gtpl0 = obsFile.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
        times0 = gtpl0['timestamps']
        peaks0 = gtpl0['peakHeights']

        hotPixCalFileName = fileName.timeMask()
        obsFile.loadHotPixCalFile(hotPixCalFileName)
        obsFile.switchOnHotPixTimeMask()

        gtpl1 = obsFile.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
        times1 = gtpl1['timestamps']
        peaks1 = gtpl1['peakHeights']

        cosmicInterval = interval()
        cosmicInterval = cosmicInterval | interval([58.25,58.75])

        beginTime = 1.234
        endTime = 5.432
        stride = 12
        threshold = 6543
        nSigma = 7.8
        populationMax = 2222
        ObsFile.ObsFile.writeCosmicIntervalToFile(cosmicInterval, 
                                                  obsFile.ticksPerSec, 
                                                  'temp.h5',
                                                  beginTime, endTime, stride,
                                                  threshold, nSigma, populationMax)

        readInterval = ObsFile.ObsFile.readCosmicIntervalFromFile("temp.h5")

        obsFile.loadCosmicMask('temp.h5')
        gtpl2 = obsFile.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
        times2 = gtpl2['timestamps']
        peaks2 = gtpl2['peakHeights']

        plt.clf()
        plt.subplot(311)
        plt.scatter(times0, peaks0, marker='o', s=1)
        plt.ylabel("no mask")
        plt.xlim(firstSec, firstSec+integrationTime)

        plt.subplot(312)
        plt.scatter(times1, peaks1, marker='o', s=1)
        plt.ylabel("hot pixel mask")
        plt.xlim(firstSec, firstSec+integrationTime)

        plt.subplot(313)
        plt.scatter(times2, peaks2, marker='o', s=1) 
        plt.ylabel("hot pix + cosmic")
        plt.xlim(firstSec, firstSec+integrationTime)

        plt.savefig(inspect.stack()[0][3]+".png")

    def testGetTimedPacketList2(self):
        """
        this seems to get a timestamp outside of beginTime, beginTime+expTime.  Matt knows about this, and since it happens only for the first photon it is not a priority.
        """
        run = 'PAL2012'
        sundownDate = '20121211'
        obsDate = '20121212'
        seq = '115722'
        fn = FileName.FileName(run, sundownDate,obsDate+"-"+seq)
        print "now open",fn.obs()
        obsFile = ObsFile.ObsFile(fn.obs())
        obsFile.loadTimeAdjustmentFile(fn.timeAdjustments())
        beginTime = 0
        integrationTime = 0.177653
        iRow = 3
        iCol = 1
        print "beginTime=",beginTime, "iRow=",iRow,"iCol=",iCol,"integrationTime=",integrationTime
        endTime = beginTime + integrationTime
        gtpl = obsFile.getTimedPacketList(iRow,iCol,beginTime,integrationTime)
        print "timestamps=",gtpl['timestamps']
        for timestamp in gtpl['timestamps']:
            msg = "beginTime=%f   timestamp=%f  endTime=%f"%(beginTime,timestamp,endTime)
            self.assertTrue(timestamp >= beginTime, msg)
            self.assertTrue(timestamp <= endTime, msg)

    def testInvertInterval(self):
        """
        test the invertInterval method of ObsFile
        Create in interval, and see that it returns the expected interval
        """
        i = interval()
        i = i | interval[1,2]
        i = i | interval[4,5]
        # all intervals
        inverted = ObsFile.ObsFile.invertInterval(i)
        self.assertEquals(len(inverted),3)
        self.assertEquals(inverted[0][0],float(-inf))
        self.assertEquals(inverted[0][1],1)
        self.assertEquals(inverted[1][0],2)
        self.assertEquals(inverted[1][1],4)
        self.assertEquals(inverted[2][0],5)
        self.assertEquals(inverted[2][1],float(inf))

        # limit min and max
        inverted = ObsFile.ObsFile.invertInterval(i,iMin=0.0, iMax=4.5)
        self.assertEquals(len(inverted),2)
        self.assertEquals(inverted[0][0],0)
        self.assertEquals(inverted[0][1],1)
        self.assertEquals(inverted[1][0],2)
        self.assertEquals(inverted[1][1],4)

def npEquals(a,b):
    if (len(a) != len(b)):
        return False

    for i in range(len(a)):
        if a[i] != b[i]:
            return False

    return True

if __name__ == '__main__':
    unittest.main()
