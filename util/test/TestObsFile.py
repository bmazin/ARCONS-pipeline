import os
import unittest
import numpy as np
from interval import interval, inf, imath
import matplotlib.pyplot as plt
import inspect
from util import ObsFile, FileName
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
        inter = interval[2,6]
        slices = ObsFile.calculateSlices(inter, secs)
        self.assertEquals(slices, ["0:2","7:10"])
        new = ObsFile.repackArray(secs, slices)
        self.assertTrue(npEquals(new,np.array([0,1,7,8,9])))

    def testCalculateSlicesBeginningAndEnd(self):
        "one interval at the beginning and one at the end"
        secs = np.arange(10)
        inter = interval[0,3] | interval[8,9]
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
        inter = interval[xs[0],xs[1]]
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

def npEquals(a,b):
    if (len(a) != len(b)):
        return False

    for i in range(len(a)):
        if a[i] != b[i]:
            return False

    return True

if __name__ == '__main__':
    unittest.main()
