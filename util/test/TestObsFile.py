import os
import unittest
import numpy as np
from util import ObsFile, FileName
from interval import interval, inf, imath
import matplotlib.pyplot as plt
import inspect
from util.FileName import FileName
from util.ObsFile import ObsFile

class TestObsFile(unittest.TestCase):
    def testNegativeTimeIssue(self):
        run = 'PAL2012'
        sundownDate = '20121211'
        obsDate = '20121212'
        seq = '112709'
        fileName = FileName(run, sundownDate, obsDate+"-"+seq)
        timeAdjustments = fileName.timeAdjustments()
        obsFile = ObsFile(fileName.obs())
        obsFile.loadTimeAdjustmentFile(timeAdjustments)
        iRow = 9
        iCol = 22
        firstSec = 0
        integrationTime = -1
        gtpl = obsFile.getTimedPacketList(iRow, iCol, 
                                          firstSec, integrationTime)
        ts = gtpl['timestamps']
        print "ts[0]=",ts[0]

    def testGetTimedPhotonPacket(self):
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
        self.assertEquals(timestamps[0],0.002847)
        self.assertEquals(timestamps[-1],297.999737)

        timestampsUint64 = (timestamps*1e6).astype(np.uint64)
        self.assertEquals(timestampsUint64.size,145542)
        self.assertEquals(timestampsUint64[0],2847)
        self.assertEquals(timestampsUint64[-1],297999737)

        del obsFile

    def testGetFrame(self):
        fn = FileName.FileName('LICK2012','20120919',  '20120920-092626')
        obsFile = ObsFile.ObsFile(fn.obs())
        frame = obsFile.getFrame(0,-1)
        shape = frame.shape
        self.assertEquals(obsFile.nRow, shape[0])
        self.assertEquals(obsFile.nCol, shape[1])
        print "frame=",frame

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
