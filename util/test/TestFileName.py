from util import FileName
import unittest
import os
class TestFileName(unittest.TestCase):
    """
    Test getting the raw and the timeMask file names

    Works by default on turk; set environment variables 
    MKID_DATA_DIR and INTERM_DIR on other machines

    Check for one file and see that it exists:
    /ScienceData/LICK2012/20120919/obs_20120920-123350.h5
    """
    def testRaw(self):
        mkidDataDir = os.getenv('MKID_DATA_DIR', default="/ScienceData")
        intermDir = os.getenv('INTERM_DIR', default="/Scratch")
        run = 'LICK2012'
        date = '20120919'
        fn = FileName.FileName(run, date, mkidDataDir, intermDir)
        flavor = 'obs'
        tstamp = '20120920-123350'
        rawFn = fn.raw(flavor, tstamp)
        #print "rawFn=",rawFn
        self.assertTrue(os.path.exists(rawFn))

    def testCalSoln(self):
        mkidDataDir = os.getenv('MKID_DATA_DIR', default="/ScienceData")
        intermDir = os.getenv('INTERM_DIR', default="/Scratch")
        run = 'LICK2012'
        date = '20120919'
        fn = FileName.FileName(run, date, mkidDataDir, intermDir)
        tstamp = '20120920-123900'
        calFn = fn.calSoln(tstamp)
        print "calFn=",calFn
        self.assertTrue(os.path.exists(calFn))

    def testTimeMask(self):
        mkidDataDir = os.getenv('MKID_DATA_DIR', default="/ScienceData")
        intermDir = os.getenv('INTERM_DIR', default="/Scratch")
        run = 'LICK2012'
        date = '20120919'
        fn = FileName.FileName(run, date, mkidDataDir, intermDir)
        tstamp = '20120920-123350'
        tmFn = fn.timeMask(tstamp)
        print "tmFn=",tmFn
        self.assertGreater(len(tmFn),0)

if __name__ == '__main__':
    unittest.main()
