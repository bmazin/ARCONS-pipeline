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
        tstamp = '20120920-123350'
        fn = FileName.FileName(run, date, tstamp, mkidDataDir, intermDir)

        obsFn = fn.obs()
        self.assertTrue(os.path.exists(obsFn), msg=obsFn+" does not exist")

        calFn = fn.cal()
        self.assertFalse(os.path.exists(calFn))


    def testCalSoln(self):
        mkidDataDir = os.getenv('MKID_DATA_DIR', default="/ScienceData")
        intermDir = os.getenv('INTERM_DIR', default="/Scratch")
        run = 'LICK2012'
        date = '20120919'
        tstamp = '20120920-123900'
        fn = FileName.FileName(run, date, tstamp, mkidDataDir, intermDir)
        calSolnFn = fn.calSoln()
        self.assertTrue(os.path.exists(calSolnFn), msg=calSolnFn)

    def testTimeMask(self):
        mkidDataDir = os.getenv('MKID_DATA_DIR', default="/ScienceData")
        intermDir = os.getenv('INTERM_DIR', default="/Scratch")
        run = 'LICK2012'
        date = '20120919'
        tstamp = '20120920-123350'
        fn = FileName.FileName(run, date, tstamp, mkidDataDir, intermDir)
        tmFn = fn.timeMask()
        self.assertGreater(len(tmFn),0)

    def testPassObsFile(self):
        run = 'LICK2012'
        date = '20120919'
        tstamp = '20120920-123350'
        fn = FileName.FileName(run, date, tstamp)
        obsFile = ObsFile(fn1.obs())
        fn2 = FileName.FileName(obsFile)
        self.assertTrue(fn1.obs()==fn2.obs())
        self.assertTrue(fn2.components==(run,date,tstamp))

if __name__ == '__main__':
    unittest.main()
