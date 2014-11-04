import unittest
import os
from util import FileName
from util import ObsFile

class TestFileName(unittest.TestCase):
    """
    Test getting the raw and the timeMask file names

    Works by default on turk; set environment variables 
    MKID_RAW_PATH and MKID_PROC_PATH on other machines

    Check for one file and see that it exists:
    /ScienceData/LICK2012/20120919/obs_20120920-123350.h5
    """
    def testRaw(self):
        mkidDataDir = os.getenv('MKID_RAW_PATH', default="/ScienceData")
        intermDir = os.getenv('MKID_PROC_PATH', default="/Scratch")
        run = 'LICK2012'
        date = '20120919'
        tstamp = '20120920-123350'
        fn = FileName.FileName(run, date, tstamp, mkidDataDir, intermDir)

        obsFn = fn.obs()
        self.assertTrue(os.path.exists(obsFn), msg=obsFn+" does not exist")

        calFn = fn.cal()
        self.assertFalse(os.path.exists(calFn))


    def testCalSoln(self):
        mkidDataDir = os.getenv('MKID_RAW_PATH', default="/ScienceData")
        intermDir = os.getenv('MKID_PROC_PATH', default="/Scratch")
        run = 'LICK2012'
        date = '20120919'
        tstamp = '20120920-123900'
        fn = FileName.FileName(run, date, tstamp, mkidDataDir, intermDir)
        calSolnFn = fn.calSoln()
        self.assertTrue(os.path.exists(calSolnFn), msg=calSolnFn)

    def testTimeMask(self):
        mkidDataDir = os.getenv('MKID_RAW_PATH', default="/ScienceData")
        intermDir = os.getenv('MKID_PROC_PATH', default="/Scratch")
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
        fn1 = FileName.FileName(run, date, tstamp)
        obsFile = ObsFile.ObsFile(fn1.obs())
        fn2 = FileName.FileName(obsFile=obsFile)
        self.assertTrue(fn1.obs()==fn2.obs())
        self.assertTrue(fn2.getComponents()==(run,date,tstamp))
        
    def testFlatCalFromObsFile(self):
        run = 'LICK2012'
        date = '20120919'
        tstamp = '20120920-123350'
        fn = FileName.FileName(run, date, tstamp)
        obsFile = ObsFile.ObsFile(fn.obs())
        flatCalFN1 = FileName.FileName(run, date, tstamp='').flatSoln()
        flatCalFN2 = FileName.FileName(obsFile=obsFile).flatSoln()        
        self.assertTrue(flatCalFN1==flatCalFN2)
        
    def testPassObsFileNameString(self):
        run = 'LICK2012'
        date = '20120919'
        tstamp = '20120920-123350'
        fn = FileName.FileName(run, date, tstamp)
        fullFileNameStr = fn.obs()
        self.assertTrue(FileName.FileName(obsFile=fullFileNameStr).obs() == fullFileNameStr)
        self.assertTrue(FileName.FileName(obsFile=fullFileNameStr).timeMask() == fn.timeMask())
        
    def testPixRemapFile(self):
        run = 'PAL2012'
        fn = FileName.FileName(run)
        pixRemapFn = fn.pixRemap()
        self.assertTrue(os.path.exists(pixRemapFn), msg=pixRemapFn)    
        
if __name__ == '__main__':
    unittest.main()
