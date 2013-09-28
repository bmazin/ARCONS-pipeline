import os
import unittest
import numpy as np
from cosmic import c_binner
from cosmic import tsBinner
from cosmic import NPBinner
import time
class TestCBinner(unittest.TestCase):
    def setUp(self):
        self.answer = np.array([0, 1, 0, 1, 0, 1, 1, 0, 0, 0])

    def testNPBinner(self):
        ts = np.array([1,3,5,6], dtype=np.uint32)
        nBins = 10
        contents = NPBinner.NPBinner.binner(ts,nBins)
        self.assertTrue((contents == self.answer).all())
        
    def testCBinner(self):
        """
        confirm that c_binner.binner gets the correct answer 
        for one data set
        """
        ts = np.array([1,3,5,6], dtype=np.uint64)
        bins = np.zeros(10, dtype=np.int64)
        c_binner.binner(ts,bins)
        # print "bins=",bins
        self.assertTrue((bins == self.answer).all())

    def testTsBinner(self):
        """
        confirm that tsBinner.tsBinner gets the correct answer 
        for one data set
        """
        ts = np.array([1,3,5,6], dtype=np.uint64)
        bins = np.zeros(10, dtype=np.int)
        tsBinner.tsBinner(ts,bins)
        # print "bins=",bins
        self.assertTrue((bins == self.answer).all())

    def testTsBinner32(self):
        """
        confirm that tsBinner.tsBinner32 gets the correct answer 
        for one data set
        """
        ts = np.array([1,3,5,6], dtype=np.uint32)
        bins = np.zeros(10, dtype=np.int)
        tsBinner.tsBinner32(ts,bins)
        # print "bins=",bins
        self.assertTrue((bins == self.answer).all())

    def testTiming(self):
        """
        For binning into nBins=20 bins, for nTs=100,000 time sample values,
        repeat c_binner.binner and then tsBinner.tsBinner nTrials=1000 times.
        Then run NPBinner.binner

        Report the elapsed time (in milliseconds) for the two ways.

        Typical speeds are tsbinner=0.17, cbinner=0.53, and NPBinner=1.94

        """
        nBins = 20
        nTs = 100000
        nTrials = 1000
        ts = np.zeros(nTs, dtype=np.uint64)
        ts += np.random.random_integers(0,1,nTs)

        bins = np.zeros(nBins, dtype=np.int)
        t0 = time.time()
        for i in range(nTrials):
            tsBinner.tsBinner(ts,bins)
        t1 = time.time()
        dtTsBinner = t1-t0
        
        bins = np.zeros(nBins, dtype=np.int)
        t0 = time.time()
        for i in range(nTrials):
            c_binner.binner(ts,bins)
        t1 = time.time()
        dtCBinner = t1-t0
        print "tsBinner=%.3f   cBinner=%.3f"%(dtTsBinner, dtCBinner)

        bins = np.zeros(nBins, dtype=np.int)
        ts32 = ts.astype(np.uint32)
        t0 = time.time()
        for i in range(nTrials):
            #print "i=",i,"/",nTrials
            contents = NPBinner.NPBinner.binner(ts32,nBins)
        t1 = time.time()
        dtNPBinner = t1-t0
        print "NPBinner=%.3f"%(dtNPBinner)

        t0 = time.time()
        for i in range(nTrials):
            contents = np.bincount(ts32,minlength=nBins)
        t1 = time.time()
        dtBincount = t1-t0
        print "bincount=%.3f"%(dtBincount)


        t0 = time.time()
        for i in range(nTrials):
            tsBinner.tsBinner32(ts32,bins)
        t1 = time.time()
        dtTsBinner32 = t1-t0
        print "tsBinner32=%.3f"%(dtTsBinner32)

if __name__ == '__main__':
    unittest.main()
