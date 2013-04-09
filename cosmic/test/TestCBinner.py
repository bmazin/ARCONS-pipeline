import os
import unittest
import numpy as np
from cosmic import c_binner
from cosmic import tsBinner
import time
class TestCBinner(unittest.TestCase):
    def setUp(self):
        self.answer = [0, 1, 0, 1, 0, 1, 1, 0, 0, 0]
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

    def testTiming(self):
        """
        For binning into nBins=20 bins, for nTs=100,000 time sample values,
        repeat c_binner.binner and then tsBinner.tsBinner nTrials=1000 times.
        Report the elapsed time (in milliseconds) for the two ways.
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
        print "tsBinner=%.2f   cBinner=%.2f"%(dtTsBinner, dtCBinner)

if __name__ == '__main__':
    unittest.main()
