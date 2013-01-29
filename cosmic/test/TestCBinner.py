import os
import unittest
import numpy as np
from cosmic import c_binner
from cosmic import tsBinner
class TestCBinner(unittest.TestCase):

    def testCBinner(self):
        
        ts = np.array([1,3,5,6], dtype=np.uint64)
        bins = np.zeros(10, dtype=np.int64)
        c_binner.binner(ts,bins)
        print "bins=",bins

    def testTsBinner(self):
        ts = np.array([1,3,5,6], dtype=np.uint64)
        bins = np.zeros(10, dtype=np.int)
        c_binner.binner(ts,bins)
        print "bins=",bins

if __name__ == '__main__':
    unittest.main()
