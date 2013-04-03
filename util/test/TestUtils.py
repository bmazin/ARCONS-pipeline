from util import utils
import numpy as np
import unittest
import os
class TestUtils(unittest.TestCase):
    """
    Test functions in utils.py
    """
    def testMakeMovie(self):
        """
        make a simple movie
        """
        nrow = 20
        ncol = 10
        listOfFrameObj = []
        frameTitles = []
        for iFrame in range(nrow):
            print "iFrame=",iFrame
            frame = []
            for iRow in range(nrow):
                if (iRow < iFrame):
                    row = [1]*ncol
                else:
                    row = [0]*ncol
                frame.append(row)
            listOfFrameObj.append(np.array(frame))
            frameTitles.append("frame=%03d"%iFrame)
        utils.makeMovie(listOfFrameObj, frameTitles)
if __name__ == '__main__':
    unittest.main()
