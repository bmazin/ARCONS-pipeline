from util import utils
import matplotlib as mpl
import numpy as np
import unittest
import os
class TestUtils(unittest.TestCase):
    """
    Test functions in utils.py
    """
    def testMakeMovie0(self):
        """
        make a simple movie all 0 and 1
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
        utils.makeMovie(listOfFrameObj, frameTitles, outName='movie0')
    def testMakeMovie1(self):
        """
        make a simple movie all 0 and 1 and values in between
        """
        nrow = 5
        ncol = 10
        listOfFrameObj = []
        frameTitles = []
        listOfPixelsToMark = []
        for iFrame in range(nrow):
            print "iFrame=",iFrame
            frame = []
            for iRow in range(nrow):
                row = [float(iRow)/nrow]*ncol
                if (iRow <= iFrame):
                    row[ncol/2] = 1
                frame.append(row)
            listOfFrameObj.append(np.array(frame))
            frameTitles.append("frame=%03d"%iFrame)
            listOfPixelsToMark.append([(iFrame,iFrame)])
        utils.makeMovie(listOfFrameObj, frameTitles, outName='movie1',
                        delay=1.0, colormap=mpl.cm.gray,
                        listOfPixelsToMark=listOfPixelsToMark, 
                        pixelMarkColor='red')
if __name__ == '__main__':
    unittest.main()
