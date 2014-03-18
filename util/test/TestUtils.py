from util import utils
import matplotlib as mpl
import numpy as np
import unittest
import os
import inspect
from util.readDict import readDict

class TestUtils(unittest.TestCase):
    """
    Test functions in utils.py
    """
    def testPlotArray(self):
        "exercise the plotArray function and make the file testPlotArray.png"
        xyarray = np.arange(20).reshape((4,5)) - 5
        fn1 = inspect.stack()[0][3]+".png"
        utils.plotArray(xyarray, showMe=False, cbar=True,
                        cbarticks=[-4, 1,2,4,8,16],
                        cbarlabels=['negative four', 'one','two','four','eight','sixteen'],
                        plotTitle='This is the Plot Title!',
                        colormap=mpl.cm.terrain,
                        pixelsToMark=[(0,1)],
                        pixelMarkColor='red',
                        plotFileName=fn1,
                        sigma=2.0)

    def testMakeMovie0(self):
        """
        make a simple movie all 0 and 1
        """
        nrow = 20
        ncol = 10
        listOfFrameObj = []
        frameTitles = []
        for iFrame in range(nrow):
            #print "iFrame=",iFrame
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
    def testGetGitStatus(self):
        """
        Test that the getStatus method does not crash and returns 
        a dictionary with more than 2 keys
        """
        gs = utils.getGitStatus()
        self.assertTrue(len(gs.keys()) > 2)

    def testReadDict(self):
        """
        Test reading geminga.dict
        """
        params = readDict()
        params.read_from_file("geminga.dict")
        self.assertTrue(params['needHotPix'])
        self.assertEqual(26, len(params['obsSequence']))

if __name__ == '__main__':
    unittest.main()
