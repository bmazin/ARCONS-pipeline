from util import utils
import matplotlib as mpl
import numpy as np
import unittest
import os
import inspect
from util.readDict import readDict
import math
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
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

    def testFitRigidRotation(self):
        """
        calculate ra,dec from a grid of x,y points

        use these mathes to calculate the WCS transformation

        test that the WCS transform works to 1e-9 degrees
        """
        nx = 3
        ny = 3
        xs = np.zeros(nx*ny,dtype=np.float)
        ys = np.zeros(nx*ny,dtype=np.float)
        ipt = 0
        for ix in range(nx):
            for iy in range(ny):
                xs[ipt] = ix*15 + 10
                ys[ipt] = iy*15 + 20
                ipt += 1
        for scale in [0.1/3600, 0.23/3600]: # degrees per pixel
            for theta in [10,123.4, -88, 196]: # rotation in degree
                ct = math.cos(math.radians(theta))
                st = math.sin(math.radians(theta))
                dra =  4.0/3600
                ddec = 5.0/3600
                ras  = scale*(xs*ct - ys*st) + dra
                decs = scale*(xs*st + ys*ct) + ddec
                w = utils.fitRigidRotation(xs,ys,ras,decs)

                wras,wdecs = w.wcs_pix2world(xs,ys,1)
                for x,y,wra,wdec,ra,dec in zip(xs,ys,wras,wdecs,ras,decs):
                    cTrue = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
                    cReco = SkyCoord(ra=wra*u.degree, dec=wdec*u.degree)
                    sep = cTrue.separation(cReco)
                    dra = (ra-wra)*3600
                    ddec = (dec-wdec)*3600
                    self.assertTrue(sep.degree < 1e-9 )
if __name__ == '__main__':
    unittest.main()
