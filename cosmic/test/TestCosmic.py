import os
import unittest
from cosmic.Cosmic import Cosmic
import matplotlib.pyplot as plt
from util import FileName
from interval import interval, inf, imath
import inspect
import numpy as np
from util import hgPlot
class TestCosmic(unittest.TestCase):
    """
    test the Cosmic class
    from command line:

    To run all tests:
    > python TestCosmic.py

    or to run just one test:
    > python TestCosmic.py TestCosmic.testBadEndTime

    """
    def setUp(self):
        self.fn = FileName.FileName('LICK2012','20120919',  '20120920-092626')
        self.cosmic = Cosmic(self.fn, beginTime= 222, endTime=228,\
                                 nBinsPerSec = 10)
        self.assertEqual(self.cosmic.file.beamImage[1][2], 
                         '/r7/p162/t1348133188')

    def testBadEndTime(self):
        try:
            cosmic = Cosmic(self.fn, 0, 3000)
        except RuntimeError, e:
            self.assertEquals("bad endTime:  endTime=3000 exptime=300", \
                                  e.message)

    def test_nPhoton(self):
        """test the nPhoton method"""
        self.assertEqual(self.cosmic.nPhoton(), 1064969)

    def test_findAndWriteFlashes(self):
        """make and plot the time hgs"""
        self.cosmic.findFlashes()
        self.cosmic.writeFlashesToHdf5()
        #self.cosmic.findCosmics()
        #self.cosmic.plotTimeHgs()
        #plt.savefig(self.fn.makeName("plot_",""))

    def test_getHgForOneSec(self):
        """test and plot results from getHgForOneSec"""

        inter = interval()
        print "now call getHgForOneSec"
        hg, hg2 = self.cosmic.getHgsForOneSec(self.cosmic.beginTime,inter)

        plt.clf()
        xp,yp = hgPlot.getPlotValues(hg2, ylog=True)
        plt.plot(xp,yp)
        plt.margins(0.1, 0.1)
        plt.savefig(self.fn.makeName(inspect.stack()[0][3]+"_",""))
if __name__ == '__main__':
    unittest.main()
