import os
import unittest
from cosmic.Cosmic import Cosmic
import matplotlib.pyplot as plt
from util import FileName

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
        self.cosmic = Cosmic(self.fn, beginTime= 223, endTime=227,\
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
        self.assertEqual(self.cosmic.nPhoton(), 815755)

    def test_makeTimeHgs(self):
        """make and plot the time hgs"""
        self.cosmic.findFlashes()
        self.cosmic.writeFlashesToHdf5()
        #self.cosmic.plotTimeHgs()
        #plt.savefig(self.fileName+".png")


if __name__ == '__main__':
    unittest.main()
