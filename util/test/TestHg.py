import unittest
import numpy as np
import time
import matplotlib.pyplot as plt
from util import hgPlot
import inspect
class TestHg(unittest.TestCase):
    """
    Test options for filling and plotting histograms
    """
    def testCompare(self):
        """
        Demonstrate that np.bincounts is 40 times faster than np.histogram
        """
        nPixels = 100
        values = []
        for i in range(nPixels):
            values.append(np.random.random_integers(0,1000000, 100))
            values[i].sort()

        # measure time for np.histogram
        begin = time.time()        
        for i in range(nPixels):
            hg = np.histogram(values[i], 1000000, range=(0,1000000)) # 3.8 sec
        end = time.time()
        deltaHg = end - begin
        #print "np.histogram elapsed time is",deltaHg

        # measure time for np.bincount
        begin = time.time()        
        for i in range(nPixels):
            hg = np.bincount(values[i],minlength=1000000) # 0.097 sec/100
        end = time.time()
        deltaBc = end - begin
        # print "np.bincount elapsed time is",deltaBc
        self.assertTrue(deltaBc*20 < deltaHg)

    def testHgPlot1(self):
        xmax = 5
        hg = np.histogram([0,1,1,1,1,1,2,2],bins=xmax, range=(-0.5,xmax-0.5))
        x,y = hgPlot.getPlotValues(hg, ylog=False)
        plt.clf()
        plt.plot(x,y)
        plt.margins(0.1, 0.1)
        tfn = inspect.stack()[0][3]
        plt.savefig(tfn)

    def testHgPlot2(self):
        hg = [0,0,1,8,4,0]
        x,y = hgPlot.getPlotValues(hg, ylog=False)
        plt.clf()
        plt.plot(x,y)
        plt.margins(0.1, 0.1)
        tfn = inspect.stack()[0][3]
        plt.savefig(tfn)

    def testHgPlot3(self):
        hg = [10000,0.1,0,1]
        x,y = hgPlot.getPlotValues(hg, ylog=False)
        plt.clf()
        plt.plot(x,y)
        plt.margins(0.1, 0.1)
        plt.yscale('log')
        tfn = inspect.stack()[0][3]
        plt.savefig(tfn)


if __name__ == '__main__':
    unittest.main()
