import os
import math
import unittest
from cosmic.Cosmic import Cosmic
import matplotlib.pyplot as plt
from util import FileName
from interval import interval, inf, imath
import inspect
import numpy as np
from util import hgPlot
from scipy.stats import poisson
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

        plt.clf()
        inter = interval()
        for stride in [1,10]:
            timeHgValues, populationHg = \
                self.cosmic.getHgsForOneSec(\
                self.cosmic.beginTime,inter,stride=stride,populationMax=20)
            xValues = populationHg[1][1:]-0.5
            xp,yp = hgPlot.getPlotValues\
                (populationHg, ylog=True)
            ypNorm = np.array(yp,dtype=np.double)/sum(populationHg[0])
            
            nEntries = sum(timeHgValues)
            nBins = timeHgValues.size
            mean = stride*nEntries/float(nBins)
            plt.plot(xp,ypNorm, label="data stride=%d"%stride)
            plt.yscale('log')
            limits = plt.axis();
            poisson = []
            for x in xValues:
                prob = (mean**x)*math.exp(-mean)/math.factorial(x)
                poisson.append(prob)
                pError = np.sqrt(poisson/sum(populationHg[0]))
            plt.errorbar(xValues, poisson, pError, \
                                 label="poisson stride=%d"%stride)
            plt.ylim(1e-8, 2)
                
            plt.margins(0.1, 0.1)
            ylim = plt.ylim()
            plt.legend()

            
        title = "%s %s %s sec=%d"%(self.fn.run,self.fn.date,self.fn.tstamp,self.cosmic.beginTime)
        plt.title(title)
        plt.xlabel("N")
        plt.ylabel("probability(N)")
        plt.savefig(self.fn.makeName(inspect.stack()[0][3]+"_",""))

    def test_populationFromTimeHgValues(self):
        minLength = 1000000
        timeHgValues = np.zeros(minLength, dtype=np.int64)
        timeHgValues[90001] = 1
        timeHgValues[90002] = 1
        timeHgValues[90003] = 1

        populationMax = 5
        stride = 10
        populationHg = Cosmic.populationFromTimeHgValues(
            timeHgValues, populationMax, stride)
        #print "populationHg[0]=",populationHg[0]
        self.assertEquals(populationHg[0][0], 199997)
        self.assertEquals(populationHg[0][1], 0)
        self.assertEquals(populationHg[0][2], 0)
        self.assertEquals(populationHg[0][3], 2)
        self.assertEquals(populationHg[0][4], 0)

    def test_FindCosmics(self):
        self.cosmic.findFlashes()
        plt.clf()
        inter = interval()
        stride = 10
        threshold = 10
        populationMax = 100
        populationHg = self.cosmic.findCosmics(stride, threshold, populationMax)
        xValues = populationHg[1][1:]-0.5
        xp,yp = hgPlot.getPlotValues\
            (populationHg, ylog=True)
        ypNorm = np.array(yp,dtype=np.double)/sum(populationHg[0])
        
        plt.plot(xp, yp, label="data stride=%d"%stride)

        #nEntries = sum(timeHgValues)
        #nBins = timeHgValues.size
        #mean = stride*nEntries/float(nBins)

        plt.yscale('log')
        plt.savefig(self.fn.makeName(inspect.stack()[0][3]+"_",""))

if __name__ == '__main__':
    unittest.main()
