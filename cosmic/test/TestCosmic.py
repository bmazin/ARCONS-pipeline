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
from scipy.stats import poisson,expon
from cosmic import tsBinner
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
        try:
            self.fn = FileName.FileName('LICK2012','20120919',  '20120920-092626')
            self.cosmic = Cosmic(self.fn, beginTime= 222, endTime=228,\
                                 nBinsPerSec = 10)
            self.assertEqual(self.cosmic.file.beamImage[1][2], 
                             '/r7/p162/t1348133188')
        except:
#            pass
            self.fail("Trouble in TestCosmic.setUp for file="+self.fn.obs())

    def testBadEndTime(self):
        try:
            cosmic = Cosmic(self.fn, 0, 3000)
        except RuntimeError, e:
            self.assertEquals("bad endTime:  endTime=3000 exptime=298", \
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

    def test_populationFromTimeHgValues(self):
        """
        test the funcion Cosmic.populationFromTimeHgValues.

        Make a time stream where in successive time values you have 1,3,1
        photons.  Look for cosmic events with a threshold of 4.  Arrange
        the time stream so that one of the large bins (starting at 89995)
        overlaps the time stamps with all five photons, but the next one
        (starting at 90000) only overlaps 4 photons.  With the threshold
        set to 4, the bin starting at 89995 is found and the bin starting
        at 90000 is not found.
        """
        minLength = 1000000
        timeHgValues = np.zeros(minLength, dtype=np.int64)
        timeHgValues[89999] = 1
        timeHgValues[90000] = 3
        timeHgValues[90001] = 1

        populationMax = 10
        stride = 10
        threshold = 4
        pfthv = Cosmic.populationFromTimeHgValues(
            timeHgValues, populationMax, stride, threshold)

        cosmicTimeList = pfthv['cosmicTimeList']
        #print "cosmicTimeList=",cosmicTimeList
        self.assertEquals(cosmicTimeList.size,1)
        self.assertEquals(cosmicTimeList[0], 89995)
        populationHg = pfthv['populationHg']
        #print "populationHg[0]=",populationHg[0]
        self.assertEquals(populationHg[0][0], 199996)
        self.assertEquals(populationHg[0][1], 1)
        self.assertEquals(populationHg[0][2], 0)
        self.assertEquals(populationHg[0][3], 0)
        self.assertEquals(populationHg[0][4], 1)
        self.assertEquals(populationHg[0][5], 1)


    def test_FindCosmics(self):
        self.cosmic.findFlashes()
        plt.clf()
        inter = interval()
        stride = 10
        threshold = 10
        populationMax = 100
        fc = self.cosmic.findCosmics(stride, threshold, populationMax)
        cosmicTimeList = fc['cosmicTimeList']
        print "number of cosmicTimeLists=",len(cosmicTimeList)
        # Value 23294 found on May 6, 2013
        self.assertEquals(23294,len(cosmicTimeList))
        #populationHg = fc['populationHg']
        #xValues = populationHg[1][1:]-0.5
        #xp,yp = hgPlot.getPlotValues\
        #    (populationHg, ylog=True)
        #ypNorm = np.array(yp,dtype=np.double)/sum(populationHg[0])
        
        #plt.plot(xp, yp, label="data stride=%d"%stride)

        #nEntries = sum(timeHgValues)
        #nBins = timeHgValues.size
        #mean = stride*nEntries/float(nBins)

        plt.yscale('log')
        plt.savefig(self.fn.makeName(inspect.stack()[0][3]+"_",""))

    def testExpon(self):
        """
        generate and fit an exponential distribution with lifetime of 25
        make a plot in testExpon.png
        """
        tau = 25.0
        nBins = 400
        size = 100
        x = range(nBins)
        timeHgValues = np.zeros(nBins, dtype=np.int64)
        timeStamps = expon.rvs(loc=0, scale=tau, size=size)
        ts64 = timeStamps.astype(np.uint64)
        tsBinner.tsBinner(ts64, timeHgValues)

        param = expon.fit(timeStamps)
        fit = expon.pdf(x,loc=param[0],scale=param[1])
        fit *= size
        tvf = timeHgValues.astype(np.double)
        tvf[tvf<1] = 1e-3 # the plot looks nicer if zero values are replaced
        plt.plot(x, tvf, label="data")
        plt.plot(x, fit, label="fit")
        plt.yscale('log')
        plt.xlim(xmax=100)
        plt.ylim(ymin=0.09)
        plt.legend()
        plt.title("true tau=%.1f   fit tau=%.1f"%(tau,param[1]))
        plt.savefig(inspect.stack()[0][3]+".png")

    def testFitDecayTime(self):
        run = 'PAL2012'
        sundownDate = '20121211'
        obsDate = '20121212'
        seq = '133303'
        t0 = 138595580
        fn = FileName.FileName(run, sundownDate, obsDate+"-"+seq)
        cosmic = Cosmic(fn, endTime='exptime')
        cosmic.fitDecayTime()

if __name__ == '__main__':
    unittest.main()
