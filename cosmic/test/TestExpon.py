import os
import math
import unittest
from cosmic.Cosmic import Cosmic
import matplotlib.pyplot as plt
from util import FileName
from interval import interval, inf, imath
import inspect
import numpy as np
import pylab as P
from util import hgPlot
from scipy.stats import poisson,expon
from cosmic import tsBinner
class TestExpon(unittest.TestCase):
    
    def test1(self):
        print 'hello'

    def testExponOneEvent(self):
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


    def testExponManyEvents(self):
        """
        generate and fit an exponential distribution with lifetime of 25
        make a plot in testExponManyEvents.png
        """
        tau = 25.0
        nBins = 400
        size = 100
        taulist = []
        for i in range(100):
            x = range(nBins)
            timeHgValues = np.zeros(nBins, dtype=np.int64)
            timeStamps = expon.rvs(loc=0, scale=tau, size=size)
            ts64 = timeStamps.astype(np.uint64)
            tsBinner.tsBinner(ts64, timeHgValues)
            
            param = expon.fit(timeStamps)
            fit = expon.pdf(x,loc=param[0],scale=param[1])
            fit *= size
            #print "i=",i," param=",param
            taulist.append(param[1])
            #mu, sigma = 
            #P.figure()
        n, bins, patches = P.hist(x, 10, normed=1, histtype='step')
        P.setp(patches, 'facecolor', 'g', 'alpha', 1)
        #y = P.normpdf(bins, mu, sigma)
        #line = P.plot(bins, y, 'k--', linewidth=1.5)
            
        #fig  = plt.figure()
        #ax = fig.add_subplot(111)
        #x = taulist
        #numBins = 10
        #ax.hist(x, numBins, color='green', alpha=0.8)
        plt.savefig(inspect.stack()[0][3]+".png")

if __name__ == '__main__':
    unittest.main()
