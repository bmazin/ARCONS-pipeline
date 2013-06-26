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
        for i in range(1000):
            x = range(nBins)
            timeHgValues = np.zeros(nBins, dtype=np.int64)
            timeStamps = expon.rvs(loc=0, scale=tau, size=size)
            ts64 = timeStamps.astype(np.uint64)
            tsBinner.tsBinner(ts64, timeHgValues)
            
            param = expon.fit(timeStamps)
            fit = expon.pdf(x,loc=param[0],scale=param[1])
            fit *= size
            print "i=",i," param[1]=",param[1]
            taulist.append(param[1]) 

        hist,bins = np.histogram(taulist, bins=20)
        width = 0.7*(bins[1]-bins[0])
        center = (bins[:-1]+bins[1:])/2
        plt.step(center, hist, where = 'mid')
        plt.savefig(inspect.stack()[0][3]+".png")

    def testExponaverage(self):

        """
        generate and fit a histogram of an exponential distribution with many
        events using the average time to find the fit. The histogram is then
        saved in testExponaverage.png
        """
        tau = 25.0
        nBins = 400
        size = 100
        taulist = []
        average = []
        for i in range(1000):
            x = range(nBins)
            timeHgValues = np.zeros(nBins, dtype=np.int64)
            timeStamps = expon.rvs(loc=0, scale=tau, size=size)
            ts64 = timeStamps.astype(np.uint64)
            tsBinner.tsBinner(ts64, timeHgValues)
            
            
            param = sum(timeStamps)/len(timeStamps)
            fit = expon.pdf(x,param)
            fit *= size
          
            taulist.append(param) 

        hist,bins = np.histogram(taulist, bins=20)
        width = 0.7*(bins[1]-bins[0])
        center = (bins[:-1]+bins[1:])/2
        #plt.bar(center, hist, align = 'center', width = width) produces bar graph
        plt.step(center, hist, where = 'mid')
        plt.savefig(inspect.stack()[0][3]+".png")
        

    def histDemo(self):
        mu, sigma = 200, 25
        x = mu + sigma*P.randn(10000)
        n, bins, patches = P.hist(x, 50, normed=1, histtype='stepfilled')
        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        
        y = P.normpdf( bins, mu, sigma)
        l = P.plot(bins, y, 'k--', linewidth=1.5)
        P.savefig(inspect.stack()[0][3]+".png")


if __name__ == '__main__':
    unittest.main()
