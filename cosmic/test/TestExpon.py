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
from scipy import optimize
from scipy.optimize import curve_fit
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

        hist,bins = np.histogram(taulist, bins=20, range=(15,25))
        width = 0.7*(bins[1]-bins[0])
        center = (bins[:-1]+bins[1:])/2
        plt.step(center, hist, where = 'post')
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

        hist,bins = np.histogram(taulist, bins=20, range=(15,35))
        width = 0.7*(bins[1]-bins[0])
        center = (bins[:-1]+bins[1:])/2
        #plt.bar(center, hist, align = 'center', width = width) produces bar graph
        plt.step(center, hist, where = 'post')
        plt.savefig(inspect.stack()[0][3]+".png")

    def displayFits(self):
        """
        generates two histograms on the same plot. One uses maximum likelihood to fit
        the data while the other uses the average time.
        """
        tau = 25.0
        nBins = 400
        size = 100
        taulist = []
        taulistavg = []
        for i in range(1000):
            x = range(nBins)
            timeHgValues = np.zeros(nBins, dtype=np.int64)
            timeStamps = expon.rvs(loc=0, scale=tau, size=size)
            ts64 = timeStamps.astype(np.uint64)
            tsBinner.tsBinner(ts64, timeHgValues)
            
            param = sum(timeStamps)/len(timeStamps)
            fit = expon.pdf(x,param)
            fit *= size
            
            taulistavg.append(param)

        for i in range(1000):
            x = range(nBins)
            timeHgValues = np.zeros(nBins, dtype=np.int64)
            timeStamps = expon.rvs(loc=0, scale=tau, size=size)
            ts64 = timeStamps.astype(np.uint64)
            tsBinner.tsBinner(ts64, timeHgValues)
            
            param = expon.fit(timeStamps)
            fit = expon.pdf(x,loc=param[0],scale=param[1])
            fit *= size
            taulist.append(param[1]) 


        hist,bins = np.histogram(taulistavg, bins=20, range=(15,35))
        width = 0.7*(bins[1]-bins[0])
        center = (bins[:-1]+bins[1:])/2
        plt.step(center, hist, where = 'post', label="averagetime", color='g')
        hist,bins = np.histogram(taulist, bins=20, range=(15,35))
        width = 0.7*(bins[1]-bins[0])
        center = (bins[:-1]+bins[1:])/2
        plt.step(center, hist, where = 'post', label="maxlikelihood")
        plt.legend()
        plt.savefig(inspect.stack()[0][3]+".png")
    
    def chisquaredDemo(self):
        """
        Demonstrate how to use curve_fit.  Generate data points with an
        exponential function, deviating with random gausian numbers.

        Use curve_fit to fit to these data points.

        Plot the data points and the fit curve to a png file
        """
            
        def func(x, a, b, c):
            return a*np.exp(-b*x) + c

        x = np.linspace(0, 4, 50)
        y = func(x, 2.5, 1.3, 0.5)
        yn = y + 0.2*np.random.normal(size=len(x))
            
        popt, pcov = curve_fit(func, x, yn)
        print 'optimal parameters: ', popt
        print 'uncertainties of parameters: ', pcov

        # plot the data points
        plt.plot(x,yn,'ro')

        # plot the fit line
        xFit = np.linspace(0,4,100)
        yFit = func(xFit, *popt)
        plt.plot(xFit,yFit)
        plt.title(inspect.stack()[0][3])
        plt.savefig(inspect.stack()[0][3]+".png")
        
    def chisquaredDemo2(self):
        """
        Demonstrate how to use curve_fit.  Generate data points with an
        exponential function, deviating with random gausian numbers.

        Use curve_fit to fit to these data points.

        Plot the data points and the fit curve to a png file
        """
            
        def func(x, a, b, c):
            return a*np.exp(-b*x) + c

        x = np.linspace(0, 4, 50)
        y = func(x, 2.5, 1.3, 0.5)
        yn = y + 0.2*np.random.normal(size=len(x))
            
        popt, pcov = curve_fit(func, x, yn)
        print 'optimal parameters: ', popt
        print 'uncertainties of parameters: ', pcov

        # plot the data points
        plt.plot(x,yn,'ro')

        # plot the fit line
        xFit = np.linspace(0,4,100)
        yFit = func(xFit, *popt)
        plt.plot(xFit,yFit)
        plt.title(inspect.stack()[0][3])
        plt.savefig(inspect.stack()[0][3]+".png")
        
    def testExponchisquared(self):
       
        def funcExpon(x, a, b, c):
            return a*np.exp(-b*x) + c

        tau = 25.0
        nBins = 400
        size = 10000
        taulist = []
        xPoints = []
        yPoints = []
        for i in range(10):
            timeHgValues = np.zeros(nBins, dtype=np.int64)
            timeStamps = expon.rvs(loc=0, scale=tau, size=size)
            ts64 = timeStamps.astype(np.uint64)
            tsBinner.tsBinner(ts64, timeHgValues)
           
        if (i == 0):
            print "first iteration"
            plt.clf()
            plt.plot(timeHgValues)
            plt.savefig("junk.png")

           
        for x in range(nBins):
            y = timeHgValues[x]
            if y > 2:
                xPoints.append(range(nBins))
                yPoints.append(timeHgValues)
                xArray = np.asarray(xPoints, dtype=list)
                yArray = np.asarray(yPoints, dtype=list)
                ySigma = (yArray ** 1/2)
                print "before call to curve_fit:  xArray",xArray
                print "before call to curve_fit:  yArray", yArray
                print "before call to curve_fit:  ySigmaArray",ySigma
                popt, pcov = curve_fit(funcExpon, xArray, yArray, sigma=ySigma)
                print "popt=",popt
                print 'optimal parameters: ', popt
                print 'uncertainties of parameters: ', pcov
        xFit = xPoints
        yFit = funcExpon(xFit, *popt)
        plt.clf()
        plt.plot(xPoints, yPoints, 'ro', label="data")
        plt.plot(xFit, yFit, color='g', label="fit")
        #plt.errorbar(yerr,  color='g', label='error')
        plt.legend()
        plt.title(inspect.stack()[0][3])
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
