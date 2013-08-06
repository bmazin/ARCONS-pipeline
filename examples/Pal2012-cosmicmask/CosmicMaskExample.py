import numpy as np
import math
import unittest
from cosmic.Cosmic import Cosmic
import matplotlib.pyplot as plt
from interval import interval, inf, imath
from util import ObsFile
from scipy.stats import poisson, gamma, expon
from util import FileName
from scipy.optimize import curve_fit
from cosmic import tsBinner
"""
run = PAL2012
sundownDate = 20121211
obsDate = 20121212
seq = 121229

Makes a plot of populationHg values, masking out cosmic ray events. The poisson
distribution is shown with red crosses. The plot is saved as 
CosmicTimeMasking.png

Makes a plot for the entire file showing the mean for every ten seconds of
data and the cosmic events. This plot is saved as MeanAndCosmicEvents.png

Makes a scatter plot of the number of cosmics masked out and the mean which is
saved as scatter.png
"""

#loads the file
run = 'PAL2012'
sundownDate = '20121211'
obsDate = '20121212'
seq = '121229'
fn = FileName.FileName(run, sundownDate, obsDate+"-"+seq)
#this particular plot is for only ten seconds of data
cosmic = Cosmic(fn, beginTime=123, endTime=133)
dictionary0 = cosmic.findCosmics(nSigma=6, threshold=15, populationMax=1000)


interval = dictionary0['interval']
ObsFile.ObsFile.writeCosmicIntervalToFile(interval, 
                                          cosmic.file.ticksPerSec,
                                          'junk.h5')

cosmic.file.loadCosmicMask('junk.h5')
dictionary1 = cosmic.findCosmics(nSigma=6, threshold=15, populationMax=1000)
nSigma = 6
threshold = 15

plt.clf()
hist0 = dictionary0['populationHg'][0]
bins0 = np.arange(len(hist0))
plt.plot(bins0, hist0, label="no mask")

hist1 = dictionary1['populationHg'][0]
bins1 = np.arange(len(hist1))
plt.plot(bins1, hist1, label="cosmic mask")


mu1 = (bins1*hist1).sum()/float(hist1.sum())
p = poisson(mu1)
xvalues = np.arange(20)
theory = hist1.sum()*p.pmf(xvalues)
plt.plot(xvalues, theory,'+',label="poisson with mu=%.3f"%mu1)

plt.xscale('log')
plt.yscale('symlog',linthreshy=0.5)
plt.legend()
plt.title("Ten Seconds of Data;"+"  "+"threshold=%.3f"%threshold+"  "+ 
          "nSigma=%.3f"%nSigma)
plt.savefig("CosmicTimeMasking.png")
plt.clf()

#now make a plot for the whole file, showing the mean and the number of
#cosmic events for every ten seconds of data


nlist = [20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
         150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270,
         280, 290, 300]
mean = []
meanErr = []
nIntervals = []
nSigma = 7
threshold = 25
stride = 10
populationMax = 1000

#nlist = nlist[0:4]

for n in nlist:
   
    beginTime = n-10
    endTime = n
    cosmic = Cosmic(fn, beginTime, endTime)
    dictionary0 = cosmic.findCosmics(nSigma, threshold, stride, populationMax)
    interval = dictionary0['interval']
    nIntervals.append(len(interval))


    hist0 = dictionary0['populationHg'][0]
    bins0 = np.arange(len(hist0))

    muNum = (bins0*hist0).sum()
    muDen = float(hist0.sum())
    mu1 = muNum/muDen
    mu1Err = np.sqrt(mu1/muDen)
    print "mu1=",mu1," mu1Err=",mu1Err," muNum=",muNum," muDen=",muDen
    mean.append(mu1)
    meanErr.append(mu1Err)

plt.clf()
plt.subplot(211)
plt.errorbar(nlist, mean, yerr=meanErr, fmt='x')
plt.ylabel("mean photons/tick")
plt.title("Mean and Interval Time"+"  "+"threshold=%.3f"%threshold+"  "+ 
          "nSigma=%.3f"%nSigma)

nIntervalsErr = np.sqrt(nIntervals)
plt.subplot(212)
plt.errorbar(nlist, nIntervals, yerr=nIntervalsErr, fmt='x')
plt.xlabel("time in exposure (sec)")
plt.ylabel("# cosmics masked")
plt.savefig("MeanAndCosmicEvents.png")

plt.clf()
plt.scatter(nIntervals, mean)
plt.xlabel("# cosmics masked")
plt.ylabel("mean photons/tick")
plt.savefig("scatter.png")
