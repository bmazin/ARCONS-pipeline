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

Also makes a plot for the entire file showing the mean for every ten seconds of
data and the cosmic events. This plot is saved as MeanAndCosmicEvents.png
"""

run = 'PAL2012'
sundownDate = '20121211'
obsDate = '20121212'
seq = '121229'
fn = FileName.FileName(run, sundownDate, obsDate+"-"+seq)
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


nlist = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140,
         150, 150, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270,
         280, 290, 300]
mean = []
cosmicevents = []
nSigma = 7
threshold = 100
stride = 10
populationMax = 1000
for n in nlist:

    cosmic = Cosmic(fn, beginTime=n-10, endTime=n)
    dictionary0 = cosmic.findCosmics(nSigma, threshold, stride, populationMax)
    interval = dictionary0['interval']
    ObsFile.ObsFile.writeCosmicIntervalToFile(interval, 
                                          cosmic.file.ticksPerSec,
                                          'junk.h5')

    cosmic.file.loadCosmicMask('junk.h5')
    dictionary1 = cosmic.findCosmics(nSigma, threshold, stride, populationMax)

    hist1 = dictionary1['populationHg'][0]
    bins1 = np.arange(len(hist1))
    hist0 = dictionary0['populationHg'][0]

    mu1 = (bins1*hist1).sum()/float(hist1.sum())
    mean.append(mu1)
    events = (hist0*hist1).sum()-hist1.sum()
    cosmicevents.append(events)

plt.subplot(211)
plt.plot(nlist, mean, 'o')
plt.ylabel("mean")
plt.title("Mean and Cosmic Events"+"  "+"threshold=%.3f"%threshold+"  "+ 
          "nSigma=%.3f"%nSigma)

plt.subplot(212)
plt.plot(nlist, cosmicevents, 'o')
plt.xlabel("seconds")
plt.ylabel("cosmics")
plt.savefig("MeanAndCosmicEvents.png")
