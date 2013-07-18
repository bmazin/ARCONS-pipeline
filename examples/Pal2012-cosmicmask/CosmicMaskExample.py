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
print "mu1=", mu1
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

        
