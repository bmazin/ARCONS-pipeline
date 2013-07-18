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
cosmic = Cosmic(fn, beginTime=123, endTime=124)
dictionary0 = cosmic.findCosmics(threshold=30, populationMax=1000)


interval = dictionary0['interval']
ObsFile.ObsFile.writeCosmicIntervalToFile(interval, 
                                          cosmic.file.ticksPerSec,
                                          'junk.h5')

cosmic.file.loadCosmicMask('junk.h5')
dictionary1 = cosmic.findCosmics(threshold=30, populationMax=1000)

plt.clf()
hist0 = dictionary0['populationHg'][0]
bins0 = np.arange(len(hist0))
plt.plot(bins0, hist0, label="no mask")

hist1 = dictionary1['populationHg'][0]
bins1 = np.arange(len(hist1))
plt.plot(bins1, hist1, label="cosmic mask")


s = np.random.poisson(np.mean(dictionary1['populationHg'][0]), 1000)
plt.plot(s, hist1)

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("testCosmicTimeMasking.png")

        
