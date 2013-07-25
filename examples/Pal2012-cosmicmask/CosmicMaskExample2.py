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

Makes a plot showing the mean of the distribution for ten seconds of data
"""


run = 'PAL2012'
sundownDate = '20121211'
obsDate = '20121212'
seq = '121229'
fn = FileName.FileName(run, sundownDate, obsDate+"-"+seq)
mean = []
nSigma = 7
threshold = 25
stride = 10
populationMax = 1000

nlist = [31, 32, 33, 34, 35, 36, 37, 38, 39, 40]

for n in nlist:
    
    cosmic = Cosmic(fn, beginTime=n-1, endTime=n)
    dictionary0 = cosmic.findCosmics(nSigma, threshold, stride, populationMax)
    interval = dictionary0['interval']
   
    hist0 = dictionary0['populationHg'][0]
    bins0 = np.arange(len(hist0))

    
    muNum = (bins0*hist0).sum()
    muDen = float(hist0.sum())
    mu1 = muNum/muDen
    print "mu1=",mu1, " muNum=",muNum," muDen=",muDen
    mean.append(mu1)


plt.plot(nlist, mean, 'x')
plt.ylabel("mean photons/tick")
plt.title("Mean"+"  "+"threshold=%.3f"%threshold+"  "+"nSigma=%.3f"%nSigma)
plt.savefig("Mean.png")
