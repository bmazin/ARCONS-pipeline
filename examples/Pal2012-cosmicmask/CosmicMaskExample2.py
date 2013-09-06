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
meanErr = []
mean2 = []
meanErr2 = []
nSigma = 7
threshold = 25
stride = 10
populationMax = 1000

nlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
         21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 
         38, 39, 40]

nlist2 = [10, 20, 30]

for n in nlist:
    
    beginTime = n-1
    endTime = n
    cosmic = Cosmic(fn, beginTime, endTime)
    dictionary0 = cosmic.findCosmics(nSigma, threshold, stride, populationMax)
    interval = dictionary0['interval']
   
    hist0 = dictionary0['populationHg'][0]
    bins0 = np.arange(len(hist0))

    
    muNum = (bins0*hist0).sum()
    muDen = float(hist0.sum())
    mu1 = muNum/muDen
    mu1Err = np.sqrt(mu1/muDen)
    print "mu1=",mu1," mu1Err=",mu1Err," muNum=",muNum," muDen=",muDen
    mean.append(mu1)
    meanErr.append(mu1Err)

for n in nlist2:
    
    beginTime = n-10
    endTime = n
    cosmic = Cosmic(fn, beginTime, endTime)
    dictionary0 = cosmic.findCosmics(nSigma, threshold, stride, populationMax)
    interval = dictionary0['interval']
   
    hist0 = dictionary0['populationHg'][0]
    bins0 = np.arange(len(hist0))

    
    muNum = (bins0*hist0).sum()
    muDen = float(hist0.sum())
    mu1 = muNum/muDen
    mu1Err = np.sqrt(mu1/muDen)
    mean2.append(mu1)
    meanErr2.append(mu1Err)

plt.errorbar(nlist, mean, yerr=meanErr, fmt= 'x', label="mean over one second")
plt.errorbar(nlist2, mean2, yerr=meanErr2, fmt='rx', 
             label="mean over ten seconds")
plt.ylabel("mean photons/tick")
plt.xlabel("seconds")
plt.legend(loc=2)
plt.title("run=PAL2012-sundownDate=20121211-obsDate=20121212-seq=121229")
plt.savefig("Mean.png")
