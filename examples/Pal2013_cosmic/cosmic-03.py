import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
from util import hgPlot
from cosmic.Cosmic import Cosmic
import tables
from hotpix import hotPixels
import pickle
from interval import interval, inf, imath
import logging, os
import pickle
run = 'PAL2013'
sundownDate = '20131204'
obsDate = '20131205'
seq = '033506'

populationMax=1000
beginTime = 0
endTime = 100

stride = 1000
fn = FileName(run, sundownDate, obsDate+"-"+seq)
cosmic = Cosmic(fn, beginTime=beginTime, endTime=endTime, loggingLevel=logging.INFO)
timeHgValues,frameSum = cosmic.getTimeHgAndFrameSum(beginTime,endTime)

lthgv=len(timeHgValues)
print "length of timeHgValues=",lthgv
nStride=lthgv/stride

binnedTimeHgValues = np.reshape(timeHgValues[0:nStride*stride],[nStride,stride]).sum(axis=1)
#plt.xkcd()
dt = cosmic.file.tickDuration*stride
times = cosmic.beginTime + dt*np.arange(nStride)
plt.plot(times,binnedTimeHgValues,
         drawstyle='steps-mid',label="all data")
#plt.xlim(xmax=1000)
#plt.xscale('symlog',linthreshx=0.9)
plt.yscale('symlog',linthreshy=0.5)
plt.title('%s %s %s %s'%(run,sundownDate,obsDate,seq))
plt.xlabel("time")
plt.ylabel("number of photons")
plt.legend()
plt.savefig("cosmic-03.png")

