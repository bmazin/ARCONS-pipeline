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
expTime = 300


fn = FileName(run, sundownDate, obsDate+"-"+seq)
cosmic = Cosmic(fn, beginTime=0, endTime=10, loggingLevel=logging.INFO)
fcd = cosmic.findCosmics(nSigma=10, stride=10, threshold=15)
print "done:  fcd.keys=",fcd.keys()
maskedTime = Cosmic.countMaskedBins(fcd['interval'])
maskedPercent = 100*maskedTime/(cosmic.endTime-cosmic.beginTime)
print "maskedPercent=",maskedPercent

cosmic.file.setCosmicMask(fcd['interval'])
fcd2 = cosmic.findCosmics()

#plt.xkcd() 
plt.plot(fcd['populationHg'][0],drawstyle='steps-mid',label="all data")
plt.plot(fcd2['populationHg'][0],drawstyle='steps-mid',label="cosmics masked")
plt.xlim(xmax=1000)
plt.xscale('symlog',linthreshx=0.9)
plt.yscale('symlog',linthreshy=0.5)
plt.title('%s %s %s %s'%(run,sundownDate,obsDate,seq))
plt.xlabel("number of photons in one tick")
plt.ylabel("number of ticks")
plt.legend()
plt.savefig("cosmic-02.png")
