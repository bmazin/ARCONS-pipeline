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
import logging
run = 'PAL2012'
sundownDate = '20121211'
obsDate = '20121212'
seq = '121229'

populationMax=1000

beginTime = 0
endTime = 'exptime'
plt.clf()
nSigma = 10
stride = 10
threshold = 15

# test the theory that that extra photons show up 23 clock ticks
# after the beginning of each second
bin0 = 23 

# find cosmics and put the results in the dictionary fc.
# fc['timeHgvalues'] is an array.  The ith value of the array
# is the number of photons on the array during the ith clock tick after 
# beginTime

fn = FileName(run, sundownDate, obsDate+"-"+seq)
cosmic = Cosmic(fn, beginTime=beginTime, endTime=endTime, 
                loggingLevel=logging.CRITICAL)

fc = cosmic.findCosmics(stride=stride, 
                        threshold=threshold, 
                        populationMax=populationMax,
                        nSigma=nSigma)

# list of x values for plotting
tList = []
for t in range(beginTime,int(cosmic.endTime)):
    tList.append(t)

# count the number of photons each second, with different offsets from bin0
offsets = [0, 100, 1234, 54321]

# dictionary of y values to plot
nPhotonLists = {}

# do actual work here
for offset in offsets:
    nPhotonLists[offset] = []
    for t in range(beginTime,int(cosmic.endTime)):
        bin = cosmic.file.ticksPerSec*(t-beginTime) + bin0 + offset
        nPhoton = fc['timeHgValues'][bin-2:bin+3].sum()
        nPhotonLists[offset].append(fc['timeHgValues'][bin])

# make a plot
plt.clf()
for offset in offsets:
    plt.plot(tList, nPhotonLists[offset], label="offset=%d"%offset)
plt.legend().get_frame().set_alpha(0.5)
plt.title("obsDate=%s seq=%s beginTime=%.1f endTime=%.1f"%(obsDate,seq,beginTime,cosmic.endTime))
plt.xlabel("time")
plt.ylabel("number of photons")
plt.ylim(ymin=-2)
plt.savefig("cosmic-second-boundary.png")

del cosmic
