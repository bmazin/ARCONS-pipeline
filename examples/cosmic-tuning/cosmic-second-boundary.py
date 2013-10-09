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
# December 11
# Final sequence, toward end of run, thin high clouds at around 12:50, moved to confirm position at '122234', also at '130752' at 125s. (16,15)

#seq5 = ['112709', '113212', '113714', '114216', '114718', '115220', '115722', '120224', '120727', '121229', '121732', '122234', '122736', '123238', '123740', '124242', '124744', '125246', '125748', '130250', '130752', '131254', '131756', '132258', '132800', '133303']

seq = '121229'

populationMax=1000

beginTime = 0
endTime = 'exptime'
plt.clf()
nSigma = 10
stride = 10
threshold = 15

bin0 = 23
print "bin0=",bin0

fn = FileName(run, sundownDate, obsDate+"-"+seq)
cosmic = Cosmic(fn, beginTime=beginTime, endTime=endTime, 
                loggingLevel=logging.CRITICAL)

fc = cosmic.findCosmics(stride=stride, 
                        threshold=threshold, 
                        populationMax=populationMax,
                        nSigma=nSigma)

tList = []
for t in range(beginTime,int(cosmic.endTime)):
    tList.append(t)

nPhotonLists = {}
offsets = [0, 100, 1234, 54321]

for offset in offsets:
    nPhotonLists[offset] = []
    for t in range(beginTime,int(cosmic.endTime)):
        bin = cosmic.file.ticksPerSec*(t-beginTime) + bin0 + offset
        nPhoton = fc['timeHgValues'][bin-2:bin+3].sum()
        nPhotonLists[offset].append(fc['timeHgValues'][bin])

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
