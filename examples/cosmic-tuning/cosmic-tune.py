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

nSigmas    = [10,   10, 1000, 10]
strides    = [10,   10, 10,   10]
thresholds = [5000, 30, 30,   15]
styles     = ['b-', 'ro', 'g+', 'cx']
populationMax=1000
ySum = np.zeros(populationMax)
frameSum = 'none'

beginTime = 123
endTime = 133
plt.clf()
for nSigma,stride,threshold,style in zip(nSigmas, strides, thresholds, styles):
    print "seq=",seq
    fn = FileName(run, sundownDate, obsDate+"-"+seq)
    #cosmic = Cosmic(fn, endTime='exptime')
    cosmic = Cosmic(fn, beginTime=beginTime, endTime=endTime, 
                    loggingLevel=logging.CRITICAL)

    fc = cosmic.findCosmics(stride=stride, 
                            threshold=threshold, 
                            populationMax=populationMax,
                            nSigma=nSigma)

    tMasked = Cosmic.countMaskedBins(fc['interval'])
    ppmMasked = 1000000*tMasked/(endTime-beginTime)
    print "ppmMasked=",ppmMasked

    cosmic.file.cosmicMask = fc['interval']
    cosmic.file.cosmicMaskIsApplied = True

    fc = cosmic.findCosmics(stride=stride, 
                            threshold=threshold, 
                            populationMax=populationMax,
                            nSigma=nSigma)
    print "this should be zero:  ",Cosmic.countMaskedBins(fc['interval'])
    populationSum = np.array(fc['populationHg'][0])

    del cosmic

    label="%4d %2d %4d %4d"%(int(nSigma), stride, threshold, int(ppmMasked))
    print "label=",label
    plt.plot(populationSum, style, drawstyle="steps-mid", label=label)

fp = matplotlib.font_manager.FontProperties(family="monospace")

plt.xscale("symlog", linthreshx=0.9)
plt.yscale("symlog", linthreshy=0.1)
plt.xlim(-0.1, populationMax)
plt.ylim(ymin=-0.1)
plt.legend(loc="center right",prop=fp).get_frame().set_alpha(0.5)
plt.text(1e3,2e3,"nSigma stride threshold ppmMasked",horizontalalignment='right', fontproperties=fp)
plt.title("obsDate=%s seq=%s beginTime=%.1f endTime=%.1f"%(obsDate,seq,beginTime,endTime))
plt.xlabel("Population -- number of photons in a time bin")
plt.ylabel("dN/dPopulation")
plt.savefig("cosmic-tune.png")
