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
from scipy.optimize import curve_fit
import pickle
from interval import interval, inf, imath
from cosmic import tsBinner
import sys
import os
run = 'PAL2012'
sundownDate = '20121211'
obsDate = '20121212'
# December 11
# Final sequence, toward end of run, thin high clouds at around 12:50, moved to confirm position at '122234', also at '130752' at 125s. (16,15)
seq5 = ['112709', '113212', '113714', '114216', '114718', '115220', '115722', '120224', '120727', '121229', '121732', '122234', '122736', '123238', '123740', '124242', '124744', '125246', '125748', '130250', '130752', '131254', '131756', '132258', '132800', '133303']

#seq5 = ['120727']

stride = 10
threshold = 100
nAboveThreshold = 0
npList = []
sigList = []
   
run = 'PAL2012'
sundownDate = '20121211'
obsDate = '20121212'

for seq in seq5:
    inFile = open("cosmicTimeList-%s.pkl"%(seq),"rb")
    cosmicTimeList = pickle.load(inFile)
    binContents = pickle.load(inFile)
    cfn = "cosmicMax-%s.h5"%seq
    intervals = ObsFile.readCosmicIntervalFromFile(cfn)
    for interval in intervals:
        print "interval=",interval
        fn = FileName(run, sundownDate,obsDate+"-"+seq)
        obsFile = ObsFile(fn.obs())
        obsFile.loadTimeAdjustmentFile(fn.timeAdjustments())
        i0=interval[0]
        i1=interval[1]
        intervalTime = i1-i0
        dt = intervalTime/2
        beginTime = max(0,i0-0.000200)
        endTime = beginTime + 0.001
        integrationTime = endTime-beginTime
        nBins = int(np.round(obsFile.ticksPerSec*(endTime-beginTime)+1))
        timeHgValues = np.zeros(nBins, dtype=np.int64)
        ymax = sys.float_info.max/100.0
        for iRow in range(obsFile.nRow):
            for iCol in range(obsFile.nCol):
                gtpl = obsFile.getTimedPacketList(iRow,iCol,
                                                  beginTime,integrationTime)
                ts = (gtpl['timestamps'] - beginTime)*obsFile.ticksPerSec
                ts64 = np.round(ts).astype(np.uint64)
                tsBinner.tsBinner(ts64, timeHgValues)
        plt.clf()
        plt.plot(timeHgValues)
        x0 = (i0-beginTime)*obsFile.ticksPerSec
        x1 = (i1-beginTime)*obsFile.ticksPerSec
        plt.fill_between((x0,x1),(0,0), (ymax,ymax), alpha=0.2, color='red')
        plt.yscale("symlog",linthreshy=0.9)
        plt.xlim(0,1000)
        #plt.ylim(-0.1,timeHgValues.max())
        plt.ylim(-0.1,300)
        tick0 = int(np.round(i0*obsFile.ticksPerSec))
        plotfn = "cp-%05d-%s-%s-%s-%09d"%(timeHgValues.sum(),run,obsDate,seq,tick0)
        plt.title(plotfn)
        plt.savefig(plotfn+".png")
        print "plotfn=",plotfn
        

print "to make movie now running this command:"
print "convert -delay 0 `ls -r cp*png` cp.gif"
os.system("convert -delay 0 `ls -r cp*png` cp.gif")
