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


# test the theory that that extra photons show up 23 clock ticks
# after the beginning of each second.  Plot properties of the photons

offsets = [23, 100, 1234, 54321]
#offsets = [23]
beginTime = 0
expTime = 300

pickleFile = open("csb2.pkl","wb")

fn = FileName(run, sundownDate, obsDate+"-"+seq)
pickle.dump(fn,pickleFile)
obsFile = ObsFile(fn.obs())
#obsFile.loadTimeAdjustmentFile(FileName(run=run).timeAdjustments(),verbose=True) # Matt's time fix
timeMaskFile = fn.timeMask()
if os.path.exists(timeMaskFile):
    obsFile.loadHotPixCalFile(timeMaskfile, switchOnMask=True)

if False:
    iRow0 = 25
    iRow1 = iRow0+1
    iCol0 = 40
    iCol1 = iCol0+1
else:
    iRow0 = 0
    iRow1 = obsFile.nRow
    iCol0 = 0
    iCol1 = obsFile.nCol


spt = obsFile.tickDuration
print "seconds per tick=",spt
masks = {}
pickle.dump(offsets,pickleFile)
for offset in offsets:
    print "begin offset=",offset
    masks[offset] = interval()

    for sec in range(beginTime,beginTime+expTime):
        tl = sec+obsFile.tickDuration*(offset-2)
        tr = sec+obsFile.tickDuration*(offset+3)
        masks[offset] = masks[offset] | interval([sec,tl]) | interval([tr,sec+1])

    obsFile.cosmicMask = masks[offset]
    obsFile.switchOnCosmicTimeMask()

    nPhotonSum = 0
    rows = np.zeros(0,dtype=np.int)
    cols = np.zeros(0,dtype=np.int)
    dts  = np.zeros(0,dtype=np.double)
    secs = np.zeros(0,dtype=np.int)
    phs  = np.zeros(0,dtype=np.double)
    for iRow in range(iRow0, iRow1):
        for iCol in range(iCol0, iCol1):
            #tpl = obsFile.getTimedPacketList(iRow, iCol, firstSec=beginTime, integrationTime=expTime)
            tpl = obsFile.getPackets(iRow, iCol, firstSec=beginTime, integrationTime=expTime, fields=['peakHeights'])
            nPhoton = len(tpl['timestamps'])
            if nPhoton>0:
                print "offset=",offset,"row=",iRow, "iCol=",iCol, "nPhotons=",nPhoton
                rows = np.append(rows, iRow*np.ones(nPhoton,dtype=np.int))
                cols = np.append(cols, iCol*np.ones(nPhoton,dtype=np.int))
                dts  = np.append(dts, tpl['timestamps']-tl)
                secs = np.append(secs, tpl['timestamps'].astype(np.int))
                phs  = np.append(phs, tpl['peakHeights'])
            nPhotonSum += nPhoton
    pickle.dump({"offset":offset,"rows":rows,"cols":cols,"dts":dts,"secs":secs,"phs":phs}, pickleFile)
    print "keys=",tpl.keys()
    print "photonSum=",nPhotonSum
#plt.savefig("csb2.png")
pickleFile.close()

del obsFile
