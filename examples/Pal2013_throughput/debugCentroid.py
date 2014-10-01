import time
import os
from util.FileName import FileName
from util.ObsFile import ObsFile
from util.utils import makeMovie
from util.utils import plotArray
import numpy as np
import astrometry.CentroidCalc as cent
import hotpix.hotPixels as hp
from photonlist import photlist
import matplotlib.pylab as plt

run = 'PAL2013'
date = '20131209'
tstamps = ['20131209-115553','20131209-120100', '20131209-121634', '20131209-122152']

for tstamp in tstamps:
    centroidTstamp = tstamp
    calTstamp='20131209-132225'
    fluxTstamp='ones'
    flatTstamp='20131209'

    obsFileName = FileName(run=run, date=date, tstamp=tstamp)
    obsFile = ObsFile(obsFileName.obs())
    frames = []
    print '-----------------------'
    print tstamp
    for i in xrange(30):
        print "making frame ", i
        frame = obsFile.getFrame(10*i,10)
        #plt.matshow(frame,vmin=10000,vmax=40000,cmap=plt.cm.gray)
        #plt.show()
        frames.append(frame)

    makeMovie(frames,frameTitles=np.arange(30)*10,cbar=True,outName='%s.gif'%(tstamp), normMin=400, normMax=4000)

    CentFile = FileName(run=run,date=date,tstamp=tstamp).centroidList()
    obsFile.loadCentroidListFile(CentFile)

    xs = obsFile.centroidListFile.getNode('/centroidlist/xPositions').read()
    ys = obsFile.centroidListFile.getNode('/centroidlist/yPositions').read()

    plt.figure()
    ax1 = plt.subplot(111)
    plt.plot(np.arange(len(xs)),xs,label="xs")
    plt.plot(np.arange(len(ys)),ys,label='ys')
    plt.legend()
    plt.title(tstamp)
    #plt.show()
    plt.savefig(tstamp)
