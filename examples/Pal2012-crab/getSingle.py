from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
import matplotlib.pyplot as plt
import numpy as np
import datetime
import tables
import time

def circ(startpx,startpy,radius=3):
    r = radius
    length = 2*r
    height = length
    allx = xrange(startpx-int(np.ceil(length/2.0)),startpx+int(np.floor(length/2.0))+1)
    ally = xrange(startpy-int(np.ceil(height/2.0)),startpy+int(np.floor(height/2.0))+1)
    pixx = []
    pixy = []
    for x in allx:
        for y in ally:
            if (np.abs(x-startpx))**2+(np.abs(y-startpy))**2 <= (r)**2:
                pixx.append(x)
                pixy.append(y)
    return pixx,pixy

def main():
    d = datetime.date(2012,10,15)#date of table value
    d2 = datetime.date(2012,12,11)#date of observation
    dt=(d2-d).total_seconds()#seconds between obs and table value
    nu=29.6957720714 #1/us
    pdot=4.2013878e-13#us/s
    period=(1.0/nu)+pdot*dt#us
    print 'period=',period,'s'

    centerRow = 30
    centerCol = 30
    circCol,circRow = circ(centerCol,centerRow,radius=5)
    integrationTime = 60
    totalTime = 300
    steps = totalTime/integrationTime
    
    firstSec = 0
    nPhaseBins = 400
    nFoldedPulses = 1.0*integrationTime/period
    #np.set_printoptions(threshold=np.nan)
    ob = ObsFile('obs_20121212-033323.h5')
    ut = ob.getFromHeader('unixtime')
    timestamps = np.array([])
    plt.ion()
    
    for i in range(len(circCol)):
        iRow = circRow[i]
        iCol = circCol[i]
        times,peaks,baselines = ob.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
        times += ut
        timestamps = np.append(timestamps,times)
        #phases = (times % period) / period
        #histPhases,phaseBinEdges = np.histogram(phases,bins=nPhaseBins)
        #plt.plot(phaseBinEdges[0:-1],histPhases)
        #plt.show()
        #time.sleep(.1)
        #plt.close()
        
    plt.ioff()
    #timestamps*=1.0e6
    phases = (timestamps % period) / period
    histPhases,phaseBinEdges = np.histogram(phases,bins=nPhaseBins)
    histPhases = np.array(histPhases,dtype='float')
    histPhases /= nFoldedPulses
    print 'avg counts per period (sky included):',sum(histPhases)
    sdev = np.std(histPhases)
    med = np.median(histPhases)
    
    sdev = np.std(histPhases[histPhases<med+sdev])
    plt.plot(phaseBinEdges[0:-1][histPhases< med+sdev],histPhases[histPhases<med+sdev])
    baseline = np.median(histPhases[histPhases<med+sdev])

    histPhases -= baseline
    plt.plot(phaseBinEdges[0:-1],histPhases)
    plt.plot([phaseBinEdges[0],phaseBinEdges[-2]],[0,0])
    print 'avg counts per period:',sum(histPhases)
    print 'baseline',baseline
        

    pulseIdx = 100
    startTime = ut
    photonsInPeriod = []
    for pulseIdx in range(int(nFoldedPulses)):
        pulseStartTime = ut+pulseIdx*period
        pulseEndTime = ut+(pulseIdx+1)*period
        selectedTimestamps = timestamps[np.logical_and(pulseStartTime < timestamps,timestamps < pulseEndTime)]
        photonsInPeriod.append(len(selectedTimestamps))

       
    plt.show()
    photonsInPeriod = np.array(photonsInPeriod)
    plt.hist(photonsInPeriod,bins=100)
    plt.show()


if __name__=='__main__':
    main()
