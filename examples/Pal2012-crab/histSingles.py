from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
from util.FileName import FileName
from util.popup import PopUp
from hotpix import hotPixelsMatt as hotPixels
import matplotlib.pyplot as plt
import numpy as np
import datetime
import tables
import ephem
import matplotlib

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
    #filenames = ['obs_20121212-033825.h5',
    fn = FileName(run='PAL2012',date='20121211',tstamp='20121212-055428').obs()
    ob = ObsFile(fn)
    frame = ob.getPixelCountImage(0,30,weighted=False)
    hotPixMask = hotPixels.findHotPixels(image=frame,firstSec=0,intTime=30,weighted=False)['badflag']
    frame[hotPixMask == 1] = 0
    def plotFrame(fig,axes):
        hMat=axes.matshow(frame,cmap=matplotlib.cm.gnuplot,origin='lower',vmax=np.mean(frame)+3*np.std(frame))
        fig.colorbar(hMat)
    PopUp(plotFunc=plotFrame)
    #d = datetime.date(2012,10,15)#date of table value
    #d2 = datetime.date(2012,12,11)#date of observation
    #dt=(d2-d).total_seconds()#seconds between obs and table value
    ##crabtime ephemeris
    #nu=29.6957720714 #1/us
    #pdot=4.2013878e-13#us/s
    #period=(1.0/nu)+pdot*dt#us

    #goldstone ephemeris
    F0=29.7414493465249770#Hz
    deltaF0=0.0000000055983574
    F1=-3.669005118205e-10#df0/dt Hz/s
    F2=-3.085573120457e-20#d^2f0/dt^2 Hz/s^2
    pEpoch = 54781.604891 #Modified Julian Date corresponding to F0
    pEpoch = pEpoch+2400000.5#convert mjd to jd
    pEpoch *= 24*3600 #in seconds
    obsDate = ob.getFromHeader('jd')
    print ob.getFromHeader('utc')
    startTime = obsDate*24*3600#in seconds
    dt = startTime-pEpoch#seconds since pepoch

    #freq = F0+F1*dt+F2/2*dt**2
    freq = F0+F1*dt
    period = 1.0/freq
    print 'period=',period,'s'


    #period=0.03367660643405371
    #period=0.03367664238573182
    #period=0.03367662440988317

    iRow = 10
    iCol = 14
    integrationTime = -1
    circCol,circRow = circ(iCol,iRow,radius=4)
    firstSec = 0
    
    dt = startTime-pEpoch + firstSec
    freq = F0+F1*dt
    period = 1.0/freq
    print 'period=',period,'s'

    nPhaseBins = 200
    #np.set_printoptions(threshold=np.nan)
    
    jdTimes = np.array([],dtype=np.float64)
    times = np.array([])
    for i in range(len(circCol)):
        iRow = circRow[i]
        iCol = circCol[i]
        timestamps,peaks,baselines = ob.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
        timestamps = np.array(timestamps,dtype=np.float64)
        jdTimestamps = obsDate+timestamps /(24.*3600.)
        jdTimes = np.append(jdTimes,jdTimestamps)
        times = np.append(times,timestamps)

    jdTimes -= 2400000.5 #convert to modified jd
    periodDays = period/(24.*3600.)
    phaseOffset = 0
    phases = (jdTimes % periodDays)/periodDays+phaseOffset
    #in case phaseOffset pushed a phase past 1. roll it back 
    phases = phases % 1.0 
    #determine which period each jd timestamp falls in
    periodIndices = (jdTimes+phaseOffset*periodDays)//periodDays
    #make folded phase hist
    histPhases,phaseBinEdges = np.histogram(phases,bins=nPhaseBins)
       

    plt.plot(phaseBinEdges[0:-1],histPhases)
    plt.show()

    firstPeriod = periodIndices.min()
    lastPeriod = periodIndices.max()
    nPeriods = lastPeriod-firstPeriod+1
    periodBins = np.arange(firstPeriod,lastPeriod+2)

    #Count how many timestamps/photons occur in each period
    singlePeriodCounts,periodBinEdges = np.histogram(periodIndices,bins=periodBins)
    plt.plot(periodBinEdges[:-1],singlePeriodCounts)
    plt.show()

    #make a histogram of counts in each period
    histSinglePeriodCounts,countBinEdges = np.histogram(singlePeriodCounts,bins=50)
    plt.plot(countBinEdges[:-1],histSinglePeriodCounts)
    plt.show()

    n = 30
    counts = []
    for i in range(n):
        jdTimesSelect = jdTimes[(periodIndices-firstPeriod)==2000+i]
        counts.append(len(jdTimesSelect))
        phases = (jdTimesSelect % periodDays)/periodDays+phaseOffset
        histPhases,phaseBinEdges = np.histogram(phases,bins=10)
        plt.plot(phaseBinEdges[0:-1],histPhases,c=matplotlib.cm.jet(1.*i/n))
    plt.show()
    plt.plot(counts)
    plt.show()
        

if __name__=='__main__':
    main()
