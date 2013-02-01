from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
import matplotlib.pyplot as plt
import numpy as np
import datetime
import tables
import ephem

def main():
    #filenames = ['obs_20121212-033825.h5',
    ob = ObsFile('obs_20121212-033825.h5')
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
    startTime = obsDate*24*3600#in seconds
    dt = startTime-pEpoch#seconds since pepoch

    #freq = F0+F1*dt+F2/2*dt**2
    freq = F0+F1*dt
    period = 1.0/freq
    print 'period=',period,'s'


    #period=0.03367660643405371
    #period=0.03367664238573182
    #period=0.03367662440988317

    iRow = 30
    iCol = 30
    
    for i in range(10):
        firstSec = 0+i*30
        dt = startTime-pEpoch + firstSec
        freq = F0+F1*dt
        period = 1.0/freq
        print 'period=',period,'s'

        integrationTime = 30
        nPhaseBins = 200
        #np.set_printoptions(threshold=np.nan)
        
        timestamps,peaks,baselines = ob.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
        #timestamps*=1.0e6
        #print '\n',timestamps[0:10]
        phases = (timestamps % period)
        #print phases[0:10]
        #print len(phases)
        histPhases,phaseBinEdges = np.histogram(phases,bins=nPhaseBins)
        c = (i/10.0,0,1.0-i/10.0)
        
       

        plt.plot(phaseBinEdges[0:-1],histPhases,color=c)
    plt.show()

if __name__=='__main__':
    main()
