import numpy as np
import astrometry.CalculateRaDec as crd
from util.FileName import FileName

def testCalcRaDec():
    run = 'PAL2012'
    sunsetDate='20121208'
    utcDate='20121209'
    centroidTimestamp = '20121209-120530'
    calTimestamp = '20121209-131132'
    centroidListFileName=FileName(run=run,date=sunsetDate,tstamp=centroidTimestamp).centroidList()

    # Test photon
    xPhotonPixel=np.linspace(-10,10, num = 500)
    yPhotonPixel=np.linspace(-10,10,num = 500)
    timestamp = np.linspace(5,5,num = 500)
    
    #tic = time()
    raDecObject = crd.CalculateRaDec(centroidListFileName)

    ras,decs,has = raDecObject.getRaDec(timestamp=timestamp,xPhotonPixel=xPhotonPixel,yPhotonPixel=yPhotonPixel)
    print 'RA:',ras
    print 'Dec: ',decs
    print 'HA:',has
    print
    #print 'Time taken (s): ',time()-tic
    return xPhotonPixel,yPhotonPixel,ras,decs,has



if __name__ == "__main__":
    testCalcRaDec()
