'''
Author: Julian van Eyken    Date: May 7 2013

Test code for photon lists. Under construction....
'''
import time
from util.FileName import FileName
from util.ObsFile import ObsFile
from util.utils import plotArray
import numpy as np
from photonlist import photlist


def testWritePhotonList(outputFileName=None,firstSec=0,integrationTime=-1):
    '''
    Test run of obsFile.writePhotonList. fileName can be used
    to specify the output file name. If not specified, default
    name/location is used.
    '''

    #Details of example obs file to run test on.
#    run = 'PAL2012'
#    date = '20121207'
#    tstamp = '20121208-074649'
#    calTstamp='20121208-070505'
#    fluxTstamp='20121208-133002'
    run = 'PAL2012'
    date = '20121208'
    tstamp = '20121209-120530'
    centroidTstamp = '20121209-120530'
    calTstamp='20121209-131132'
    fluxTstamp='20121209-020416'
    flatTstamp='20121209-021036'
    
    
    #Load up the obs file
    obsFileName = FileName(run=run, date=date, tstamp=tstamp)
    obsFile = ObsFile(obsFileName.obs())
    
    #Load up associated calibrations
    obsFile.loadWvlCalFile(FileName(run=run,date=date,tstamp=calTstamp).calSoln())
    obsFile.loadFlatCalFile(FileName(run=run,date=date,tstamp=flatTstamp).flatSoln())
    obsFile.loadFluxCalFile(FileName(run=run,date=date,tstamp=fluxTstamp).fluxSoln())
    obsFile.loadTimeAdjustmentFile(FileName(run=run,date=date,tstamp=tstamp).timeAdjustments())
    obsFile.loadHotPixCalFile(FileName(run=run,date=date,tstamp=tstamp).timeMask())
    
    #And write out the results....
    obsFile.writePhotonList(outputFileName,firstSec,integrationTime,
                            astrometryFileName=FileName(run=run,date=date,tstamp=tstamp).centroidList())
    
    #Read the results back in....
    #photFile = photList.PhotFile(outputFilename)
    
    

    
if __name__ == "__main__":
    
    import sys
    #hotPixelsTest(sys.argv[1])
    #hotPixelsTest2(startTime=2.5, endTime=6.3)
    TestWritePhotonList()
    
    
