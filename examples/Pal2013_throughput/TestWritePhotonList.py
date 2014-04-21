'''
Author: Julian van Eyken    Date: May 7 2013

Edited 3/2/14 SRM
altered for 2013 throughput calculation files' photonlist generation

Test code for photon lists. Under construction....
'''
import time
import os
from util.FileName import FileName
from util.ObsFile import ObsFile
from util.utils import plotArray
import numpy as np
import astrometry.CentroidCalc as cent
import hotpix.hotPixels as hp
from photonlist import photlist


def testWritePhotonList(outputFileName=None,firstSec=0,integrationTime=-1,doPixRemap=False):
    '''
    Test run of obsFile.writePhotonList. fileName can be used
    to specify the output file name. If not specified, default
    name/location is used.
    
    Now includes test for pixel remapping....
    '''

    #Details of example obs file to run test on.
#    run = 'PAL2012'
#    date = '20121207'
#    tstamp = '20121208-074649'
#    calTstamp='20121208-070505'
#    fluxTstamp='20121208-133002'

    run = 'PAL2013'
    date = '20131209'
    tstamp = '20131209-115553'
    centroidTstamp = tstamp
    calTstamp='20131209-132225'
    fluxTstamp='ones'
    flatTstamp='20131209'
    if doPixRemap==True:
        pixRemapFileName = FileName(run=run).pixRemap()
    else:
        pixRemapFileName = None
    
    #Load up the obs file
    obsFileName = FileName(run=run, date=date, tstamp=tstamp)
    obsFile = ObsFile(obsFileName.obs())

    obsFile.loadWvlCalFile(FileName(run=run,date=date,tstamp=calTstamp).calSoln())
    obsFile.loadFlatCalFile(FileName(run=run,date=date,tstamp=flatTstamp).flatSoln())
    obsFile.loadFluxCalFile(FileName(run=run,date=date,tstamp=fluxTstamp).fluxSoln())

    HotPixFile = FileName(run=run,date=date,tstamp=tstamp).timeMask()
    if not os.path.exists(HotPixFile): #check if hot pix file already exists
        hp.findHotPixels(obsFileName.obs(),HotPixFile)
        print "Flux file pixel mask saved to %s"%(HotPixFile)
    obsFile.loadHotPixCalFile(HotPixFile)
    print "Hot pixel mask loaded %s"%(HotPixFile)   

    #hz21 location alpha(2000) = 12h 13m 56.42s , delta(2000) = +32d 56' 30.8''
    centroid_RA = '12:13:56.42'
    centroid_DEC = '32:56:30.8'
    CentFile = FileName(run=run,date=date,tstamp=tstamp).centroidList()
    #print "checking centroid file", CentFile
    #if not os.path.exists(CentFile): #check if cent pix file already exists
    cent.centroidCalc(obsFile, centroid_RA, centroid_DEC, outputFileName = CentFile, guessTime=10, integrationTime=10)
    #    print "Flux file centroid cal saved to %s"%(CentFile)
    obsFile.loadCentroidListFile(CentFile)
    print "Centroid File loaded %s"%(CentFile) 

    #Load up associated calibrations

    #obsFile.loadTimeAdjustmentFile(FileName(run=run,date=date,tstamp=tstamp).timeAdjustments())
    #obsFile.loadHotPixCalFile(FileName(run=run,date=date,tstamp=tstamp).timeMask())
    #obsFile.loadCentroidListFile(FileName(run=run,date=date,tstamp=tstamp).centroidList())
    
    #And write out the results....
    obsFile.writePhotonList(outputFileName,firstSec,integrationTime,
                            pixRemapFileName=pixRemapFileName)
    
    #Read the results back in....
    #photFile = photList.PhotFile(outputFilename)
    
    print "Wrote photon list file: ", outputFileName

    
if __name__ == "__main__":
    
    import sys
    #hotPixelsTest(sys.argv[1])
    #hotPixelsTest2(startTime=2.5, endTime=6.3)
    testWritePhotonList()
    
    
