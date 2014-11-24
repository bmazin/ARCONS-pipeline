'''
To produce .pkl files for the plots demonstrating hot pixel masking for the pipeline paper.
JvE 11/7/2014
'''

import matplotlib.pyplot as mpl
import hotpix.hotPixels as hp
from hotpix import quantifyBadTime as qbt
from util import FileName as fn
from util import ObsFile

def hotPixPlot():
    
    if 1==1:
    
        #File #1 - 2012, SDSS J0926
        run = 'PAL2012'
        date = '20121208'
        tstamp='20121209-120530'
        inputFile = fn.FileName(run=run,date=date,tstamp=tstamp).obs()
        timeMaskFile = fn.FileName(run=run,date=date,tstamp=tstamp).timeMask()
        timeAdjFileName = fn.FileName(run='PAL2012').timeAdjustments()
        #wvlCalFileName = fn.FileName(run=run, date=date, tstamp='20121209-060704').calSoln()    #'calsol_20121209-060704.h5'
        #flatCalFileName = fn.FileName(run=run, date='20121210').flatSoln()    #'flatsol_20121210.h5'
        #fluxCalFileName = fn.FileName(run=run, date='20121211', tstamp='absolute_021727').fluxSoln()
        
        obsFile = ObsFile.ObsFile(inputFile)
        print inputFile
        obsFile.loadTimeAdjustmentFile(timeAdjFileName)
        print timeAdjFileName
        obsFile.loadBestWvlCalFile()
        print obsFile.wvlCalFileName
        #print flatCalFileName
        #obsFile.loadFlatCalFile(flatCalFileName)
        #print fluxCalFileName
        #obsFile.loadFluxCalFile(fluxCalFileName)
        obsFile.setWvlCutoffs()   #Use default wavelength cutoffs
        
        hp.findHotPixels(obsFile=obsFile, outputFileName='badPix1.h5', startTime=54,
                        endTime=55, display=True, dispMinPerc=0, dispMaxPerc=98, maxIter=10,
                        nSigmaHot=4.0, bkgdPercentile=50, fwhm=2.5, diagnosticPlots=False,
                        boxSize=5, weighted=False,fluxWeighted=False, useRawCounts=True,
                        dispToPickle='hotPixData-'+tstamp+'.pkl')
        mpl.title(obsFile.getFromHeader('target'))
        
        #Get the stats for the same image
        qbt.quantifyBadTime(inputFile,startTime=0,endTime=-1,maxIter=10,nSigmaHot=4.0,
                            bkgdPercentile=50,fwhm=2.0,boxSize=5,useRawCounts=True,
                            defaultTimeMaskFileName=timeMaskFile)
        

    
    #File 2
    run = 'PAL2014'
    date = '20141021'  #'20141022' 
    tstamp = '20141022-092629'   #'20141022-032331'               # '20141023-050637'
    inputFile = fn.FileName(run=run,date=date,tstamp=tstamp).obs()
    timeMaskFile = fn.FileName(run=run,date=date,tstamp=tstamp).timeMask()

    #Time adjustment file not needed
    #flatCalFileName = fn.FileName(run=run, date=date).flatSoln()
    #Ignore flux cal for this one (shouldn't make much difference, and don't know if new one is set up yet)
    #fluxCalFileName = fn.FileName(run=run, date=', tstamp='absolute_021727').fluxSoln()
    
    obsFile = ObsFile.ObsFile(inputFile)
    print inputFile
    obsFile.loadBestWvlCalFile()
    print obsFile.wvlCalFileName
    #print flatCalFileName
    #obsFile.loadFlatCalFile(flatCalFileName)
    obsFile.setWvlCutoffs()   #Use default wavelength cutoffs
    
    hp.findHotPixels(obsFile=obsFile, outputFileName='badPix2.h5',
                     startTime=24, endTime=25, display=True, dispMinPerc=0, dispMaxPerc=99, maxIter=10,
                     nSigmaHot=4.0, bkgdPercentile=50, fwhm=2.0, diagnosticPlots=False,
                     boxSize=5, weighted=False,fluxWeighted=False, useRawCounts=True,
                     dispToPickle='hotPixData-'+tstamp+'.pkl')
        
    mpl.title(obsFile.getFromHeader('target'))
    
    #Get the stats for the same image
    qbt.quantifyBadTime(inputFile,startTime=0,endTime=-1,maxIter=10,nSigmaHot=4.0,
                        bkgdPercentile=50,fwhm=2.0,boxSize=5,useRawCounts=True,
                        defaultTimeMaskFileName=timeMaskFile)
        
    return obsFile
