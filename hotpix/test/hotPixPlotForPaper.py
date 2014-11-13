'''
To produce .pkl files for the plots demonstrating hot pixel masking for the pipeline paper.
JvE 11/7/2014
'''

import hotpix.hotPixels as hp
from util import FileName as fn
from util import ObsFile

def hotPixPlot():
    
    if 1==1:
    
        #File #1 - 2012, SDSS J0926
        run = 'PAL2012'
        date = '20121208'
        inputFile = fn.FileName(run=run,date=date,tstamp='20121209-120530').obs()
        timeAdjFileName = fn.FileName(run='PAL2012').timeAdjustments()
        #wvlCalFileName = fn.FileName(run=run, date=date, tstamp='20121209-060704').calSoln()    #'calsol_20121209-060704.h5'
        flatCalFileName = fn.FileName(run=run, date='20121210').flatSoln()    #'flatsol_20121210.h5'
        fluxCalFileName = fn.FileName(run=run, date='20121211', tstamp='absolute_021727').fluxSoln()
        
        obsFile = ObsFile.ObsFile(inputFile)
        print inputFile
        obsFile.loadTimeAdjustmentFile(timeAdjFileName)
        print timeAdjFileName
        obsFile.loadBestWvlCalFile()
        print obsFile.wvlCalFileName
        print flatCalFileName
        obsFile.loadFlatCalFile(flatCalFileName)
        print fluxCalFileName
        obsFile.loadFluxCalFile(fluxCalFileName)
        obsFile.setWvlCutoffs()   #Use default wavelength cutoffs
        
        hp.findHotPixels(obsFile=obsFile, outputFileName='hpTest.h5', startTime=54,
                        endTime=55, display=True, dispMinPerc=0, dispMaxPerc=98, maxIter=10,
                        nSigmaHot=4.0, bkgdPercentile=50, fwhm=2.0, diagnosticPlots=False,
                        boxSize=5,weighted=True,fluxWeighted=True, useRawCounts=True,
                        dispToPickle=True)
    
    
    #File 2
    run = 'PAL2014'
    date = '20141021'  #'20141022' 
    tstamp = '20141022-092629'               # '20141023-050637'
    inputFile = fn.FileName(run=run,date=date,tstamp=tstamp).obs()
    #Time adjustment file not needed
    flatCalFileName = fn.FileName(run=run, date=date).flatSoln()
    #Ignore flux cal for this one (shouldn't make much difference, and don't know if new one is set up yet)
    #fluxCalFileName = fn.FileName(run=run, date=', tstamp='absolute_021727').fluxSoln()
    
    obsFile = ObsFile.ObsFile(inputFile)
    print inputFile
    obsFile.loadBestWvlCalFile()
    print obsFile.wvlCalFileName
    print flatCalFileName
    obsFile.loadFlatCalFile(flatCalFileName)
    obsFile.setWvlCutoffs()   #Use default wavelength cutoffs
    
    hp.findHotPixels(obsFile=obsFile, outputFileName='hpTest.h5', startTime=24,
                    endTime=25, display=True, dispMinPerc=0, dispMaxPerc=99, maxIter=10,
                    nSigmaHot=4.0, bkgdPercentile=50, fwhm=2.0, diagnosticPlots=False,
                    boxSize=5,weighted=True,fluxWeighted=False, useRawCounts=True,
                    dispToPickle=True)
    
    return obsFile
