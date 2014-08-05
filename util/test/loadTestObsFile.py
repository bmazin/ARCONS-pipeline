import os.path
from util import FileName as fn
from util import ObsFile as of

def loadTestObsFile(obsFileName=None, wvlCalFileName=None, flatCalFileName=None,
                    fluxCalFileName=None, hotPixFileName=None, timeAdjFileName=None):
    '''
    Quick function to load up a test obs file and all its calibrations, and return
    the resulting ObsFile object.
    INPUTS:
        obsFileName, wvlCalFileName, flatCalFileName, fluxCalFileName, hotPixFileName - 
            set any of these to override the default values (see below).
    '''
    
    run = 'PAL2012'
    date = '20121208'
    if obsFileName is None: obsFileName = fn.FileName(run=run,date=date,tstamp='20121209-044636').obs()   #'obs_20121209-044636.h5'
    if timeAdjFileName is None: timeAdjFileName = fn.FileName(run='PAL2012').timeAdjustments()
    if wvlCalFileName is None: wvlCalFileName = fn.FileName(run=run, date=date, tstamp='20121209-060704').calSoln()    #'calsol_20121209-060704.h5'
    if flatCalFileName is None: flatCalFileName = fn.FileName(run=run, date='20121210').flatSoln()    #'flatsol_20121210.h5'
    if fluxCalFileName is None: fluxCalFileName = fn.FileName(run=run, date=date, tstamp='20121209-020416').fluxSoln()
    #if hotPixFileName is None: hotPixFileName = fn.FileName(run=run, date=date, tstamp='20121209-044636')
    if hotPixFileName is None: hotPixFileName = os.path.abspath('./test-calibratedHotPix_20121209-044636.h5')
    paramFile = os.path.join(os.path.dirname(__file__),'../../params/hotPixels.dict')
      
    print 'Loading obs file and calibrations:'
    print obsFileName
    obsFile = of.ObsFile(obsFileName)
    obsFile.loadTimeAdjustmentFile(timeAdjFileName)
    print timeAdjFileName
    obsFile.loadBestWvlCalFile()
    print obsFile.wvlCalFileName
    print flatCalFileName
    obsFile.loadFlatCalFile(flatCalFileName)
    print fluxCalFileName
    obsFile.loadFluxCalFile(fluxCalFileName)
    print hotPixFileName
    obsFile.loadHotPixCalFile(hotPixFileName)
    print 'Setting wavelength cutoffs (default)'
    obsFile.setWvlCutoffs()

    return obsFile