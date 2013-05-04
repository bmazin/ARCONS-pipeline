from util.FileName import FileName
from util.ObsFile import ObsFile


def TestWritePhotonList(outputFilename=None):
    '''
    Test run of obsFile.writePhotonList. fileName can be used
    to specify the output file name. If not specified, default
    name/location is used.
    '''

    #Details of example obs file to run test on.
    run = 'PAL2012'
    date = '20121207'
    tstamp = '20121208-074649'
    calTstamp='20121208-070505'
    fluxTstamp='20121208-133002'
    
    #Load up the obs file
    obsFileName = FileName(run=run, date=date, tstamp=tstamp)
    obsFile = ObsFile(obsFileName.obs())
    
    #Load up associated calibrations
    obsFile.loadWvlCalFile(FileName(run=run,date=date,tstamp=calTstamp).calSoln())
    obsFile.loadFlatCalFile(FileName(run=run,date=date).flatSoln())
    obsFile.loadFluxCalFile(FileName(run=run,date=date,tstamp=fluxTstamp).fluxSoln())
    obsFile.loadTimeAdjustmentFile(FileName(run=run,date=date,tstamp=tstamp).timeAdjustments())
    obsFile.loadHotPixCalFile(FileName(run=run,date=date,tstamp=tstamp).timeMask())
    
    #And write out the results....
    obsFile.writePhotonList(outputFilename)
    
    

    
if __name__ == "__main__":
    
    import sys
    #hotPixelsTest(sys.argv[1])
    #hotPixelsTest2(startTime=2.5, endTime=6.3)
    TestWritePhotonList()
    
    
