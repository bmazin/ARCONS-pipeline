import hotpix.hotPixels as hp
import numpy as np
from interval import interval
import sys
import util.ObsFile as of
import matplotlib.pylab as mpl

'''Internal consistency checks for hotPixels.py'''


def hotPixelsTest(testFileName='obs_20121211-024511.h5'):
    '''
    Runs some basic checks for consistency between intermediate output
    masks and the final 'bad time lists'.
    
    To run
        - from Python:
            hotPixelsTest('someObsFile.h5')
        
        - from command line:
            python hotPixelsTest.py someObsFile.h5
    
    
    (No parameter file need be supplied).
    
    ''' 
    
    outputFile = 'testoutput.h5'
    paramFile = '../../params/hotPixels.dict'
    testStartTime = 2   #In seconds
    testEndTime = 4     #In seconds
    timeStep = 2        #In seconds
    fwhm = 3.0
    boxSize = 5
    nSigma = 3.0

    hp.findHotPixels(paramFile=paramFile, inputFileName=testFileName,
                     outputFileName=outputFile, timeStep=timeStep,
                     startTime=testStartTime, endTime=testEndTime,
                     fwhm=fwhm, boxSize=boxSize, nSigma=nSigma, display=True)
    
    intermediateOutput = hp.checkInterval(inputFileName=testFileName, display=True,
                                          firstSec=testStartTime,
                                          intTime=testEndTime - testStartTime,
                                          fwhm=fwhm, boxSize=boxSize,
                                          nSigma=nSigma)
    
    hpOutput = hp.readHotPixels(outputFile)

    intMask = intermediateOutput['mask']
    intervals = hpOutput['intervals']
    reasons = hpOutput['reasons']

    #Find the number of entries for each pixel in both the 'intervals' and the
    #'reasons' arrays.
    nIntervals = np.reshape([len(x) for x in intervals.flat], np.shape(intervals))
    nReasons = np.reshape([len(x) for x in reasons.flat], np.shape(reasons))

    #Union the interval lists for each pixel to give an array of single (multi-component) interval objects:
    uIntervals = np.reshape(np.array([interval.union(x) for x in intervals.flat],
                                     dtype='object'), np.shape(intervals))


    #Create a boolean mask that should be True for all hot pixels within the test time range
    finalMask = np.reshape([(interval(testStartTime, testEndTime) in x) 
                            for x in uIntervals.flat], np.shape(uIntervals))   

    assert np.all(np.equal(intMask, finalMask))
    assert np.all(np.equal(nIntervals, nReasons))

    print
    print "All seems okay. The two plots shown should look identical."



def hotPixelsTest2(startTime=12.3, endTime=23.1, getRawCount=False):
    '''
    Runs a check of the bad pixel time masking. Checks:
        - Number of photons removed from returned image is consistent
        with time masks.
        - That timestamps of all photons in the time-masked image are outside
        the time intervals defined for each pixel in the hot pixel file.
        
    INPUTS:
        startTime - time from beginning of obs file at which to start the check.
        endTime - time from beginning to obs file at which to end (both in seconds).
        getRawCounts - if True, use raw, non-wavelength calibrated photon counts
                        with no wavelength cutoffs applied.
    '''
    
    dir = '/Users/vaneyken/Data/UCSB/ARCONS/Palomar2012/hotPixTest2/'
    obsFileName = 'obs_20121209-044636.h5'
    wvlCalFileName = 'calsol_20121209-060704.h5'
    flatCalFileName = 'flatsol_20121210.h5'
    hotPixFileName = 'hotPix_20121209-044636.h5'
    startTime = float(startTime)    #Force these values to floats to make sure
    endTime = float(endTime)        #that getPixelCountImage calls getPixelSpectrum
                                    #and applies wavelength cutoffs consistently.
    
    intTime = endTime - startTime
    
    obsFile = of.ObsFile(dir + obsFileName)
    obsFile.loadWvlCalFile(dir + wvlCalFileName)
    obsFile.loadFlatCalFile(dir + flatCalFileName)
    print 'Loading hot pixel file into obsFile...'
    obsFile.loadHotPixCalFile(dir + hotPixFileName)
    obsFile.setWvlCutoffs()
    print 'Getting image with masking...'
    imhp = obsFile.getPixelCountImage(startTime, intTime, weighted=False, getRawCount=getRawCount)
    print 'Getting image without masking...'
    obsFile.switchOffHotPixTimeMask()
    im = obsFile.getPixelCountImage(startTime, intTime, weighted=False, getRawCount=getRawCount)
    
    diffim = im - imhp #Should end up containing the total number of photons masked from each pixel.

    
    print 'Displaying images...'
    mpl.ion()
    mpl.matshow(imhp)
    mpl.title('Hot-pixel masked')
    mpl.colorbar()
    mpl.matshow(im)
    mpl.title('Unmasked')
    mpl.colorbar()
    
    print 'Loading local version of hot pixel file for direct inspection...'
    hotPix = hp.readHotPixels(dir + hotPixFileName)
    
    if True:
        print 'Checking for consistency in number of photons removed....'
        for iRow in range(np.shape(diffim)[0]):
            for iCol in range(np.shape(diffim)[1]):
                nMaskedPhotons = 0
                for eachInter in hotPix['intervals'][iRow, iCol]:
                    #Make sure there is only one component per interval
                    assert len(eachInter) == 1
                    #If the interval overlaps with the time range in question then 
                    #add the number of photons in the interval to our running tally.
                    if eachInter[0][0] < endTime and eachInter[0][1] > startTime:
                        firstSec = max(eachInter[0][0], startTime)
                        lastSec = min(eachInter[0][1], endTime)
                        nMaskedPhotons += obsFile.getPixelCount(iRow, iCol, firstSec=firstSec,
                                                                integrationTime=lastSec - firstSec,
                                                                weighted=False, fluxWeighted=False,
                                                                getRawCount=getRawCount)
                assert nMaskedPhotons == diffim[iRow, iCol]
        print 'Okay.'
    
    print 'Checking timestamps of remaining photons for consistency with exposure masks'
    obsFile.switchOnHotPixTimeMask()       #Switch on hot pixel masking
    for iRow in range(np.shape(diffim)[0]):
        for iCol in range(np.shape(diffim)[1]):
            timeStamps, dummy1, dummy2 = obsFile.getTimedPacketList(iRow, iCol, firstSec=startTime, integrationTime=intTime)
            #timeStamps = timeStamps[np.logical_and(timeStamps<=endTime, timeStamps>=startTime)]
            badInterval = interval.union(hotPix['intervals'][iRow, iCol])
           
            #The following check would be nice, but doesn't work because getTimedPacketList doesn't do a wavelength cut like getPixelCount does.
            #assert len(timeStamps) == im[iRow,iCol]     #Double check that the number of photons returned matches the number in the masked image
            #
            
            for eachTimestamp in timeStamps:
               assert eachTimestamp not in badInterval   #Check that none of those photons' timestamps are in the masked time range

    print 'Okay.'
    print       
    print 'Done. All looks good.'



if __name__ == "__main__":
    
    #hotPixelsTest(sys.argv[1])
    hotPixelsTest2(startTime=2.5, endTime=6.3)
    
