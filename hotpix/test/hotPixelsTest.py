import hotpix.hotPixels as hp
import numpy as np
from interval import interval
import sys

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

    hp.findHotPixels(paramFile, testFileName, outputFile, timeStep=timeStep,
                     startTime=testStartTime, endTime=testEndTime,
                     fwhm=fwhm, boxSize=boxSize, nSigma=nSigma, display=True)
    
    intermediateOutput = hp.checkInterval(inputFileName=testFileName, display=True,
                                          firstSec=testStartTime,
                                          intTime=testEndTime - testStartTime,
                                          fwhm=fwhm, boxSize=boxSize,
                                          nSigma=nSigma)
    
    hpOutput = hp.readHotPixels('testoutput.h5')

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
    print "All seems okay."




if __name__ == "__main__":
    
    hotPixelsTest(sys.argv[1])
    
    