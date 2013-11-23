import os.path
import sys
import tables
import numpy as np
from util import utils
from hotpix import hotPixels as hp


def quantifyHotTime(inputFileName, startTime=0, endTime=-1, 
                    defaultTimeMaskFileName='./testTimeMask.hd5',
                    timeStep=1, fwhm=3.0, boxSize=5,
                    nSigmaHot=2.5, nSigmaCold=2.5,
                    maxIter=5):
    '''
    Function to calculate a metric for the degree of hot pixel behaviour in a raw 
    raw obs file. Calculates the mean total 'hot time' per good pixel (i.e., per 
    pixel which is not permanently dead, hot, or cold).
    JvE Nov 20 2013.
    
    The parameters for the finding algorithm may need to be tuned a little,
    but the defaults should basically work, in principle. nSigmaHot and nSigmaCold
    are good places to start if things need playing around with.
    
    e.g. call, in principle:
    
        from hotpix import quantifyHotTime as qht
        qht.quantifyHotTime('/Users/vaneyken/Data/UCSB/ARCONS/turkDataCopy/ScienceData/PAL2012/20121208/obs_20121209-120530.h5')
    
    - should be all it needs....
    
    
    INPUTS:
    
        inputFileName - either a raw obs. file or a time mask file.
        
        startTime, endTime - start and end times within the obs file to
                        calculate the hot pixels for. (endTime =-1 means to end of file).
        
        defaultTimeMaskFileName - use this filename to output new time mask to
                                (if inputFileName is an obs file)
    
        Parameters passed on to findHotPixels routine (see also documentation
        for that function):
    
        timeStep        #Check for hot pixels every timeStep seconds
    
        fwhm            #Expected full width half max of PSF in pixels. Any pixel
                        #with flux much tighter than this will be flagged as bad.
                        #Larger value => more sensitive hot pixel flagging.
    
        boxSize         #Compare flux in each pixel with median flux in a 
                        #surrounding box of this size on a side.
    
        nSigmaHot       #Require flux in a pixel to be > nSigmaHot std. deviations
                        #above the max expected for a Gaussian PSF in order for it
                        #to be flagged as hot. Larger value => less sensitive flagging.
    
        nSigmaCold      #Require flux to be < nSigmaCold std. deviations below the median 
                        #in a surrounding box in order to be flagged as cold (where std.
                        #deviation is estimated as the square root of the median flux).
    
        maxIter=5       #Max num. of iterations for the bad pixel detection algorithm.


    OUTPUTS:
        
        A bunch of statistics on the different kinds of bad pixel behaviour, and a 
        'heat' plot showing how much time each pixel was bad in the array, for each
        type of behaviour. 
    

    '''
    
    #Check whether the input file is a time mask or a regular obs file.
    absInputFileName = os.path.abspath(inputFileName)   #To avoid weird issues with the way findHotPixels expands paths....
    hdffile = tables.openFile(absInputFileName)
    inputIsTimeMaskFile = '/timeMasks' in hdffile   
    hdffile.close()
    
    #Decide whether to generate a new time mask file or not
    if inputIsTimeMaskFile:
        print 'Input file is a time mask file'
        timeMaskFileName = absInputFileName
    else:
        print 'Assuming input file is an obs. file'
        timeMaskFileName = os.path.abspath(defaultTimeMaskFileName)
        if os.path.exists(timeMaskFileName):
            response=''
            while response != 'u' and response !='r':
                response = raw_input(timeMaskFileName+' already exists - "u" to use this (default) or "r" to regenerate? ')
                response = response.strip().lower()
                if response == '': response = 'u'
        else:
            response = 'r'         #If the file *didn't already exist, pretend the user entered 'r' despite not having been asked.
            
        if response == 'r':
                #Make/regenerate the hot pixel file.
                print 'Making new hot pixel time mask file '+timeMaskFileName
                hp.findHotPixels(inputFileName=absInputFileName,
                                 outputFileName=timeMaskFileName,
                                 startTime=startTime,
                                 endTime=endTime,
                                 timeStep=timeStep, fwhm=fwhm,
                                 boxSize=boxSize, nSigmaHot=nSigmaHot, nSigmaCold=nSigmaCold,
                                 display=True,
                                 weighted=False, maxIter=maxIter)
    
    
    #Read in the time mask file and calculate hot, cold, and 'other' bad time per pixel.
    timeMask = hp.readHotPixels(timeMaskFileName)
    hotTime = np.zeros((timeMask['nRow'],timeMask['nCol']))
    coldTime = np.zeros((timeMask['nRow'],timeMask['nCol']))
    deadTime = np.zeros((timeMask['nRow'],timeMask['nCol']))
    otherTime = np.zeros((timeMask['nRow'],timeMask['nCol']))
    reasonStringMap = timeMask['reasonEnum']
    for iRow in range(timeMask['nRow']):
        for iCol in range(timeMask['nCol']):
            for eachInterval,eachReasonCode in zip(timeMask['intervals'][iRow,iCol], timeMask['reasons'][iRow,iCol]):
                eachReasonString = reasonStringMap(eachReasonCode)   #Convert integer code to human readable string
                if eachReasonString == 'hot pixel':
                    hotTime[iRow,iCol] += utils.intervalSize(eachInterval)
                elif eachReasonString == 'cold pixel':
                    coldTime[iRow,iCol] += utils.intervalSize(eachInterval)
                elif eachReasonString == 'dead pixel':
                    deadTime[iRow,iCol] += utils.intervalSize(eachInterval)
                else:
                    otherTime[iRow,iCol] += utils.intervalSize(eachInterval)
    
    maskDuration = timeMask['endTime']-timeMask['startTime']
    
    #Figure out which pixels are hot, cold, permanently hot, temporarily cold, etc. etc.
    nPix = timeMask['nRow'] * timeMask['nCol']
    hotPix = hotTime > 0
    coldPix = coldTime > 0
    deadPix = deadTime > 0
    otherPix = otherTime > 0
    multiBehaviourPix = (hotPix+coldPix+deadPix+otherPix > 1)
    #assert np.all(deadTime[deadPix] == maskDuration)      #All dead pixels should be permanently dead....
    
    permHotPix = hotTime >= maskDuration
    permColdPix = coldTime >= maskDuration
    permDeadPix = deadTime >= maskDuration
    permOtherPix = otherTime >= maskDuration
    permGoodPix = (hotTime+coldTime+deadTime+otherTime == 0)
    permBadPix = permHotPix | permColdPix | permDeadPix | permOtherPix
    
    tempHotPix = (hotTime < maskDuration) & (hotTime > 0)
    tempColdPix = (coldTime < maskDuration) & (coldTime > 0)
    tempDeadPix = (deadTime < maskDuration) & (deadTime > 0)
    tempOtherPix = (otherTime < maskDuration) & (otherTime > 0)
    tempGoodPix = tempHotPix | tempColdPix | tempDeadPix | tempOtherPix     #Bitwise or should work okay with boolean arrays
    tempBadPix = tempGoodPix        #Just to be explicit about it....
    
    nGoodPix = np.sum(permGoodPix | tempGoodPix)        #A 'good pixel' is either temporarily or permanently good
    
    #assert np.sum(tempDeadPix) == 0             #Shouldn't be any temporarily dead pixels. APPARENTLY THERE ARE....
    assert np.all((tempHotPix & permHotPix)==False)     #A pixel can't be permanently AND temporarily hot
    assert np.all((tempColdPix & permColdPix)==False)   #... etc.
    assert np.all((tempOtherPix & permOtherPix)==False)
    assert np.all((tempDeadPix & permDeadPix)==False)

    
    #Print out the results
    print '----------------------------------------------'
    print
    print '# pixels total: ', nPix
    print 'Mask duration (sec): ', maskDuration
    print
    print '% permanently hot pixels: ', float(np.sum(permHotPix))/nPix * 100.
    print '% permanently cold pixels: ', float(np.sum(permColdPix))/nPix * 100.
    print '% permanently dead pixels: ', float(np.sum(permDeadPix))/nPix * 100.
    print '% permanently "other" bad pixels: ', float(np.sum(permOtherPix))/nPix * 100.
    print '% temporarily hot pixels: ', float(np.sum(tempHotPix))/nPix * 100.
    print '% temporarily cold pixels: ', float(np.sum(tempColdPix))/nPix * 100.
    print '% temporarily dead pixels: ', float(np.sum(tempDeadPix))/nPix * 100.
    print '% temporarily "other" bad pixels: ', float(np.sum(tempOtherPix))/nPix * 100.
    print
    print '% pixels showing multiple bad behaviours: ', float(np.sum(multiBehaviourPix))/nPix * 100.
    print
    print '% permanently bad pixels: ', float(np.sum(permBadPix))/nPix * 100.
    print '% temporarily bad pixels: ', float(np.sum(tempGoodPix))/nPix * 100.      #Temp. good == temp. bad!
    print '% permanently good pixels: ', float(np.sum(permGoodPix))/nPix * 100.
    print
    print 'Mean temp. hot pixel time per good pixel: ', np.sum(hotTime[tempHotPix])/nGoodPix
    print 'Mean temp. hot pixel time per temp. hot pixel: ', np.sum(hotTime[tempHotPix])/np.sum(tempHotPix)
    print
    print 'Mean temp. cold pixel time per good pixel: ', np.sum(coldTime[tempColdPix])/nGoodPix
    print 'Mean temp. cold pixel time per temp. cold pixel: ', np.sum(coldTime[tempColdPix])/np.sum(tempColdPix)
    print
    print 'Mean temp. "other" bad pixel time per good pixel: ', np.sum(otherTime[tempOtherPix])/nGoodPix
    print 'Mean temp. "other" bad pixel time per temp. "other" pixel: ', np.sum(otherTime[tempOtherPix])/np.sum(tempOtherPix)
    print
    print '(All times in seconds)'
    print
    print 'Done.'
    print
    

    #Display contour plots of the array of total bad times for each pixel for each type of behaviour
    utils.plotArray(hotTime, plotTitle='Hot time per pixel (sec)', fignum=None, cbar=True)
    utils.plotArray(coldTime, plotTitle='Cold time per pixel (sec)', fignum=None, cbar=True)
    utils.plotArray(otherTime, plotTitle='Other bad time per pixel (sec)', fignum=None, cbar=True)
    utils.plotArray(deadTime, plotTitle='Dead time per pixel (sec)', fignum=None, cbar=True)
    
    
 
 
 
 
if __name__ == "__main__":
    
    inputFileName = None
    startTime = 0
    endTime = -1
    nArg = len(sys.argv)
    if nArg == 1: print 'Need to provide input file name'
    if nArg > 1: inputFileName = sys.argv[1]
    if nArg > 2: startTime = float(sys.argv[2])
    if nArg > 3: endTime = float(sys.argv[3])
    if nArg > 4: print 'Too many arguments provided - ignoring everything after third'
    if inputFileName is not None: quantifyHotTime(inputFileName, startTime=startTime, endTime=endTime)
       
       