import os.path
import sys
import pickle
import tables
import numpy as np
import matplotlib.pyplot as mpl
from util import utils
from hotpix import hotPixels as hp


def quantifyBadTime(inputFileName, startTime=0, endTime=-1, 
                    defaultTimeMaskFileName='./testTimeMask.h5',
                    timeStep=1, fwhm=3.0, boxSize=5, nSigmaHot=3.0,
                    nSigmaCold=2.5,maxIter=5,binWidth=3,
                    dispToPickle=False, showHist=False):
    '''
    Function to calculate various metrics for the degree of bad pixel behaviour in a raw 
    raw obs file. Calculates the mean total hot/cold/dead time per good pixel (i.e., per 
    pixel which is not permanently dead, hot, or cold).
    
    Makes a couple of heat maps showing time spent bad in each way for each pixel,
    as well as a histogram of times spent bad for the *temporarily* bad pixels.
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
    
        binWidth - width of time bins for plotting the bad-time histogram (seconds)

        dispToPickle - if not False, saves the data for the histogram plot to a pickle file.
                       If a string, then uses that as the name for the pickle file. Otherwise
                       saves to a default name. Saves a dictionary with four entries, the first
                       three of which are each a flat array of total times spent bad for every
                       pixel (in sec) ('hotTime','coldTime','deadTime'). The fourth, 'duration',
                       is the duration of the input time-mask in seconds.
        
        showHist - if True, plot up a histogram of everything. Currently fails though 
                   if there are no bad-flagged intervals in any of the type categories
                   (or their intersections, hot+cold, cold+dead).
    
    
        Parameters passed on to findHotPixels routine if called (see also documentation
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
    
        maxIter         #Max num. of iterations for the bad pixel detection algorithm.
        


    OUTPUTS:
        
        A bunch of statistics on the different kinds of bad pixel behaviour, and a 
        'heat' plot showing how much time each pixel was bad in the array, for each
        type of behaviour. In theory it can plot a histogram of everything too, but
        currently it breaks if there are no bad intervals within any of the type
        categories (hot only, hot and cold, cold only, cold and dead, dead only...)
    
    '''
    
    defaultPklFileName = 'badPixTimeHist.pickle'
    
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
                                 display=True, dispToPickle=dispToPickle,
                                 weighted=False, maxIter=maxIter)
    
    
    #Read in the time mask file and calculate hot, cold, and 'other' bad time per pixel.
    timeMask = hp.readHotPixels(timeMaskFileName)
    hotTime = np.zeros((timeMask['nRow'],timeMask['nCol']))
    coldTime = np.zeros((timeMask['nRow'],timeMask['nCol']))
    deadTime = np.zeros((timeMask['nRow'],timeMask['nCol']))
    otherTime = np.zeros((timeMask['nRow'],timeMask['nCol']))
    hotIntervals = np.array([])
    coldIntervals = np.array([])
    deadIntervals = np.array([])
    otherIntervals = np.array([])

    reasonStringMap = timeMask['reasonEnum']
    for iRow in range(timeMask['nRow']):
        for iCol in range(timeMask['nCol']):
            for eachInterval,eachReasonCode in zip(timeMask['intervals'][iRow,iCol], timeMask['reasons'][iRow,iCol]):
                eachReasonString = reasonStringMap(eachReasonCode)   #Convert integer code to human readable string
                intSize = utils.intervalSize(eachInterval)
                if eachReasonString == 'hot pixel':
                    hotTime[iRow,iCol] += intSize
                    hotIntervals = np.append(hotIntervals, intSize)
                elif eachReasonString == 'cold pixel':
                    coldTime[iRow,iCol] += intSize
                    coldIntervals = np.append(coldIntervals, intSize)
                elif eachReasonString == 'dead pixel':
                    deadTime[iRow,iCol] += intSize
                    deadIntervals = np.append(deadIntervals, intSize)
                else:
                    otherTime[iRow,iCol] += intSize
                    otherIntervals = np.append(otherIntervals, intSize)

    
    if np.size(hotIntervals) == 0: hotIntervals = np.array([-1])
    if np.size(coldIntervals) == 0: coldIntervals = np.array([-1])
    if np.size(deadIntervals) == 0: deadIntervals = np.array([-1])
    if np.size(otherIntervals) == 0: otherIntervals = np.array([-1])
        
    totBadTime = hotTime+coldTime+deadTime+otherTime
    
    maskDuration = timeMask['endTime']-timeMask['startTime']
    
    #Figure out which pixels are hot, cold, permanently hot, temporarily cold, etc. etc.
    nPix = timeMask['nRow'] * timeMask['nCol']
    hotPix = hotTime > 0.1
    coldPix = coldTime > 0.1
    deadPix = deadTime > 0.1
    otherPix = otherTime > 0.1
    multiBehaviourPix = ( (np.array(hotPix,dtype=int)+np.array(coldPix,dtype=int)
                         +np.array(deadPix,dtype=int)+np.array(otherPix,dtype=int)) > 1)
    #assert np.all(deadTime[deadPix] == maskDuration)      #All dead pixels should be permanently dead....
    
    tol = timeStep/10. #Tolerance for the comparison operators below.
    
    permHotPix = hotTime >= maskDuration-tol
    permColdPix = coldTime >= maskDuration-tol
    permDeadPix = deadTime >= maskDuration-tol
    permOtherPix = otherTime >= maskDuration-tol
    permGoodPix = (hotTime+coldTime+deadTime+otherTime < tol)
    permBadPix = permHotPix | permColdPix | permDeadPix | permOtherPix
    
    tempHotPix = (hotTime < maskDuration-tol) & (hotTime > tol)
    tempColdPix = (coldTime < maskDuration-tol) & (coldTime > tol)
    tempDeadPix = (deadTime < maskDuration-tol) & (deadTime > tol)
    tempOtherPix = (otherTime < maskDuration-tol) & (otherTime > tol)
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
    print '% hot pixels: ', float(np.sum(permHotPix+tempHotPix))/nPix * 100.
    print '% cold pixels: ', float(np.sum(permColdPix+tempColdPix))/nPix * 100.
    print '% dead pixels: ', float(np.sum(permDeadPix+tempDeadPix))/nPix * 100.
    print '% other pixels: ', float(np.sum(permOtherPix+tempOtherPix))/nPix * 100.
    print
    print '% permanently hot pixels: ', float(np.sum(permHotPix))/nPix * 100.
    print '% permanently cold pixels: ', float(np.sum(permColdPix))/nPix * 100.
    print '% permanently dead pixels: ', float(np.sum(permDeadPix))/nPix * 100.
    print '% permanently "other" bad pixels: ', float(np.sum(permOtherPix))/nPix * 100.
    print
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
    if np.sum(tempOtherPix) > 0 or np.sum(permOtherPix) > 0:
        print '--------------------------------------------------------'
        print 'WARNING: Pixels flagged for "other" reasons detected - '
        print 'Histogram plot will not account for these!!'
        print '--------------------------------------------------------'

    #Display contour plots of the array of total bad times for each pixel for each type of behaviour
    utils.plotArray(hotTime, plotTitle='Hot time per pixel (sec)', fignum=None, cbar=True)
    utils.plotArray(coldTime, plotTitle='Cold time per pixel (sec)', fignum=None, cbar=True)
    utils.plotArray(otherTime, plotTitle='Other bad time per pixel (sec)', fignum=None, cbar=True)
    utils.plotArray(deadTime, plotTitle='Dead time per pixel (sec)', fignum=None, cbar=True)
    
    #Make histogram of time spent 'bad' for the temporarily bad pixels.
    #Ignore pixels flagged as bad for 'other' reasons (other than hot/cold/dead),
    #of which there should be none at the moment.
    assert np.all(otherPix == False)
    #Find total bad time for pixels which go only one of hot, cold, or dead
    onlyHotBadTime = totBadTime[hotPix & ~coldPix & ~deadPix]
    onlyColdBadTime = totBadTime[~hotPix & coldPix & ~deadPix]
    onlyDeadBadTime = totBadTime[~hotPix & ~coldPix & deadPix]
    #Find total bad time for pixels which alternate between more than one bad state
    hotAndColdBadTime = totBadTime[hotPix & coldPix & ~deadPix]
    hotAndDeadBadTime = totBadTime[hotPix & ~coldPix & deadPix]
    coldAndDeadBadTime = totBadTime[~hotPix & coldPix & deadPix]
    hotAndColdAndDeadBadTime = totBadTime[hotPix & coldPix & deadPix]
    allGoodBadTime = totBadTime[~hotPix & ~coldPix & ~deadPix]
    assert np.sum(allGoodBadTime) == 0

    if dispToPickle is not False:
        #Save to pickle file to feed into a separate plotting script, primarily for 
        #the pipeline paper.
        if type(dispToPickle) is str:
            pklFileName = dispToPickle
        else:
            pklFileName = defaultPklFileName
        pDict = {'hotTime':hotTime,
                 'coldTime':coldTime,
                 'deadTime':deadTime,
                 'onlyHotBadTime':onlyHotBadTime,
                 'onlyColdBadTime':onlyColdBadTime, 
                 'onlyDeadBadTime':onlyDeadBadTime,
                 'hotAndColdBadTime':hotAndColdBadTime,
                 'hotAndDeadBadTime':hotAndDeadBadTime,
                 'coldAndDeadBadTime':coldAndDeadBadTime,
                 'hotAndColdAndDeadBadTime':hotAndColdAndDeadBadTime,
                 'maskDuration':maskDuration}
        #pDict = {"hotTime":hotTime.ravel(),"coldTime":coldTime.ravel(),"deadTime":deadTime.ravel(),
        #         "duration":maskDuration}
        print 'Saving to file: ',pklFileName
        output = open(pklFileName, 'wb')
        pickle.dump(pDict,output)
        output.close()

    assert np.size(hotTime)==nPix and np.size(coldTime)==nPix and np.size(deadTime)==nPix
    assert (len(onlyHotBadTime)+len(onlyColdBadTime)+len(onlyDeadBadTime)+len(hotAndColdBadTime)
            +len(coldAndDeadBadTime)+len(hotAndDeadBadTime)+len(hotAndColdAndDeadBadTime)
            +len(allGoodBadTime))==nPix
       
    if showHist is True: 
        #Be sure it's okay to leave hot+dead pixels out, and hot+cold+dead pixels.
        #assert len(hotAndDeadBadTime)==0 and len(hotAndColdAndDeadBadTime)==0 
        mpl.figure()
        norm = 1. #1./np.size(hotTime)*100.
        dataList = [onlyHotBadTime,hotAndColdBadTime,onlyColdBadTime,coldAndDeadBadTime,onlyDeadBadTime]
        dataList2 = [x if np.size(x)>0 else np.array([-1.0]) for x in dataList]     #-1 is a dummy value for empty arrays, so that pyplot.hist doesn't barf.
        weights = [np.ones_like(x)*norm if np.size(x)>0 else np.array([0]) for x in dataList]
        mpl.hist(dataList2, range=None, #[-0.1,maskDuration+0.001],         #Eliminate data at 0sec and maskDuration sec. (always good or permanently bad)
                 weights=weights, 
                 label=['Hot only','Hot/cold','Cold only','Cold/dead','Dead only'], #,'Hot/dead','Hot/cold/dead'],
                 color=['red','pink','lightgray','lightblue','blue'], #,'green','black'],
                 bins=maskDuration/binWidth,histtype='stepfilled',stacked=True,log=False)         
        mpl.title('Duration of Bad Pixel Behaviour - '+os.path.basename(inputFileName))
        mpl.xlabel('Total time "bad" (sec)')
        mpl.ylabel('Percantage of pixels')
        mpl.legend()
        
        mpl.figure()
        mpl.hist(hotIntervals, bins=maskDuration/binWidth)
        print 'Median hot interval: ', np.median(hotIntervals)
        mpl.xlabel('Duration of hot intervals')
        mpl.ylabel('Number of intervals')
 
        mpl.figure()
        mpl.hist(coldIntervals, bins=maskDuration/binWidth)
        print 'Median cold interval: ', np.median(coldIntervals)
        mpl.xlabel('Duration of cold intervals')
        mpl.ylabel('Number of intervals')
 
        mpl.figure()
        mpl.hist(deadIntervals, bins=maskDuration/binWidth)
        print 'Median dead interval: ', np.median(deadIntervals)
        mpl.xlabel('Duration of dead intervals')
        mpl.ylabel('Number of intervals')
 
        mpl.figure()
        mpl.hist(otherIntervals, bins=maskDuration/binWidth)
        print 'Median "other" interval: ', np.median(otherIntervals)
        mpl.xlabel('Duration of "other" intervals')
        mpl.ylabel('Number of intervals')
        
        
 
    print 'Mask duration (s): ',maskDuration
    print 'Number of pixels: ',nPix
    print 'Fraction at 0s (hot,cold,dead): ', np.array([np.sum(hotTime<tol),np.sum(coldTime<tol),
                                                        np.sum(deadTime<tol)]) / float(nPix)                            
    print 'Fraction at '+str(maskDuration)+'s (hot,cold,dead): ', np.array([np.sum(hotTime>maskDuration-tol),
                                                                        np.sum(coldTime>maskDuration-tol),
                                                                        np.sum(deadTime>maskDuration-tol)])/float(nPix)
 
 
if __name__ == "__main__":
    '''
    To call from command line:
    
        python quantifyHotTime inputfilename [startTime endTime]
    
    '''
    
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
       
       