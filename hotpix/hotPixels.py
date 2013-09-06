'''
Author: Julian van Eyken
Date: Jan 25 2013

Routines for checking for hot-pixels in (uncalibrated) obs. files. In principle
should be relatively easily extendable to include other checks (e.g. for 'cold'
pixels).

The supplied obs. file is stepped through in short time steps, and a 2D mask
made for each step. The results are then converted to lists of bad time periods
for each pixel. The current detection algorithm compares the flux in each pixel
with that in each pixel in a surrounding box. If the difference is significantly
higher than could be expected for a stellar (point-source) PSF (allowing for
noise), then the flux level cannot be real and the pixel is flagged for the 
time step under examination.


Main routines of interest are:



-------------
findHotPixels: currently the main routine. Takes an obs. file as input and 
                outputs an .h5 file containing lists of bad time intervals for
                all the pixels, along with reason indicators for each of those
                intervals (*Note currently only one reason - i.e. hot pixel - 
                is covered; but this should be easily extendible to cover e.g.
                cold pixels*).
                
readHotPixels: reads in a hot-pixel .h5 file into somewhat sensible structures
                for use by external routines.

checkInterval: creates a 2D mask for a given time interval within a given 
                exposure.
                
getEffIntTimeImage: Once hot pixels have been found, can use this to return the
                effective integration times for each pixel in a given period after
                masking out bad times. 
                
getHotPixels: Similar to getEffIntTimeImage, but just returns a boolean array 
                indicating which pixels went bad at any point during the specified
                period.
-------------



Dependencies: pytables; pyinterval; headers.TimeMask; util.ObsFile; numpy;
              matplotlib; util.readDict

To do:
    - Extend to check for 'cold' pixels
    - Suppress annoying invalid value warnings that come up (and can be ignored)
    - Speed up (most time is spent reading the h5 file; significant overhead can
        probably be saved by reading in the whole file at once instead of
        second-by-second).
    - At some point, perhaps update to allow sub-second time steps.
    - Time intervals are currently stored in units of seconds - need to convert
        to clock ticks for consistency with TimeInterval class. FIXED - 02/12/2013.
    - If necessary, apply algorithm only at the red end, where hot-pixel
        contrast should be higher.
    - Output files seem unnecessarily large. May be some way of making this more
        efficient.
    - Current image display in checkInterval is wrong way up - need to update
      to use Danica's image display routine.

See individual routines for more detail.

'''

import os.path
import sys
from math import *
from interval import interval
import tables
import numpy as np
import numpy.ma as ma
import matplotlib.pylab as mpl
import util.ObsFile as ObsFile
import util.utils as utils
import headers.TimeMask as tm
import util.readDict as readDict


headerGroupName = 'header'  #Define labels for the output .h5 file.
headerTableName = 'header'
dataGroupName = 'timeMasks' #data *table* name is constructed on a per-pixel basis (see constructDataTableName)
nRowColName = 'nRow'        #Name of the no.-of-rows column in the header for the output .h5 file.
nColColName = 'nCol'        #Ditto for no.-of-columns column.
obsFileColName = 'obsFileName' #Ditto for the column containing the name of the original obs. file.
ticksPerSecColName = 'ticksPerSec'  #Ditto for column containing number of clock ticks per second for original obs. file.
expTimeColName = 'expTime'  #Ditto for column containing total exposure time

def constructDataTableName(x, y):
    '''Construct the name for the pytables table for a given pixel at location x,y'''
    return 'x' + str(int(x)) + 'y' + str(int(y))



class headerDescription(tables.IsDescription):
    """Description of header info structure for the output hot-pixel .h5 file."""
    obsFileName = tables.StringCol(40)  #To record associated obs file name
    nCol = tables.UInt32Col(dflt=-1)   #To record how many pixel columns/rows were in
    nRow = tables.UInt32Col(dflt=-1)   #the data used to construct the mask.
    ticksPerSec = tables.Float64Col(dflt=np.nan)    #So that code which reads in the file can back out intervals in seconds if needed.
    exptime = tables.Float64Col(dflt=np.nan)        #So that per-pixel effective exposure times can be calculated any time without having to refer to the original obs file. Added 6/21/2013, JvE, although not yet used elsewhere in the code.
 
 
def checkInterval(firstSec=None, intTime=None, fwhm=4.0, boxSize=5, nSigma=3.0,
                  obsFile=None, inputFileName=None, image=None,
                  display=False, weighted=False, maxIter=5):
    '''
    To find the hot pixels in a given time interval for an observation file.
    Compares the ratio of flux in each pixel to the median of the flux in an
    enclosing box. If the ratio is too high -- i.e. the flux is too tightly 
    distributed compared to a Gaussian PSF of the expected FWHM -- then the 
    pixel is flagged as bad.
    
    Main return value is a 2D integer array of flags corresponding to 
    the input image (see OUTPUTS below).
    
    INPUTS:
        Must provide one of:-
            obsFile - an ObsFile instance for an already-open observation file.
            inputFileName: String - input observation file name if not already opened.
            image - a 2D image array of photon counts.
        
        Other:-
        
        firstSec: Scalar integer - start time in seconds within obs. file from 
                    which to integrate when looking for hot pixels.
        intTime: Scalar integer - integration time for hot pixel search 
                    (seconds).
        fwhm: Scalar float. Estimated full-width-half-max of the PSF (in 
                    pixels).
        boxSize: Scalar integer. Size of edge of box used for calculating median 
                    flux in the region surrounding each pixel.
        nSigma: Scalar float. If the flux ratio for a pixel is (nSigma x expected
                    error) above the max expected given the PSF FWHM, then flag it 
                    as hot.
        display: Boolean. If true, then display the input image and mark those 
                    flagged as hot with a dot.
        maxIter: Scalar integer. Maximum number of iterations allowed.
        

    OUTPUTS:
        A dictionary containing the result and various diagnostics. Keys are:
        
        'mask': the main output. Contains a 2D array of integers of the same
                shape as the input image, where:
                        0 = Good pixel
                        1 = Hot pixel
                        (Other numbers reserved for later use...)
        'image': 2D array containing the input image
        'medfiltimage': The median-filtered image
        'maxratio': 2D array - the maximum allowed ratio between pixel flux
                    and the median-filtered image
        'diff': the difference between the input image and an image representing
                the max allowed flux in each pixel.
        'differr': the expected error in the difference calculation.
        'weighted': boolean, set to True to use flat cal weights (see flatcal/ 
                    and util.obsFile.getPixelCountImage() )
        'niter': number of iterations performed.

    HISTORY:
        2/8/2013: Added iteration to the algorithm, to help in cases where some 
                  pixels are missed because of other nearby hot pixels. Added
                  input parameter 'maxIter' (also present in parameter file).
                  Added output dictionary key 'niter'.

    '''
    
    if image is not None:
        im = np.copy(image)      #So that we pass by value instead of by reference (since we will change 'im').
    else:
        im = None
    
    #Approximate peak/median ratio for an ideal (Gaussian) PSF sampled at 
    #pixel locations corresponding to the median kernal used with the real data.  
    gauss_array = utils.gaussian_psf(fwhm, boxSize)
    maxRatio = np.max(gauss_array) / np.median(gauss_array)

    if obsFile is None and im is None:
        obsFile = ObsFile.ObsFile(inputFileName)

    if im is None:
        print 'Counting photons per pixel'
        im = (obsFile.getPixelCountImage(firstSec=firstSec, integrationTime=intTime,
                                           weighted=weighted, getRawCount=True))['image']
        print 'Done'
    
    #Now im definitely exists, make a copy for display purposes later (before we change im).
    imOriginal = np.copy(im)
        
    #For now, assume 0 counts in a pixel means the pixel is dead.
    #Turn such pixel values into NaNs.
    im[im < 1] = np.nan
    
    oldHotMask = np.zeros(shape=np.shape(im), dtype=bool)  #Initialise a mask (all False) for comparison on each iteration.
    
    for iIter in range(maxIter):
        #Calculate median filtered image
        #Each pixel takes the median of itself and the surrounding boxSize x boxSize box.
        #(Not sure what edge effects there may be...)
        medFiltImage = utils.median_filterNaN(im, boxSize, mode='mirror') #Note - 'reflect' mode looks like it would repeat the edge row/column in the 'reflection'; 'mirror' does not, and makes more sense for this application.
    
        #Calculate difference between flux in each pixel and maxRatio * the median in the enclosing box.
        #Also calculate the error in the same quantity.
        diff = im - maxRatio * medFiltImage
        diffErr = np.sqrt(im + (maxRatio ** 2) * medFiltImage)
        
        #Any pixel that has a peak/median ratio more than nSigma above the maximum ratio should be flagged
        #(True=>bad pixel; False=> good pixel).
        hotMask = (diff > (nSigma * diffErr)) | oldHotMask

        #If no change between between this and the last iteration then stop iterating
        if np.all(hotMask == oldHotMask): break

        #Otherwise update 'oldHotMask' and set all detected bad pixels to NaN for the next iteration
        oldHotMask = np.copy(hotMask)
        im[hotMask] = np.nan
                
    
    if display:
        imForDisplay = np.copy(imOriginal)
        imForDisplay[np.isnan(imOriginal)] = 0  #Just because it looks prettier
        
        #utils.plotArray(im, cbar=True)
        fig = mpl.figure(figsize=(10, 10))
        mpl.matshow(imForDisplay, vmax=np.percentile(imForDisplay, 99.5), fignum=False)
        mpl.colorbar()
        x = np.arange(np.shape(imForDisplay)[1])
        y = np.arange(np.shape(imForDisplay)[0])
        xx, yy = np.meshgrid(x, y)
        
        if np.sum(hotMask) > 0: mpl.scatter(xx[hotMask], yy[hotMask], c='y')
        #if obsFile is None:
        #    plotTitle = ('im' + ' ' + str(firstSec) + '-' + str(firstSec + intTime) + 's')
        if obsFile is not None:
            plotTitle = (obsFile.fileName + ' ' + str(firstSec) + '-' + str(firstSec + intTime) + 's')
            mpl.suptitle(plotTitle)

    return {'mask':hotMask, 'image':im, 'medfiltimage':medFiltImage,
            'maxratio':maxRatio, 'diff':diff, 'differr':diffErr, 'niter':iIter + 1}








def findHotPixels(inputFileName=None, outputFileName=None,
                  paramFile=None,
                  timeStep=1, startTime=0, endTime= -1, fwhm=3.0,
                  boxSize=5, nSigma=3.0, display=False, weighted=False,
                  maxIter=5):
    '''
    To find hot pixels (and possibly cold pixels too at some point...).
    This routine is the main code entry point.
    
    
    INPUTS:
        paramFile - string, name of input parameter file. Other parameters 
                    provided as arguments will override values in the parameter
                    file. See example parameter file hotPixels.dict in pipeline
                    'params' directory.
        inputFileName - string, pathname of input observation file.
        outputFileName - string, pathname of output .h5 file.
        timeStep - integer (for now), check for hot pixels in intervals of 
                   timeStep seconds.
        startTime - integer (for now), number of seconds into exposure to begin
                    searching for hot pixels (default is 0, start of exposure).
        endTime - integer (for now), number of seconds into exposure to end at.
                  If endTime=-1, will continue to end of exposure.
        
        The following are as for checkInterval() and are passed on to that function:
        
        fwhm: Scalar float. Estimated full-width-half-max of the PSF (in 
                    pixels).
        boxSize: Scalar integer. Size of edge of box used for calculating median 
                    flux in the region surrounding each pixel. Larger
                    value => more robust against outliers, but probably worse
                    at accounting for real local flux variations. (5 is a good
                    value).
        nSigma: Scalar float. If the flux ratio for a pixel is (nSigma x expected
                    error) above the max expected given the PSF FWHM, then flag it 
                    as hot. Larger value => lower sensitivity.
        display: Boolean. If true, then display the input image and mark those 
                    flagged as hot with a dot.
        'weighted': boolean, set to True to use flat cal weights (see flatcal/ 
                    and and util.obsFile.getPixelCountImage() )
        
    OUTPUTS:
        Writes an hdf (.h5) file to outputFile containing tables of start/end
        times and flags for bad periods; one table for each pixel. This file 
        can most conveniently be read in using the readHotPix() function.
        
        However, the h5 structure is as follows:
        
        /header/header - 1-row table containing header info (currently just
                         nRow and nCol for the input obs File image)
        /timeMasks/ - group containing one table for each pixel.
    
        Each pixel table in /timeMasks/ is named for the x, y location of the 
        pixel, e.g., for the pixel at x=12, y=34, the table name is 'x12y34'
    
        The structure of each pixel table is as specified in
        headers.TimeMask.TimeMask - i.e., three columns, one row per 'bad' time
        interval. The columns are tBegin, tEnd, reason,
        representing the start and end times for each 'bad' interval within the
        exposure, and a pytables Enum representing indicating the reason for 
        which the pixel was tagged as bad for each interval.
 
        Note - times recorded in the output .h5 file are now in clock ticks.
                
        
    EXAMPLES: 
    
        findHotPixels('obs_20121211-024511.h5', 'testOutput.h5', 
                      'hotPixels.dict', startTime=0, endTime=5)
        
        - Writes time-mask data for 'obs_2012...' to file 'testOutput.h5', 
            using parameter file 'hotPixels.dict', and considering only 
            the first 5 seconds of the input file (here the start/end time 
            values override those in the parameter file).
        
        
    HISTORY:
        2/8/2013: Added iteration to the algorithm, now takes parameter 
                    'maxIter' (also added to parameter file). See
                    checkInterval() for more info.
    '''

    if paramFile is not None:
        print 'Reading parameter file...'
        params = readDict.readDict(paramFile)
        params.readFromFile(paramFile)
        if inputFileName is None: inputFileName = params['inputFileName']
        if outputFileName is None: outputFileName = params['outputFileName']
        if timeStep is None: timeStep = params['timeStep']
        if startTime is None: startTime = params['startTime']
        if endTime is None: endTime = params['endTime']
        if fwhm is None: fwhm = params['fwhm']
        if boxSize is None: boxSize = params['boxSize']
        if nSigma is None: nSigma = params['nSigma']
        if display is None: display = params['display']
        if maxIter is None: maxIter = params['maxIter']
    else:
        print 'No parameter file provided - using defaults/input params'
    
    #A few defaults that will be used in the absence of parameter file or 
    #arguments provided by the caller.
    if timeStep is None: timeStep = 1
    if startTime is None: startTime = 0
    if endTime is None: pass
    if maxIter is None: maxIter = 5
    
    obsFile = ObsFile.ObsFile(inputFileName)
    expTime = obsFile.getFromHeader('exptime')
    if endTime < 0: endTime = expTime
    stepStarts = np.arange(startTime, endTime, timeStep)  #Start time for each step (in seconds).
    stepEnds = stepStarts + timeStep                      #End time for each step
    nSteps = len(stepStarts)
    stepStartsTicks = stepStarts * obsFile.ticksPerSec
    stepEndsTicks = stepEnds * obsFile.ticksPerSec
        
    #Initialise stack of masks, one for each time step
    masks = np.zeros([obsFile.nRow, obsFile.nCol, nSteps], dtype=np.int8)

    #Get the mask for each time step
    for i, eachTime in enumerate(stepStarts):
        print str(eachTime) + ' - ' + str(eachTime + timeStep) + 's'
        masks[:, :, i] = checkInterval(obsFile=obsFile, firstSec=eachTime, intTime=timeStep,
                                     fwhm=fwhm, boxSize=boxSize, nSigma=nSigma,
                                     display=display, weighted=weighted,
                                     maxIter=maxIter)['mask']
 
    
    #Initialise an empty list that will eventually be a list of entries corresponding to:
    #  [xpos, ypos, badTimeList]
    #
    #where badTimeList is a nested list of tuples indicating bad time intervals and 
    #reasons for flagging:
    #   
    #  [(start time 1, end time1, flag number 1), (Start time 2, end time 2, flag 2), ...] 
    #
    #(with times in seconds).
    
    timeMaskData = []  


    #Convert the masks to bad-time lists.
    for iRow in range(obsFile.nRow):
        for iCol in range(obsFile.nCol):
            flagSequence = masks[iRow, iCol, :]

            #What (unique) non-zero flags are listed for this pixel?
            uniqueFlags = [x for x in set(flagSequence) if x != 0]   #What non-zero flags are listed for this pixel?

            #Initialise a list for bad times for current pixel
            badTimeList = []
            
            for eachFlag in uniqueFlags:
                #Make a list of intervals for each bad timestep - e.g. one interval for every second if timestep is seconds and the pixel is always bad.
                badStepIntervals = [interval([stepStartsTicks[i], stepEndsTicks[i]]) 
                                    for i in range(nSteps) if flagSequence[i] == eachFlag]  #In units of ticks (not seconds).
                #Find the union of those intervals (concatenate adjacent intervals)
                badIntervals = interval().union(badStepIntervals)
                #And then separate out the individual components after concatenating...
                for eachComponent in badIntervals.components:  #Slightly annoying way you have to do it with the Interval class
                    badTimeList.append((eachComponent[0][0], eachComponent[0][1], eachFlag)) #Begin time, end time, flag number   

            #Now we've gone through all the flags for this pixel, sort all the
            #entries in badTimeList by start time (there shouldn't be any overlap between
            #intervals with different flags, since this is all one pixel!)
            badTimeList.sort(key=lambda x: x[0])

            #Create a new entry in the 'bad times' table with all the bad times for this pixel. 
            timeMaskData.append([iCol, iRow, badTimeList])

    #End looping through pixels.
    
    
    #Write it all out to the .h5 file
    writeHotPixels(timeMaskData, obsFile, outputFileName)
    




def writeHotPixels(timeMaskData, obsFile, outputFileName):
    '''
    Write the output hot-pixel time masks table to an .h5 file. Called by
    findHotPixels().
    
    INPUTS:
        timeMaskData - list structure as constructed by findHotPixels()
        obsFile - the ObsFile object from which the data to be written
                        was derived.
        outputFileName - string, the pathname of the .h5 file to write to.
    
    OUTPUTS:
        Writes an .h5 file to filename outputFileName. See findHotPixels()
        for full description of the output data structure.
        
    HISTORY:
        2/15/2013: Updated so that behaviour of outputFileName is consistent with the
        behaviour of the input file name for an ObsFile instance. (i.e., 
        unless the path provided is absolute, $MKID_DATA_DIR is prepended
        to the file name.
    '''
    
    if (os.path.isabs(outputFileName)):
        #self.fileName = os.path.basename(fileName)
        fullFileName = outputFileName
    else:
        #self.fileName = fileName
        # make the full file name by joining the input name 
        # to the MKID_DATA_DIR (or . if the environment variable 
        # is not defined)
        dataDir = os.getenv('MKID_DATA_DIR', '/')
        fullFileName = os.path.join(dataDir, outputFileName)
    
    fileh = tables.openFile(fullFileName, mode='w')
    
    try:    
        #Construct groups for header and time mask data.
        headerGroup = fileh.createGroup("/", headerGroupName, 'Header')
        timeMaskGroup = fileh.createGroup("/", dataGroupName, 'Time masks for temporarily bad pixels')


        #Fill in header info (just a one-row table)
        headerTable = fileh.createTable(headerGroup, headerTableName, headerDescription,
                                        'Header Info')
        header = headerTable.row
        header[obsFileColName] = obsFile.fileName
        header[nColColName] = max([x[0] for x in timeMaskData]) + 1  #Assume max. x value represents number of columns
        header[nRowColName] = max([x[1] for x in timeMaskData]) + 1  #Same for rows.
        header[ticksPerSecColName] = obsFile.ticksPerSec
        header[expTimeColName] = obsFile.getFromHeader('exptime')    #Newly implemented - SHOULD DOUBLE CHECK! Should automatically account for any of Matt's time corrections JvE 7/08/2013.
        header.append()
        headerTable.flush()
        headerTable.close()

        #Fill in the time mask info          
        for eachPixelEntry in timeMaskData:
            
            #One table for every pixel:
            tableName = constructDataTableName(eachPixelEntry[0], eachPixelEntry[1]) #Table name from x,y pos.    
            timeMaskTable = fileh.createTable(timeMaskGroup, tableName, tm.TimeMask,
                                      'Time mask for pix. ' + tableName)
            row = timeMaskTable.row
            for eachPeriod in eachPixelEntry[2]:
                row['tBegin'] = eachPeriod[0]
                row['tEnd'] = eachPeriod[1]
                if eachPeriod[2] == 1: 
                    row['reason'] = tm.timeMaskReason['hot pixel']
                else:
                    row['reason'] = tm.timeMaskReason['unknown'] #For now - other reasons (e.g. cold pixels) can be added later.
                row.append()
                
            timeMaskTable.flush()
            timeMaskTable.close()
    
    finally:
        fileh.close()




    
def readHotPixels(inputFile,nodePath=None):
    '''
    To read in a hot-pixels HDF file as written by findHotPixels(). 
    (Note 'hot pixels' may later include cold pixels and possibly other
     such things as well).
    
    INPUTS:
        inputFile - either: a pytables hot-pixel file instance; a pytables 
                    photon-list file instance from which to extract the same
                    information; or the pathname of either of these. If a file
                    instance, the file will be left open on completion, while if 
                    a pathname is given, the file will be opened and then closed
                    on completion.
        nodePath - optionally supply the location of the PyTables group
                   within the HDF hierarchy that contains the hot pixel 
                   information. Default is just the root ('/') of the .h5
                   file, which is appropriate for standard hot pixel files.
                   However, by setting this appropriately, you can instead 
                   supply a photon list file in 'inputFile', since the hot pixel
                   data hierarchy is copied directly into a sub-group within the
                   hot pixel file.
    
    OUTPUTS:
        Returns a dictionary with the following info:
        
            'nRow' - number of rows in original obs File.
            'nCol' - number of columns.
            'intervals' - nRow x nCol array of lists of pyinterval 'interval' objects.
                          The intervals represent time intervals where the pixel went bad
                          in *seconds* (although the values are stored in clock ticks
                          in the .h5 file read in). Where there are no bad intervals
                          for a pixel, its list is empty. Note the list for a given 
                          interval is not unioned into a single 'interval' object
                          since there may be different 'reasons' for the different
                          intervals!
            'reasons' - nRow x nCol array of lists of 'timeMaskReason' enums (see 
                        headers/TimeMask.py). Entries in these lists correspond directly 
                        to entries in the 'intervals' array.            
            'reasonEnum' - the enum instance for the 'reasons', so the 'concrete values'
                           stored in that array can be parsed back to their enum labels.
            'obsFileName' - the name of the original obs. file from which the
                            hot pixel masking info was derived.
            'ticksPerSecond' - number of clock ticks per second assumed in creating
                               the output file.
                               
        Note - would probably make more sense to return an object here at some point, instead
        of a dictionary....
    
    
    EXAMPLES:
        
        Make a hot pixel file and read it in to 'hpData':
        
        import hotPixels as hp
        findHotPixels('obs_20121211-024511.h5', 'testOutput.h5', 
                      'hotPixels.dict', startTime=0, endTime=5)
        hpData = hp.readHotPixels('testoutput.h5')
        

        Find out how many discrete time intervals were flagged for pixel 44,5 
        (which happens to be bad at two times between 0-5s for obs file 
        'obs_20121211-024511.h5')
            
            >>> len(hpData['reasons'][44,5])
            2
        

        Find out the reason for the *first* time that this pixel was flagged:
        
            >>> enum = hpData['reasonEnum']
            >>> enum(hpData['reasons'][44,5][0])
            'hot Pixel'
        

        Find the time range for which the same pixel was flagged the *second*
        time:
            
            >>> hpData['intervals'][44,5][1]
            interval([2.0], [3.0])
        
            (i.e. from time=2sec to 3sec).
        

        Make an array containing the number of bad time intervals for each pixel:
            
            >>> nIntervals = np.vectorize(len)(hpData['intervals'])
        

        Make a boolean mask with True for all pixels where there was ANY bad time:
        
            >>> mask = nIntervals > 0
            

        Find the coordinates of all pixels which were bad at any time:
        
            >>> np.where(nIntervals > 0)
            (array([ 0,  0,  3, ...]),
             array([26, 34,  5, ...]))


        Check if pixel was bad at time 2.5sec, 1.5sec, and 0.8sec:
                    
            >>> 2.5 in interval.union(hpData['intervals'][44,5])
            True
        
            >>> 1.5 in interval.union(hpData['intervals'][44,5])
            False
        
            >>> 0.8 in interval.union(hpData['intervals'][44,5])
            True
            
            
    '''
    
    if type(inputFile) is str:
        fileh = tables.openFile(inputFile, mode='r')
    elif type(inputFile) is tables.file.File:
        fileh = inputFile
    else: raise ValueError('inputFile must be either a pathname string or a PyTables file instance.')
    
    if nodePath is None:
        if '/timemask/'+dataGroupName in fileh:
            nodePath = '/timemask/'
        elif '/timeMask/'+dataGroupName in fileh:
            nodePath = '/timeMask/'
        elif '/'+dataGroupName in fileh:
            nodePath = '/'
        
    if nodePath+'/'+dataGroupName not in fileh:
        raise RuntimeError, 'Time mask data node not found in HDF file: '+fileh.filename
    
    try:        
        #Get the header info
        headerInfo = fileh.getNode(nodePath + headerGroupName, headerTableName)[0] #The one row from the header
        nRow = headerInfo[nRowColName]
        nCol = headerInfo[nColColName]
        obsFileName = headerInfo[obsFileColName]
        ticksPerSec = headerInfo[ticksPerSecColName]
    
        #Intialise two ragged object arrays, one to take lists of Interval objects
        timeIntervals = np.empty((nRow, nCol), dtype='object')
        timeIntervals.fill([])
        #And one to take lists of corresponding flags
        reasons = np.empty((nRow, nCol), dtype='object')
        reasons.fill([])
        
        #Read in the data and fill in the arrays
        for iRow in range(nRow):
            for iCol in range(nCol):
                tableName = constructDataTableName(iCol, iRow)   #'x'+str(iCol)+'y'+str(iRow)
                #Get the list for the current pixel
                eventListTable = fileh.getNode(nodePath + dataGroupName, name=tableName)
                reasonEnum = eventListTable.getEnum('reason')  #Gets read in once for every table, but they should all be the same...
                timeIntervals[iRow, iCol] = \
                    [interval([eachRow['tBegin'], eachRow['tEnd']]) / ticksPerSec for eachRow
                      in eventListTable]        #Get the times in seconds (not ticks). No doubt this can be sped up if necessary...
                reasons[iRow, iCol] = [eachRow['reason'] for eachRow in eventListTable]
                    
        #Return a simple dictionary
        return {"intervals":timeIntervals, "reasons":reasons,
                "reasonEnum":reasonEnum, "nRow":nRow, "nCol":nCol,
                "obsFileName":obsFileName, "ticksPerSecond":ticksPerSec}


    finally:
        #If a filename was passed as a string (as opposed to a file instance) then close the file.
        if type(inputFile) is str: fileh.close()
    


def getEffIntTimeImage(hotPixDict,integrationTime,firstSec=0):
    '''    
    Get the total effective exposure time for each pixel after subtracting 
    any intervals where a pixel was masked as hot or bad.
    
    INPUTS:
        hotPixDict -  a hot pixels dictionary as returned by hotPixels.readHotPixels()
        firstSec - Start time (sec) to start calculations from, starting from
                    the beginning of the exposure to which timeMask refers.
        integrationTime - Length of integration time (sec) from firstSec to include
                    in the calculation. NOTE - Don't give an integration time
                    that goes beyond the end of the exposure! Currently does not (always)
                    have direct access to the total exposure time, so you can't set 
                    integrationTime=-1 and hope to integrate to the end of the 
                    exposure for this routine. As it stands, this function just
                    subtracts off the hot-pixel intervals from whatever integration
                    start/length you provide it without regard to the exposure length.

    RETURNS:
        A 2D array representing the total effective exposure time
        in seconds for each pixel in the detector array.
    '''
    
    #Figure out what time represents the end of the integration
    #if integrationTime == -1 or integerIntTime > len(pixelData):
    #    lastSec = int(np.floor(firstSec))+len(pixelData)
    #else:
    lastSec = firstSec + integrationTime
    outsideIntegration = interval([-np.inf, firstSec], [lastSec, np.inf])
    integrationInterval = interval([firstSec, lastSec])
    effectiveIntTimes = np.zeros((hotPixDict['nRow'],hotPixDict['nCol']),dtype=float)
    effectiveIntTimes.fill(np.nan)
    
    for iRow in np.arange(hotPixDict['nRow']):
        for iCol in np.arange(hotPixDict['nCol']):
            #Get the unioned (possibly multi-component) bad interval for this pixel.
            #(As in ObsFile.getPixelBadTimes)
            allBadIntervals = interval.union(hotPixDict['intervals'][iRow, iCol])
            #Get intersection of integration time interval and the bad time intervals.
            maskedIntervals = allBadIntervals & integrationInterval
            effectiveIntTimes[iRow,iCol] = integrationTime - utils.intervalSize(maskedIntervals)
    
    return effectiveIntTimes


def getHotPixels(hotPixDict,integrationTime=-1,firstSec=0):
    '''
    Return a boolean array indicating which pixels went bad at any point
    during the specified integration time.
    
    INPUTS:
        hotPixDict -  a hot pixels dictionary as returned by hotPixels.readHotPixels()
        firstSec - Start time (sec) to start calculations from, starting from
                    the beginning of the exposure to which timeMask refers.
        integrationTime - Length of integration time (sec) from firstSec to include
                    in the calculation. NOTE - currently does not (always)
                    have direct access to the total exposure time, so if you set 
                    integrationTime=-1, it'll hopefully give results good to the
                    end of the exposure, assuming that the measured hot pixel time mask
                    doesn't somehow extend over the end of the exposure. But that's not
                    totally 100% guaranteed at this point (at least as far as I can think right now).
                    
    RETURNS:
        A 2D Boolean array matching the size/shape of the detector image. True indicates
        a pixel that went bad between firstSec and firstSec+integrationTime, and False 
        indicates that the pixel was okay during that time.
    '''
    
    if integrationTime == -1:
        intTimeInternal = np.Inf
    else:
        intTimeInternal = integrationTime
    lastSec = firstSec + intTimeInternal
    
    #outsideIntegration = interval([-np.inf, firstSec], [lastSec, np.inf])
    integrationInterval = interval([firstSec, lastSec])
    badPix = np.zeros((hotPixDict['nRow'],hotPixDict['nCol']),dtype=bool)
    badPix.fill(False)
    
    for iRow in np.arange(hotPixDict['nRow']):
        for iCol in np.arange(hotPixDict['nCol']):
            #Get the unioned (possibly multi-component) bad interval for this pixel.
            #(As in ObsFile.getPixelBadTimes)
            allBadIntervals = interval.union(hotPixDict['intervals'][iRow, iCol])
            #Get intersection of integration time interval and the bad time intervals.
            maskedIntervals = allBadIntervals & integrationInterval
            #Figure out if there were any bad intervals during the integration time or not.
            if len(maskedIntervals) != 0:
                badPix[iRow,iCol]=True
    
    return badPix



if __name__ == "__main__":
    
    paramFile = None
    inputFile = None
    outputFile = None
    
    nArg = len(sys.argv)
    if nArg == 1: paramFile = 'hotPixels.dict'
    if nArg > 1: paramFile = sys.argv[1]
    if nArg > 2: inputFile = sys.argv[2]
    if nArg > 3: outputFile = sys.argv[3]
    findHotPixels(paramFile=paramFile, inputFileName=inputFile,
                  outputFileName=outputFile)
    
    









#===============================================================================
# def rebinND(a, *args):
#    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
#    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
#    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
#    example usages:
#    >>> a=rand(6,4); b=rebin(a,3,2)
#    >>> a=rand(6); b=rebin(a,2)
# 
#    Taken directly from SciPy cookbook - http://www.scipy.org/Cookbook/Rebinning
#    '''
#    shape = a.shape
#    lenShape = len(shape)
#    factor = np.asarray(shape) / np.asarray(args)
#    evList = ['a.reshape('] + \
#       ['args[%d],factor[%d],' % (i, i) for i in range(lenShape)] + \
#       [')'] + ['.sum(%d)' % (i + 1) for i in range(lenShape)] + \
#       ['/factor[%d]' % i for i in range(lenShape)]
#    print ''.join(evList)
#    return eval(''.join(evList))
#===============================================================================


