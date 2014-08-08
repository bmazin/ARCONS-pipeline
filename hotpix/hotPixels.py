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
                intervals.
                
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


History/notes:
    - COLD PIXEL MASKING SWITCHED OFF FOR NOW - 5/6/2014. 


To do:
    - Extend to check for 'cold' pixels. DONE/IN-PROGRESS 11/22/2013.
    - Suppress annoying invalid value warnings that come up (and can be ignored)
    - Speed up (most time is spent reading the h5 file; significant overhead can
        maybe be saved by reading in the whole file at once instead of
        second-by-second).
    - At some point, perhaps update to allow sub-second time steps. THINK THIS SHOULD WORK OKAY NOW - BUT NEED TO TEST TO BE SURE. 11/27/2013.
    - Time intervals are currently stored in units of seconds - need to convert
        to clock ticks for consistency with TimeInterval class. FIXED - 02/12/2013.
    - If necessary, apply algorithm only at the red end, where hot-pixel
        contrast should be higher. MAYBE WORTH REVISITING, BUT LOOKS LIKE NOT *ALL* HOT PIXELS ARE RED....
    - Output files seem unnecessarily large. May be some way of making this more
        efficient.
    
    - NOTE - HAVE NOT PROPERLY CODED TO MAKE SURE THAT TIME STEPS ARE PROPERLY
      CONSECUTIVE IN INTEGER 'TICKS' FOR THE RECORDED INTERVALS - SHOULD PROBABLY
      LOOK INTO THIS.

See individual routines for more detail.

'''

import os.path
import sys
import pickle
from math import *
from interval import interval
import tables
import ds9
import numpy as np
import numpy.ma as ma
import scipy.ndimage.filters as spfilters
import matplotlib.pylab as mpl
from matplotlib.colors import LogNorm
import astropy.stats
import util.ObsFile as ObsFile
import util.utils as utils
import headers.TimeMask as tm
from headers import pipelineFlags as pflags
import util.readDict as readDict

headerGroupName = 'header'  #Define labels for the output .h5 file.
headerTableName = 'header'
dataGroupName = 'timeMasks' #data *table* name is constructed on a per-pixel basis (see constructDataTableName)
nRowColName = 'nRow'        #Name of the no.-of-rows column in the header for the output .h5 file.
nColColName = 'nCol'        #Ditto for no.-of-columns column.
startTimeColName = 'startTime'  #Header column name for recording start time of time-mask data, in sec from beginning of obs. file.
endTimeColName = 'endTime'      #Ditto for end time.
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
    expTime = tables.Float64Col(dflt=np.nan)        #So that per-pixel effective exposure times can be calculated any time without having to refer to the original obs file. Added 6/21/2013, JvE, although not yet used elsewhere in the code.
    startTime = tables.Float64Col(dflt=np.nan)      #To record start and end times within the obs file for which the time masks were created.
    endTime = tables.Float64Col(dflt=np.nan)
 
def checkInterval(firstSec=None, intTime=None, fwhm=4.0, boxSize=5, nSigmaHot=3.0,
                  nSigmaCold=3.0, obsFile=None, inputFileName=None, image=None,
                  display=False, ds9display=False, dispToPickle=None, weighted=False,
                  fluxWeighted=False, maxIter=5, dispMinPerc=0.0, dispMaxPerc=98.0, 
                  diagnosticPlots=False, useLocalStdDev=False, useRawCounts=True,
                  bkgdPercentile=10.0):
    '''
    To find the hot, cold, or dead pixels in a given time interval for an observation file.
    This is the guts of the bad pixel finding algorithm, but only works on a single time
    interval. Compares the ratio of flux in each pixel to the median of the flux in an
    enclosing box. If the ratio is too high -- i.e. the flux is too tightly 
    distributed compared to a Gaussian PSF of the expected FWHM -- then the 
    pixel is flagged as bad.
    
    Cold pixels are indentified as those with a flux significantly below the median
    flux of the surrounding pixels.
    
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
        nSigmaHot: Scalar float. If the flux ratio for a pixel is (nSigmaHot x expected
                    error) above the max expected given the PSF FWHM, then flag it 
                    as hot.
        nSigmaCold: If the flux for a pixel is more than nSigmaCold*expected error
                    below the median in the surrounding box, then flag it as cold.
        display: Boolean. If true, then display the input image and mark those 
                    flagged as hot with a dot.
        ds9display: Boolean, as for 'display', but displays the output in ds9 instead.
        maxIter: Scalar integer. Maximum number of iterations allowed.
        dispMinPerc: Lower percentile for image stretch if display=True
        dispMaxPerc: Upper percentile "     "     "      "       "
        dispToPickle: If not None (default), then save a pickle file of the data fo
                    the main plot shown with display=True and ds9display=True. If a string
                    value, then this value is used as the name for the output pickle file.
                    Otherwise uses default output name 'badPixels.pickle'.
                    Saves a dictionary with four entries:
                        "image" - the integrated image for the interval requested
                        "hotMask" - mask of hot pixels detected (False=good, True=hot)
                        "coldMask" - similar mask for the cold pixels (False=good)
                        "deadMask" - similar mask for the dead pixels (False=good)
        diagnosticPlots: if True, shows a bunch of diagnostic plots for debug purposes.
        useLocalStdDev - if True, use the local (robust) standard deviation within the 
                    moving box for the sigma value in the hot pixel thresholding
                    instead of Poisson statistics. Mainly intended for situations where
                    you know there is no astrophysical source in the image (e.g. flatfields,
                    laser calibrations), where you can also set fwhm=np.inf
        weighted: boolean, set to True to use flat cal weights (see flatcal/ 
                    and util.ObsFile.getPixelCountImage() ).
        fluxWeighted: boolean, if True, flux cal weights are applied (also see fluxcal
                    and util.ObsFile.getPixelCountImage() ) Added JvE 7/16/2014
        useRawCounts - Boolean. if True, creates the mask on the basis of the raw photon counts summed
                        across all wavelengths in the obs file. **If False, input must be provided
                        in obsFile** - in which case, whatever calibrations are already loaded and
                        (and switched on) in the obsFile instance will be applied before looking
                        for badly behaved pixels. *Overrides weighted and fluxWeighted*. Added JvE 7/16/2014
        bkgdPercentile - percentile level (in %) in image to use as an estimate of the background.
                         ***SHOULD BE ADDED AS A PARAMETER TO THE PARAMETER FILE...!!****

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
        'niter': number of iterations performed.

    HISTORY:
        2/8/2013: Added iteration to the algorithm, to help in cases where some 
                  pixels are missed because of other nearby hot pixels. Added
                  input parameter 'maxIter' (also present in parameter file).
                  Added output dictionary key 'niter'.
        5/6/2014: Switched off cold pixel flagging for now - cutoff between
                  genuinely bad cold pixels and those with low QE is not
                  unambiguous enough.
        6/7/2014: Added option to use local box robust estimate of sigma instead
                  of sigma from photon noise estimate.

    '''
    
    doColdFlagging = False      #Switch off cold flagging for now, 5/6/2014. Hard coded in as this is
                                #not a user option at this point....
    #assert 1==0
    defaultPklFileName = 'badPixels.pickle'

    if useRawCounts is False and obsFile is None:
        raise ValueError, 'Must provide obsFile object if you want to use calibrated (not raw) counts'

    if image is not None:
        im = np.copy(image)      #So that we pass by value instead of by reference (since we will change 'im').
    else:
        im = None
    
    #Approximate peak/median ratio for an ideal (Gaussian) PSF sampled at 
    #pixel locations corresponding to the median kernal used with the real data.  
    gauss_array = utils.gaussian_psf(fwhm, boxSize)     #***MAYBE BETTER IF WE OVERSAMPLE THIS?****
    maxRatio = np.max(gauss_array) / np.median(gauss_array)

    if obsFile is None and im is None:
        obsFile = ObsFile.ObsFile(inputFileName)

    #Useful generic subtitle for various plots
    if obsFile is not None:
        plotTitle = (obsFile.fileName + ' ' + str(firstSec) + '-' + str(firstSec + intTime) + 's')
    else:
        plotTitle = ''

    if im is None:
        print 'Getting image time-slice'
        im = (obsFile.getPixelCountImage(firstSec=firstSec, integrationTime=intTime,
                                           weighted=weighted, fluxWeighted=fluxWeighted, 
                                           getRawCount=useRawCounts))['image']
        print 'Done'
    
    #Now im definitely exists, make a copy for display purposes later (before we change im).
    imOriginal = np.copy(im)
        
    #For now, assume 0 counts in a pixel means the pixel is dead.
    #Turn such pixel values into NaNs.
    deadMask = (im<1)
    im[im < 1] = np.nan
    
    oldHotMask = np.zeros(shape=np.shape(im), dtype=bool)   #Initialise a mask for hot pixels (all False) for comparison on each iteration.
    oldColdMask = np.zeros(shape=np.shape(im), dtype=bool)  #Same for cold pixels
    hotMask = np.zeros(shape=np.shape(im), dtype=bool)      
    coldMask = np.zeros(shape=np.shape(im), dtype=bool)     
    
    #Initialise some arrays with nan's in case we don't get to fill them out for real,
    #just so that things don't go awry when we try to return the arrays at the end.
    medFiltImage = np.zeros_like(im)
    medFiltImage.fill(np.nan)
    diff = np.zeros_like(im)
    diff.fill(np.nan)
    diffErr = np.zeros_like(im)
    diffErr.fill(np.nan)
    #Ditto for number of iterations
    iIter=-1 
    
    if not np.all(deadMask):  #Check to make sure not *all* the pixels are dead before doing further masking.
        for iIter in range(maxIter):
            print 'Iteration: ',iIter
            #Calculate median filtered image
            #Each pixel takes the median of itself and the surrounding boxSize x boxSize box.
            #(Not sure what edge effects there may be, ignore for now...)
            #Note - 'reflect' mode looks like it would repeat the edge row/column in the 'reflection';
            #'mirror' does not, and makes more sense for this application.
            #Do the median filter on a NaN-fixed version of im.
            nanFixedImage = utils.replaceNaN(im, mode='nearestNmedian', boxsize=boxSize**2-1)
            assert np.all(np.isfinite(nanFixedImage))  #Just make sure there's nothing weird still in there.
            medFiltImage = spfilters.median_filter(nanFixedImage, boxSize, mode='mirror')
            #medFiltImage = utils.median_filterNaN(im, boxSize, mode='mirror')  #Original version without interpolating the NaNs
            
            overallMedian = np.median(im[~np.isnan(im)])
            overallBkgd = np.percentile(im[~np.isnan(im)],bkgdPercentile)
        
            if doColdFlagging is True or useLocalStdDev is True:
                stdFiltImage = utils.nearestNrobustSigmaFilter(im, n=boxSize**2-1)    
            
            #-------------- Cold flagging switched off for now, May 6 2014-----------------    
            if doColdFlagging is True:
                nrstNbrMedFiltImage = utils.nearestNmedFilter(im, n=boxSize**2-1)  #Possibly useful with cold pixel flagging.
                overallStdDev = astropy.stats.median_absolute_deviation(im[~np.isnan(im)])*1.4826
                #Calculate the standard-deviation filtered image,
                #using a kernel footprint that will miss out the central pixel:
                #footprint = np.ones((boxSize,boxSize))
                #footprint[boxSize/2,boxSize/2] = 0      #Can move this definition outside the loop if needed, but prob. okay at the moment. 
                #stdFiltImage = utils.stdDev_filterNaN(im, boxSize, mode='mirror', footprint=footprint)
                #stdFiltImage = utils.nearestNstdDevFilter(im, n=boxSize**2-1)       #TRYING OUT USING A NEAREST-NEIGHBOUR STD. DEV FILTER....
            #----------------------------------------------------------------------------------
    
            #Calculate difference between flux in each pixel and maxRatio * the median in the enclosing box.
            #Also calculate the error that would exist in a measurment of a pixel that *was* at the peak of a real PSF
            print 'overallMedian: ',overallMedian
            print 'overallBkgd: ',overallBkgd
            print 'maxRatio: ',maxRatio
            diff = im - maxRatio * medFiltImage + (maxRatio-1.)*overallBkgd   #****TESTING ACCOUNTING FOR BACKGROUND - SEE LAST TERM****
    
            #Simple estimate, probably makes the most sense: photon error in the max value allowed. Neglect errors in the median itself here.
            if useLocalStdDev is False:
                diffErr = np.sqrt(maxRatio * medFiltImage)       
                #Alternatively, the corrected version of what I was trying to do before - i.e., the error in diff, which
                #seems bogus because if you have a very high value in im, then you'll have a large error, which is
                #not what you're looking for.
                #diffErr = np.sqrt(im + (maxRatio ** 2) * stdFiltImage**2 / (boxSize**2))        #/boxSize**2 because it's like error in the mean, i.e., divide by sqrt(n). Def. not rigorous to apply that to a median, but, better than nothing....
                #
                #Originally was something like:
                #diffErr = np.sqrt(im + (maxRatio ** 2) * medFiltImage)      #which is really not very sensible.
            else:
                diffErr = stdFiltImage
            
            if iIter == 0:
                diffOriginal = np.copy(diff)
                diffErrOriginal = np.copy(diffErr)
                
            #Any pixel that has a peak/median ratio more than nSigma above the maximum ratio should be flagged as hot:
            #(True=>bad pixel; False=> good pixel).
            hotMask = (diff > (nSigmaHot * diffErr)) | oldHotMask
            
            #-------------- Cold flagging switched off for now, May 6 2014-----------------        
            if doColdFlagging is True:
                #And any pixel that is more than nSigma *below* the std. dev. of the surrounding box (not including itself)
                #should be flagged as cold:
                coldMask = ((nrstNbrMedFiltImage - im) > nSigmaCold * stdFiltImage) | oldColdMask 
                #coldMask = ((overallMedian - im) > nSigmaCold * overallStdDev) | oldColdMask 
            #------------------------------------------------------------------------------------------
            
            if diagnosticPlots is True and iIter==0:
                #Display a histogram of fluxes by pixel
                mpl.figure()
                mpl.hist(im[~np.isnan(im)],bins=400)
                                                #(np.nanmax(im)-np.nanmin(im)/np.median(im[~np.isnan(im)]))*10.)
                                                 #np.sqrt(np.sum(~np.isnan(im))))
                mpl.xlabel('Photon counts')
                mpl.ylabel('Number of pixels')
                mpl.title('Flux by Pixel, Before Flagging')
                mpl.suptitle(plotTitle)
                    
            #If no change between between this and the last iteration then stop iterating
            if np.all(hotMask == oldHotMask) and np.all(coldMask == oldColdMask): break
    
            #Otherwise update 'oldHotMask' and set all detected bad pixels to NaN for the next iteration
            oldHotMask = np.copy(hotMask)
            im[hotMask] = np.nan
            oldColdMask = np.copy(coldMask)
            im[coldMask] = np.nan
    
    #Finished with loop, now combine the hot and cold masks:
    assert np.all(coldMask & hotMask == False)  #Presumably a pixel can't be both hot AND cold....
    assert np.all(hotMask & deadMask == False)  #Ditto hot and dead. (But *cold* and dead maybe okay at this point).
    mask = np.empty_like(hotMask,dtype=int)
    mask.fill(pflags.badPixCal['good'])
    mask[hotMask] = pflags.badPixCal['hot']
    mask[coldMask] = pflags.badPixCal['cold']
    mask[deadMask] = pflags.badPixCal['dead']    
    
    
    if display or ds9display or (dispToPickle is not False):
        imForDisplay = np.copy(imOriginal)
        cleanImForDisplay = np.copy(imForDisplay)
        imForDisplay[np.isnan(imOriginal)] = 0  #Just because it looks prettier
        cleanImForDisplay[mask!=pflags.badPixCal['good']] = 0   #An image with only good pixels
    
        vmin=np.percentile(imForDisplay,dispMinPerc)
        vmax=np.percentile(imForDisplay,dispMaxPerc)

        x = np.arange(np.shape(imForDisplay)[1])
        y = np.arange(np.shape(imForDisplay)[0])
        xx, yy = np.meshgrid(x, y)
        
        if display:
            fig = mpl.figure(figsize=(5,5))
            mpl.matshow(imForDisplay,vmin=vmin,vmax=vmax,
                        fignum=False,origin='lower',cmap=mpl.cm.hot)     #, norm=LogNorm())  #, cmap=mpl.cm.hot)
            mpl.colorbar()
            if np.sum(hotMask) > 0: mpl.scatter(xx[hotMask], yy[hotMask], marker='x', c='b', label='Hot', linewidths=2)
            if np.sum(coldMask) > 0: mpl.scatter(xx[coldMask], yy[coldMask], marker='o', c='w', label='Cold')
            if np.sum(deadMask) > 0: mpl.scatter(xx[deadMask], yy[deadMask], marker='o', c='b', label='Dead')
            mpl.legend()
            #if obsFile is None:
            #    plotTitle = ('im' + ' ' + str(firstSec) + '-' + str(firstSec + intTime) + 's')
            mpl.title('Image + mask')
            mpl.suptitle(plotTitle)
            utils.showzcoord()

        if ds9display:
            utils.ds9Array(imForDisplay, normMin=vmin, normMax=vmax, colormap='bb')
            d=ds9.ds9()     #Get reference to the now (hopefully) open ds9 instance.
            for i in range(np.sum(hotMask)):
                d.set("regions command {point "+str(xx[hotMask][i]+1)+" "+str(yy[hotMask][i]+1)
                      +"  #point=x color=blue}")
            for i in range(np.sum(coldMask)):
                d.set("regions command {circle "+str(xx[coldMask][i]+1)+" "+str(yy[coldMask][i]+1)
                      +" 0.3 #color=white}")
            for i in range(np.sum(deadMask)):
                d.set("regions command {circle "+str(xx[deadMask][i]+1)+" "+str(yy[deadMask][i]+1)
                      +" 0.3 #color=blue}")

        if dispToPickle is not False:
            #Save to pickle file to feed into a separate plotting script, primarily for 
            #the pipeline paper.
            if type(dispToPickle) is str:
                pklFileName = dispToPickle
            else:
                pklFileName = defaultPklFileName
            pDict = {"image":imForDisplay,"hotMask":hotMask,"coldMask":coldMask,
                     "deadMask":deadMask}
            print 'Saving to file: ',pklFileName
            output = open(pklFileName, 'wb')
            pickle.dump(pDict,output)
            output.close()

        #Show diagnostic images (not in ds9).
        if diagnosticPlots is True:
            #fig = mpl.figure(figsize=(5,5))
            #mpl.matshow(cleanImForDisplay, vmax=np.percentile(imForDisplay, 98.0),
            #            cmap=mpl.cm.hot, fignum=False, origin='lower')
            #mpl.colorbar()
            #mpl.title('Cleaned image')

            print 'Max ratio: ', maxRatio

            if doColdFlagging is True:
                print 'Overall median: ',overallMedian
                print 'Overall std. dev.: ',overallStdDev

            fig = mpl.figure(figsize=(5,5))
            imToPlot = np.copy(medFiltImage)
            imToPlot[np.isnan(imToPlot)] = 0
            mpl.matshow(imToPlot, vmax=np.percentile(imToPlot[np.isfinite(imToPlot)], 100.0),
                        fignum=False, origin='lower', cmap=mpl.cm.hot)
            mpl.colorbar()
            mpl.title('Median filtered image')
            mpl.suptitle(plotTitle)
            utils.showzcoord()

            
            fig = mpl.figure(figsize=(5,5))
            imToPlot = np.copy(diffOriginal)
            imToPlot[np.isnan(imToPlot)] = np.nanmin(imToPlot)
            mpl.matshow(imToPlot, vmax=np.percentile(imToPlot[np.isfinite(imToPlot)], 98.5),
                        fignum=False, origin='lower', cmap=mpl.cm.hot)
            mpl.colorbar()
            mpl.title('Difference image')
            mpl.suptitle(plotTitle)
            utils.showzcoord()

            
            fig = mpl.figure(figsize=(5,5))            
            imToPlot = np.copy(diffErrOriginal)
            imToPlot[np.isnan(imToPlot)] = 0
            mpl.matshow(imToPlot, vmax=np.percentile(imToPlot[np.isfinite(imToPlot)], 100.0),
                        fignum=False, origin='lower', cmap=mpl.cm.hot)
            mpl.colorbar()
            mpl.title('Difference Error')
            mpl.suptitle(plotTitle)
            utils.showzcoord()

            
            fig = mpl.figure(figsize=(5,5))            
            imToPlot = utils.replaceNaN(np.copy(im),mode='nearestNmedian',boxsize=24)
            imToPlot[np.isnan(imToPlot)] = 0
            mpl.matshow(imToPlot, vmax=np.percentile(imToPlot[np.isfinite(imToPlot)], 100.0),
                        fignum=False, origin='lower', cmap=mpl.cm.hot)
            mpl.colorbar()
            mpl.title('Cleaned+interpolated image')
            mpl.suptitle(plotTitle)
            utils.showzcoord()

            
            #--------------------------------------------------------------------

            if doColdFlagging is True:
                fig = mpl.figure(figsize=(5,5))            
                imToPlot = np.copy(stdFiltImage)
                imToPlot[np.isnan(imToPlot)] = 0
                mpl.matshow(imToPlot, vmax=np.percentile(imToPlot[np.isfinite(imToPlot)], 100.0),
                            fignum=False, origin='lower')  #, cmap=mpl.cm.hot)
                mpl.colorbar()
                mpl.title('Std. Dev. Image')
                mpl.suptitle(plotTitle)

    if not doColdFlagging:
        assert np.sum(coldMask)==0
        assert np.all(mask != pflags.badPixCal['cold'])

    return {'mask':mask, 'image':im, 'medfiltimage':medFiltImage,
            'maxratio':maxRatio, 'diff':diff, 'differr':diffErr, 'niter':iIter + 1}




def findHotPixels(inputFileName=None, obsFile=None, outputFileName=None,
                  paramFile=None, timeStep=1, startTime=0, endTime= -1, fwhm=3.0,
                  boxSize=5, nSigmaHot=3.0, nSigmaCold=2.5, display=False,
                  ds9display=False, dispToPickle=False, weighted=False, fluxWeighted=False,
                  maxIter=5, dispMinPerc=0.0, dispMaxPerc=98.0, diagnosticPlots=False,
                  useLocalStdDev=None, useRawCounts=True, bkgdPercentile=10.0):
    '''
    To find hot (and cold/dead) pixels. This routine is the main code entry point.
    Takes an obs. file as input and outputs an .h5 file containing lists of bad time
    intervals for all the pixels, along with reason indicators for each of those
    intervals. Defaults should be somewhat reasonable for a typical on-sky image.
    
    Note that for calibrations frames where there should be no point sources (e.g.
    laser calibrations, flatfields), can set fwhm=np.inf, and useLocalStdDev=True
    (and/or set nSigma appropriately) in order to just do sigma-clipping without
    bothering to try to account for real astrophysical PSFs. Should be a bit more
    aggressive. Can also set useLocalStdDev=True in this case, which may also help.
    
    NB - AT THE MOMENT, I THINK THIS WILL OVERWRITE PRE-EXISTING HOT PIXEL FILES
    WITHOUT WARNING (should update this behaviour....)
    
    INPUTS:
        inputFileName - string, pathname of input observation file.
        obsFile - instead of providing inputFileName, can optionally directly pass
                  an ObsFile instance here instead (overrides inputFileName). (Added JvE 7/16/2014)
        outputFileName - string, pathname of output .h5 file.
        paramFile - string, name of input parameter file. Other parameters 
                    provided as arguments will override values in the parameter
                    file. See example parameter file .dict in pipeline
                    'params' directory.
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
        nSigmaHot: Scalar float. If the flux ratio for a pixel is (nSigma x expected
                    error) above the max expected given the PSF FWHM, then flag it 
                    as hot. Larger value => lower sensitivity.
        nSigmaCold: Scalar float. Similar to nSigmaHot, but for cold pixels: if flux
                    in a pixel is lower than that in the median filtered image by more
                    than nSigma times the std. dev. in a surrounding box, flag it 
                    as cold.
        display: Boolean. If true, then display the input image and mark those 
                    pixels flagged as hot/cold/dead with a coloured dot. NB - UPDATED
                    TO ONLY SHOW FIRST AND LAST TIME SLICES FOR NOW.
        ds9display: Boolean, as for 'display', but displays the output in ds9 instead.
        dispToPickle: If not False, save whatever would/will be the data for the main plot
                    output by by display=True or ds9display=True to a dictionary in a 
                    pickle file. If dispToPickle is a string, use this as the output
                    file name. Otherwise use a default. See checkInterval for more info.
        weighted: boolean, set to True to use flat cal weights (see flatcal/ 
                    and util.obsFile.getPixelCountImage() )
        fluxWeighted: boolean, if True, flux cal weights are applied (also see fluxcal
                    and util.ObsFile.getPixelCountImage() ). Added JvE 7/16/2014.
        maxIter: Max. number of iterations to do on the bad pixel detection before
                    stopping trying to improve it.
        dispMinPerc: Lower percentile for image stretch if display=True
        dispMaxPerc: Upper percentile "     "     "      "       "
        diagnosticPlots: if True, and display is also True, then show a bunch of
                         diagnostic plots for debug purposes.
        useLocalStdDev - if True, use the local (robust) standard deviation within the 
                    moving box for the sigma value in the hot pixel thresholding
                    instead of Poisson statistics. Mainly intended for situations where
                    you know there is no astrophysical source in the image (e.g. flatfields,
                    laser calibrations), where you can also set fwhm=np.inf
        useRawCounts - Boolean. if True, creates the mask on the basis of the raw photon counts summed
                    across all wavelengths in the obs file. If False, input must be provided
                    in obsFile - in which case, whatever calibrations are already loaded and
                    (and switched on) in the obsFile instance will be applied before looking
                    for badly behaved pixels. *Overrides weighted and fluxWeighted*.
                    Added JvE 7/16/2014
        bkgdPercentile - percentile level (in %) in image to use as an estimate of the background.
                         Added JvE 8/1/2014
                         ***SHOULD BE ADDED AS A PARAMETER TO THE PARAMETER FILE...!!****

        
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
            
        To get a quick view of how the chosen parameters perform, can
        just run on a one-second time-slice and display the results:
    
        findHotPixels(...., startTime=0, endTime=1, display=True)
        
        
    HISTORY:
        2/8/2013: Added iteration to the algorithm, now takes parameter 
                  'maxIter' (also added to parameter file). See
                  checkInterval() for more info.
        5/6/2014: Switched off cold pixel flagging for now - cutoff between
                  genuinely bad cold pixels and those with low QE is not
                  unambiguous enough.
        6/7/2014: Added option to use local box robust estimate of sigma instead
                  of sigma from photon noise estimate.
        8/1/2014: Added proper handling of sky background level in statistics
                  for hot pixel code.  
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
        if nSigmaHot is None: nSigmaHot = params['nSigmaHot']
        if nSigmaCold is None: nSigmaCold = params['nSigmaCold']
        if display is None: display = params['display']
        if maxIter is None: maxIter = params['maxIter']
        if useLocalStdDev is None: useLocalStdDev = params['useLocalStdDev']
    else:
        print 'No parameter file provided - using defaults/input params'
    
    #A few defaults that will be used in the absence of parameter file or 
    #arguments provided by the caller.
    if timeStep is None: timeStep = 1
    if startTime is None: startTime = 0
    if endTime is None: pass
    if maxIter is None: maxIter = 5
    if outputFileName is None: outputFileName = 'badPixTimeMask.h5'
    if useLocalStdDev is None: useLocalStdDev = False
    
    if useRawCounts is False and obsFile is None:
        raise ValueError, ("Must provide obsFile instance if you want to use" +
                            "calibrated (not raw) counts")
        
    if obsFile is None:
        obsFile = ObsFile.ObsFile(inputFileName)
    
    expTime = obsFile.getFromHeader('exptime')    
    if endTime < 0: endTime = expTime
    stepStarts = np.arange(startTime, endTime, timeStep)  #Start time for each step (in seconds).
    stepEnds = stepStarts + timeStep                      #End time for each step
    assert(np.sum(stepEnds>endTime)<=1)                   #Shouldn't be more than one end time that runs over
    stepEnds[stepEnds>endTime] = endTime                  #Clip any time steps that run over the end of the requested time range.
    assert(np.all((stepStarts >= startTime) & (stepStarts <= endTime)))
    assert(np.all((stepEnds >= startTime) & (stepEnds <= endTime)))

    nSteps = len(stepStarts)
    stepStartsTicks = stepStarts * obsFile.ticksPerSec
    stepEndsTicks = stepEnds * obsFile.ticksPerSec
        
    #Initialise stack of masks, one for each time step
    masks = np.zeros([obsFile.nRow, obsFile.nCol, nSteps], dtype=np.int8)

    #Get the mask for each time step
    for i, eachTime in enumerate(stepStarts):
        print 'Processing time slice: ', str(eachTime) + ' - ' + str(eachTime + timeStep) + 's'
        displayThisOne = (display or diagnosticPlots) and (i==0 or i==len(stepStarts)-1)
        ds9ThisOne = ds9display and (i==0 or i==len(stepStarts)-1)
        dispToPickleThisOne = dispToPickle and (i==0 or i==len(stepStarts)-1)
        masks[:, :, i] = checkInterval(obsFile=obsFile, firstSec=eachTime, intTime=timeStep,
                                     fwhm=fwhm, boxSize=boxSize, nSigmaHot=nSigmaHot,
                                     nSigmaCold=nSigmaCold, display=displayThisOne, ds9display=ds9ThisOne, 
                                     dispToPickle=dispToPickleThisOne, weighted=weighted, fluxWeighted=fluxWeighted,
                                     maxIter=maxIter, dispMinPerc=dispMinPerc, dispMaxPerc=dispMaxPerc,
                                     useLocalStdDev=useLocalStdDev, diagnosticPlots=diagnosticPlots 
                                     and displayThisOne, useRawCounts=useRawCounts, bkgdPercentile=bkgdPercentile)['mask']
                                     #Note checkInterval call should automatically clip at end of obsFile,
                                     #so don't need to worry about endTime.
    
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
    writeHotPixels(timeMaskData, obsFile, outputFileName, startTime=startTime, endTime=endTime)
    




def writeHotPixels(timeMaskData, obsFile, outputFileName, startTime=None, endTime=None):
    '''
    Write the output hot-pixel time masks table to an .h5 file. Called by
    findHotPixels().
    
    ** AT THE MOMENT, I THINK THIS WILL OVERWRITE PRE-EXISTING HOT PIXEL FILES
    WITHOUT WARNING - SHOULD UPDATE THIS BEHAVIOUR **
    
    INPUTS:
        timeMaskData - list structure as constructed by findHotPixels()
        obsFile - the ObsFile object from which the data to be written
                        was derived.
        outputFileName - string, the pathname of the .h5 file to write to.
        startTime, endTime - start and end times of the time mask data, in 
                        seconds from the start of the obsFile. These values
                        are written into the header information.
    
    OUTPUTS:
        Writes an .h5 file to filename outputFileName. See findHotPixels()
        for full description of the output data structure.
        
    HISTORY:
        2/15/2013: Updated so that behaviour of outputFileName is consistent with the
        behaviour of the input file name for an ObsFile instance. (i.e., 
        unless the path provided is absolute, $MKID_RAW_PATH is prepended
        to the file name.)
        
        11/22/2013 (actually a few days before): added startTime, and endTime to saved
        header information. Had also added expTime to saved header info some time ago....
    '''
    
    #**** - THIS MECHANISM IS KIND OF REDUNDANT NOW, PRETTY MEANINGLESS *****
    #Don't try and do anything fancy with interpreting the file name here any more.
    #
    #if (os.path.isabs(outputFileName)):
    #    #self.fileName = os.path.basename(fileName)
    #    fullFileName = outputFileName
    #else:
    #    #self.fileName = fileName
    #    # make the full file name by joining the input name 
    #    # to the MKID_RAW_PATH (or . if the environment variable 
    #    # is not defined)
    #    dataDir = os.getenv('MKID_RAW_PATH', '/')
    #    fullFileName = os.path.join(dataDir, outputFileName)
    
    
    fullFileName = os.path.abspath(outputFileName)

    fileh = tables.openFile(fullFileName, mode='w')
    
    try:    
        #Construct groups for header and time mask data.
        headerGroup = fileh.createGroup("/", headerGroupName, 'Header')
        timeMaskGroup = fileh.createGroup("/", dataGroupName, 'Time masks for bad pixels')

        #Fill in header info (just a one-row table)
        headerTable = fileh.createTable(headerGroup, headerTableName, headerDescription,
                                        'Header Info')
        header = headerTable.row
        header[obsFileColName] = obsFile.fileName
        header[nColColName] = max([x[0] for x in timeMaskData]) + 1  #Assume max. x value represents number of columns
        header[nRowColName] = max([x[1] for x in timeMaskData]) + 1  #Same for rows.
        header[ticksPerSecColName] = obsFile.ticksPerSec
        header[expTimeColName] = obsFile.getFromHeader('exptime')    #Newly implemented - SHOULD DOUBLE CHECK! Should automatically account for any of Matt's time corrections JvE 7/08/2013.
        if startTime is not None:
            header[startTimeColName] = startTime
        else:
            header[startTimeColName] = np.nan
        if endTime is not None:
            header[endTimeColName] = endTime
        else:
            header[endTimeColName] = np.nan
            
        header.append()
        headerTable.flush()
        headerTable.close()

        
        #Establish a mapping indicating which bad pixel map flag numbers (uesd in this bad pixel masking
        #code) correspond to which timeMaskReason enumeration values (used by the headers.timeMask 
        #table format in the output files):
        flagMap = {
                   pflags.badPixCal['hot']: tm.timeMaskReason['hot pixel'],
                   pflags.badPixCal['cold']: tm.timeMaskReason['cold pixel'],
                   pflags.badPixCal['dead']: tm.timeMaskReason['dead pixel']
                   }
                  
 
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
                row['reason'] = flagMap.get(eachPeriod[2],'unknown') #Map flag number to time mask reason using flagMap dictionary (default is 'unknown').
                #if eachPeriod[2] == pflags.badPixCal['hot']: 
                #    row['reason'] = tm.timeMaskReason['hot pixel']
                #elif eachPeriod[2] == pflags.badPixCal['cold']:
                #    row['reason'] = tm.timeMaskReason['cold pixel']
                #elif eachPeriod[2] == pflags.badPixCal['dead']:
                #    row['reason'] = tm.timeMaskReason['dead pixel']
                #else:
                #    row['reason'] = tm.timeMaskReason['unknown'] #For now - other reasons (e.g. cold pixels) can be added later.
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
                   photon list file.
    
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
            'expTime' - duration of original obs file (sec).
            'startTime' - start time of the time mask file within the original obs file (sec).
            'endTime' - end time of the time mask file within the original obs file (sec).
            
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
            
    
    NOTES:
        Added expTime, startTime, and endTime return values, JvE, 11/22/2013
        
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
        headerTable = fileh.getNode(nodePath + headerGroupName, headerTableName)
        headerInfo = headerTable[0] #The one row from the header
        nRow = headerInfo[nRowColName]
        nCol = headerInfo[nColColName]
        obsFileName = headerInfo[obsFileColName]
        ticksPerSec = headerInfo[ticksPerSecColName]
        expTime = headerInfo[expTimeColName] if expTimeColName in headerTable.colnames else None 
        startTime = headerInfo[startTimeColName] if startTimeColName in headerTable.colnames else None 
        endTime = headerInfo[endTimeColName] if endTimeColName in headerTable.colnames else None 
            
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
                "obsFileName":obsFileName, "ticksPerSecond":ticksPerSec,
                "expTime":expTime, "startTime":startTime, "endTime":endTime}


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
    
    paramFile = '/Users/vaneyken/ARCONS/pipeline/github/ARCONS-pipeline/params/hotPixels.dict'
    inputFile = '/Users/vaneyken/Data/UCSB/ARCONS/turkDataCopy/ScienceData/PAL2012/20121208/obs_20121209-120530.h5'
    outputFile = None
    
    nArg = len(sys.argv)
    #if nArg == 1: paramFile = 'hotPixels.dict'
    if nArg > 1: paramFile = sys.argv[1]
    if nArg > 2: inputFile = sys.argv[2]
    if nArg > 3: outputFile = sys.argv[3]
    #findHotPixels(paramFile=paramFile, inputFileName=inputFile,
                  #outputFileName=outputFile, endTime=4.5)
    fn = '/Users/vaneyken/Data/UCSB/ARCONS/turkDataCopy/ScienceData/PAL2013/20131209/cal_20131209-120649.h5'
    findHotPixels(fn,outputFileName='hptest.h5',startTime=61, endTime=62,display=True,fwhm=np.inf,useLocalStdDev=True,nSigmaHot=3.0,maxIter=10)
    





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


