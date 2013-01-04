'''
Created on Dec 12, 2012

@author: vaneyken

Documentation here
'''


from math import *
import numpy as np
import numpy.ma as ma
import matplotlib.pylab as mpl
import util.ObsFile as ObsFile
import util.utils as utils




def findHotPixels(fwhm=3.0, boxsize=5, nsigma=3.0, display=False,
                  firstSec=0, intTime=1, inputFileName='c') :
    '''
    To find the hot pixels in a given time interval for an observation file.
    Compares the ratio of flux in each pixel to the median of the flux in an
    enclosing box. If the ratio is too high -- i.e. the flux is too tightly 
    distributed compared to a Gaussian PSF of the expected FWHM -- then the pixel
    is flagged as bad.
    
        *** Under construction!! ***   JvE Jan 2 2013.
    

    
    For now, main return value is a 2D boolean array of flags corresponding to the
    input image, where 'True'=bad pixel, 'False'=good pixel.
    
    INPUTS:
        fwhm: Scalar float. Estimated full-width-half-max of the PSF (in pixels).
        boxsize: Scalar integer. Size of edge of box used for calculating median 
                 flux in the region surrounding each pixel.
        nsigma: Scalar float. If the flux ratio for a pixel is (nsigma x expected
                error) above the max expected given the PSF FWHM, then flag it as hot.
        display: Boolean. If true, then display the input image and mark those 
                 flagged as hot with a dot.
        inputFileName: String -- input file name...
        firstSec: Scalar integer - start time in seconds within obs. file from 
                  which to integrate when looking for hot pixels.
        intTime: Scalar integer - integration time for hot pixel search (seconds).
 
    OUTPUTS:
        Dictionary containing various diagnostics; main output key is "badflag",
        which contains a 2D array of booleans of the same size as the input image, 
        where True = hot pixel.

    '''
    
    #FirstSec and intTime need to be integers for now.

    #NaN-ignoring median function, to be supplied to median filter later
    #def myMedianFunc(footprint):
    #    return np.median(footprint[~np.isnan(footprint)])
    
    
    
    #Approximate peak/median ratio for an ideal (Gaussian) PSF sampled at 
    #pixel locations corresponding to the median kernal used with the real data.  
    gauss_array = utils.gaussian_psf(fwhm, boxsize)
    maxRatio = np.max(gauss_array) / np.median(gauss_array)

    obsFile = ObsFile.ObsFile(inputFileName)
    #image = np.zeros((obsFile.nRow, obsFile.nCol))
    print 'Counting photons per pixel'
    image = obsFile.getPixelCountImage()
    
    #for i in xrange(obsFile.nRow):
    #    for j in xrange(obsFile.nCol):
    #        image[i, j] = obsFile.getPixelCount(i, j, firstSec, integrationTime=intTime)
    #print 'Done'

    #For now, assume 0 counts in a pixel means the pixel is dead.
    #Make a numpy masked array, and also turn such pixel values into NaNs.
    image = ma.masked_less(image, 1)
    image.data[image.mask] = np.nan
    
    #background = ma.median(image)  #Assume median is a reasonable estimate of the sky background.
    #bsubImage = image - background   #Background subtracted image
    #bsubImageClipped = ma.clip(bsubImage,0,np.Infinity)     #Set all negative values to zero.
    #backgroundSigma1 = (scipy.stats.scoreatpercentile(image.data[~bsubImage.mask], 84.13)
    #                      - scipy.stats.scoreatpercentile(bsubImage.data[~image.mask], 15.87)) / 2.0     #Difference between upper and lower 1-sigma percentiles -- robust estimate of std. dev.
    #absBsubImage = np.abs(bsubImage)
    #backgroundMAD = ma.median(absBsubImage)  #Median absolute deviation (a robust measure of spread).
    #backgroundSigma2 = 1.4826 * backgroundMAD    #See Wikipedia definition of MAD...

    #bsubImageErr = np.sqrt(bsubImageClipped + backgroundSigma ** 2)   #Estimate sigma of flux for each pixel (shot noise + background error)

    #Calculate median filtered image
    #Each pixel takes the median of itself and the surrounding boxsize x boxsize box.
    #(Not sure what edge effects there may be...)
    medFiltImage = utils.median_filterNaN(image.data, boxsize, mode='mirror') #Note - 'reflect' mode looks like it would repeat the edge row/column in the 'reflection'; 'mirror' does not, and makes more sense for this application.

    #Calculate difference between flux in each pixel and maxRatio * the median in the enclosing box.
    #Also calculate the error in the same quantity.
    diff = image - maxRatio * medFiltImage
    diffErr = np.sqrt(image + (maxRatio ** 2) * medFiltImage)

    #Any pixel that has a peak/median ratio more than nsigma above the maximum ratio should be flagged
    #(True=>bad pixel; False=> good pixel).
    badFlag = diff > (nsigma * diffErr)    

    if display:
        imForDisplay = np.copy(image.data)
        imForDisplay[image.mask] = 0  #Just because it looks prettier
        
        #utils.plotArray(image, cbar=True)
        mpl.matshow(imForDisplay)
        mpl.colorbar()
        x = np.arange(np.shape(image)[1])
        y = np.arange(np.shape(image)[0])
        xx, yy = np.meshgrid(x, y)
        mpl.scatter(xx[badFlag], yy[badFlag], c='y')

    return {'badflag':badFlag, 'image':image,
            'medfiltimage':medFiltImage,
            'maxratio':maxRatio,
            'diff':diff, 'differr':diffErr}





if __name__ == "__main__":
    
    findHotPixels()
    
    











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


