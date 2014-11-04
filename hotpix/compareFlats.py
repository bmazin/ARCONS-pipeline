'''
For investigating the comparison of flats of different flux
levels in order to find badly behaved pixels. Not part of the
main pipeline for now, just for testing/development purposes.
'''

import os.path
import numpy as np
import astropy.stats
import matplotlib.pyplot as mpl
from util import ObsFile,utils

def divideObsFiles(fileName1='/Users/vaneyken/Data/UCSB/ARCONS/turkDataCopy/ScienceData/PAL2012/20121211/flat_20121212-134024.h5',
                   fileName2='/Users/vaneyken/Data/UCSB/ARCONS/turkDataCopy/ScienceData/PAL2012/20121211/flat_20121212-134637.h5',
                   firstSec=0, integrationTime=10., nbins=None):
    
    '''
    Divide the raw image from one obs file by another, and display and return the result.
    This works with raw counts, so the obs file need not be calibrated.
    
    INPUTS:
        fileName1, fileName2 -- names of two obs files to divide by each other, OR two ObsFile
                                objects can be passed directly.
        firstSec - time during the obs files at which to start integration of images (sec)
        integrationTime - time to integrate for to make the images (sec)
        nbins - number of bins for histogram plot (if None, makes a semi-reasonable guess)
        
    OUTPUTS:
        Displays the divided result to the screen and to ds9.
        Returs a tuple of image arrays:
            divided image, input image 1, input image 2, obsFile 1, obsFile 2
    '''

    if type(fileName1) is str:
        obsf1 = ObsFile.ObsFile(fileName1)
        fn1 = fileName1
    else:
        obsf1=fileName1       #Otherwise just assume it's an ObsFile instance.
        fn1=obsf1.fileName
        
    if type(fileName2) is str:
        obsf2 = ObsFile.ObsFile(fileName2)
        fn2 = fileName2
    else:
        obsf2=fileName2
        fn2=obsf2.fileName
    
    print 'Reading '+os.path.basename(fn1)
    pci1 = obsf1.getPixelCountImage(firstSec=firstSec,integrationTime=integrationTime,getRawCount=True)
    print 'Reading '+os.path.basename(fn2)
    pci2 = obsf2.getPixelCountImage(firstSec=firstSec,integrationTime=integrationTime,getRawCount=True)
    
    im1 = pci1['image']
    im2 = pci2['image']
    divIm = im1/im2
    
    med = np.median(divIm[~np.isnan(divIm)])
    #Approximate std. dev from median abs. dev. (more robust)
    sdev = astropy.stats.median_absolute_deviation(divIm[~np.isnan(divIm)]) *1.4826 
    badFlag = np.abs(divIm - med) > 3.0*sdev
    
    print 'Displaying'
    toDisplay = np.copy(divIm)
    toDisplay[np.isnan(toDisplay)]=0
    utils.plotArray(toDisplay, cbar=True, normMin=med-4.*sdev, normMax=med+4.*sdev, colormap=mpl.cm.hot, fignum=None)
    yy,xx = np.indices(np.shape(divIm))
    mpl.scatter(xx[badFlag], yy[badFlag], marker='o', c='r')
    utils.ds9Array(divIm)
    mpl.figure()
    if nbins is None: nbins=np.sqrt(np.sum(np.isfinite(divIm)))
    mpl.hist(divIm[np.isfinite(divIm)].flatten(),bins=nbins)
    
    print 'Median: ',med
    print 'Approx std. dev. (M.A.D * 1.4826): ',sdev
    print 'Done.'
    return divIm,im1,im2,obsf1,obsf2,badFlag


    
    