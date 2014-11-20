import numpy as np
import matplotlib.pyplot as pl
from util.test import loadTestObsFile
from util import utils


def testGetPixelCountImage(bins=250, integrationTime=1):
    '''
    Do two runs of getPixelCountImage and compare the results
    to check the repeatability (i.e., test the degree of 
    effect of the random dithering in the wavelength handling.)
    
    INPUTS:
        bins - set the number of bins for the output histogram
        
    OUTPUTS:
        Displays the two images in ds9, as well as image1 divided
        by image2. Also shows the latter in a regular plot, as well
        as a histogram of the image1/image2 ratios over all pixels.
    '''
    
    obsfile = loadTestObsFile.loadTestObsFile()
    obsfile.setWvlCutoffs(4000,8000)
    
    #Get first image
    im1 = obsfile.getPixelCountImage(firstSec=0, integrationTime=30, weighted=True,
                                     fluxWeighted=False, getRawCount=False, 
                                     scaleByEffInt=False)['image']
    
    #Get supposedly identical image
    im2 = obsfile.getPixelCountImage(firstSec=0, integrationTime=30, weighted=True,
                                     fluxWeighted=False, getRawCount=False,
                                     scaleByEffInt=False)['image']
    
    utils.ds9Array(im1,frame=1)
    utils.ds9Array(im2,frame=2)
    
    divim = im1/im2
    
    utils.ds9Array(divim,frame=3)
    utils.plotArray(divim, colormap=pl.cm.hot, cbar=True, normMax=np.mean(divim)+2*np.std(divim))
    
    toHist = divim.flatten()
    toHist = toHist[np.isfinite(toHist)]
    pl.figure()
    pl.hist(toHist,bins=bins)
    pl.title('Ratio of image1/image2, wavelength range '+str(obsfile.wvlLowerLimit)
             +'-'+str(obsfile.wvlUpperLimit)+'Ang')
    pl.xlabel('Ratio')
    pl.ylabel('Number of pixels')
    
    print 'Mean image1/image2: ',np.mean(toHist)
    print 'Std. dev image1/image2: ',np.std(toHist)
    