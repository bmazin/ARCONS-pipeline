import time,calendar
import string
import numpy
import scipy
import sys, os
import tables
import pylab
import glob
import matplotlib as mpl
import matplotlib.pylab as plt
import binascii
import math
from scipy.signal import convolve
import scipy.ndimage

"""
Modules:

bin12_9ToRad(binOffset12_9)
confirm(prompt,defaultResponse=True)
convertDegToHex(ra, dec)
convertHexToDeg(ra, dec)
linearFit(x, y, err=None)
plotArray( x, y, z, colormap=mpl.cm.gnuplot2, normMin=None, normMax=None, showMe=True,
              plotFileName='arrayPlot.png', plotTitle='')
printCalFileDescriptions( dir_path )
printObsFileDescriptions( dir_path )

"""


def bin12_9ToRad(binOffset12_9):
   """
   To convert one of the raw 12-bit unsigned values from the photon packet
   into a signed float in radians
   """
   x = binOffset12_9/2.**9. - 4.
   return x


def confirm(prompt,defaultResponse = True):
    """
    Displays a prompt, accepts a yes or no answer, and returns a boolean
    defaultResponse is the response returned if no response is given
    if an ill-formed response is given, the prompt is given again
    """
    if defaultResponse == True:
        optionsString = '[y]|n'
    else:
        optionsString = 'y|[n]'
    goodResponse = False
    while goodResponse == False:
        try:
            responseString = raw_input('%s %s: '%(prompt,optionsString))
            if responseString in ['y','Y','yes','Yes','YES']:
                response = True
                goodResponse = True
            elif responseString in ['n','N','no','No','NO']:
                response = False
                goodResponse = True
            elif responseString == '':
                response = defaultResponse
                goodResponse = True
            else:
                goodResponse = False
        except:
            goodResponse = False
        if goodResponse == False:
            print 'Unrecognized response. Try again.'
    return response
 
def convertDegToHex(ra, dec):
   """
   Convert RA, Dec in decimal degrees to (hh:mm:ss, dd:mm:ss)
   """

   if(ra<0):
      sign = -1
      ra   = -ra
   else:
      sign = 1
      ra   = ra

   h = int( ra/15. )
   ra -= h*15.
   m = int( ra*4.)
   ra -= m/4.
   s = ra*240.

   if(sign == -1):
      outra = '-%02d:%02d:%06.3f'%(h,m,s)
   else: outra = '+%02d:%02d:%06.3f'%(h,m,s)

   if(dec<0):
      sign = -1
      dec  = -dec
   else:
      sign = 1
      dec  = dec

   d = int( dec )
   dec -= d
   dec *= 100.
   m = int( dec*3./5. )
   dec -= m*5./3.
   s = dec*180./5.

   if(sign == -1):
      outdec = '-%02d:%02d:%06.3f'%(d,m,s)
   else: outdec = '+%02d:%02d:%06.3f'%(d,m,s)

   return outra, outdec


def convertHexToDeg(ra, dec):
   """
   Convert RA, Dec in ('hh:mm:ss', 'dd:mm:ss') into floating point degrees.
   """

   try :
      pieces = ra.split(':')
      hh=int(pieces[0])
      mm=int(pieces[1])
      ss=float(pieces[2])
   except:
      raise
   else:
      pass
   
   Csign=dec[0]
   if Csign=='-':
      sign=-1.
      off = 1
   elif Csign=='+':
      sign= 1.
      off = 1
   else:
      sign= 1.
      off = 0

   try :
      parts = dec.split(':')
      deg=int(parts[0][off:len(parts[0])])
      arcmin=int(parts[1])
      arcsec=float(parts[2])
   except:
      raise
   else:
      pass

   return(hh*15.+mm/4.+ss/240., sign*(deg+(arcmin*5./3.+arcsec*5./180.)/100.) )

   
def linearFit( x, y, err=None ):
    """
    Fit a linear function y as a function of x.  Optional parameter err is the
    vector of standard errors (in the y direction).

    Returns:  solution - where y = solution[0] + solution[1]*x
    """
    x = numpy.copy(x)
    y = numpy.copy(y)
    N = len(x)
    A = numpy.ones((2, N), x.dtype)
    A[1] = x
    if err!=None: A /= err
    A = numpy.transpose(A)
    if err!=None: y /= err

    solution, residuals, rank, s = scipy.linalg.lstsq(A, y)
    return solution


def makeMovie( listOfFrameObj, frameTitles=None, outName='Test_movie',
              delay=0.1, **plotArrayKeys):
    """
    Makes a movie out of a list of frame objects (2-D arrays)
    
    """
    if len(listOfFrameObj) == 1:
        raise ValueError, "I cannot make movie out of a list of one object!"

    if frameTitles != None:
        assert len(frameTitles) == len(listOfFrameObj), "Number of Frame titles\
        must equal number of frames"

    if os.path.exists("./.tmp_movie"):
        os.system("rm -rf .tmp_movie")

    os.mkdir(".tmp_movie")
    iFrame = 0
    print 'Making individual frames ...'
    
    for frame in listOfFrameObj:

       if frameTitles!= None:
          plotTitle = frameTitles[iFrame]
       else:
          plotTitle=''
       fp = plotArray(frame, showMe=False, plotFileName='.tmp_movie/mov_'+repr(iFrame+10000)+'.png', \
                      plotTitle=plotTitle, **plotArrayKeys)
       iFrame += 1
       del fp

    os.chdir('.tmp_movie')

    if outName[-4:-1]+outName[-1] != '.gif':
        outName += '.gif'

    delay *= 100
    delay = int(delay)
    print 'Making Movie ...'

    if '/' in outName:
        os.system('convert -delay %s -loop 0 mov_* %s'%(repr(delay),outName))
    else:
        os.system('convert -delay %s -loop 0 mov_* ../%s'%(repr(delay),outName))
    os.chdir("../")
    os.system("rm -rf .tmp_movie")
    print 'done.'



def plotArray( xyarray, colormap=mpl.cm.gnuplot2, normMin=None, normMax=None, showMe=True,
              cbar=False, cbarticks=None, cbarlabels=None, plotFileName='arrayPlot.png',
              plotTitle='', sigma=None):
    """
    Plots the 2D array to screen or if showMe is set to False, to file.  If normMin and
    normMax are None, the norm is just set to the full range of the array.
    """
    if sigma != None:
       meanVal = np.mean(accumulatePositive(xyarray))
       stdVal = np.std(accumulatePositive(xyarray))
       normMin = meanVal - sigma*stdVal
       normMax = meanVal + sigma*stdVal
    if normMin == None:
       normMin = xyarray.min()
    if normMax == None:
       normMax = xyarray.max()
    norm = mpl.colors.Normalize(vmin=normMin,vmax=normMax)

    figWidthPt = 550.0
    inchesPerPt = 1.0/72.27                 # Convert pt to inch
    figWidth = figWidthPt*inchesPerPt       # width in inches
    figHeight = figWidth*1.0                # height in inches
    figSize =  [figWidth,figHeight]
    params = {'backend': 'ps',
              'axes.labelsize': 10,
              'axes.titlesize': 12,
              'text.fontsize': 10,
              'legend.fontsize': 10,
              'xtick.labelsize': 10,
              'ytick.labelsize': 10,
              'figure.figsize': figSize}
    plt.rcParams.update(params)

    plt.matshow(xyarray, cmap=colormap, origin='lower',norm=norm)

    if cbar:
        if cbarticks == None:
           cbar = plt.colorbar(shrink=0.8)
        else:
           cbar = plt.colorbar(ticks=cbarticks, shrink=0.8)
        if cbarlabels != None:
           cbar.ax.set_yticklabels(cbarlabels)
    
    plt.ylabel('Row Number')
    plt.xlabel('Column Number')
    plt.title(plotTitle)

    if showMe == False:
        plt.savefig(plotFileName)
    else:    
        plt.show()
 

def printCalFileDescriptions( dir_path ):
    """
    Prints the 'description' and 'target' header values for all calibration
    files in the specified directory
    """
    for obs in glob.glob(os.path.join(dir_path,'cal*.h5')):
       f=tables.openFile(obs,'r')
       hdr=f.root.header.header.read()
       print obs,hdr['description'][0]
       target = f.root.header.header.col('target')[0]
       print target
       f.close()
    

def printObsFileDescriptions( dir_path ):
    """
    Prints the 'description' and 'target' header values for all observation
    files in the specified directory
    """
    for obs in glob.glob(os.path.join(dir_path,'obs*.h5')):
       f=tables.openFile(obs,'r')
       hdr=f.root.header.header.read()
       print obs,hdr['description'][0]
       target = f.root.header.header.col('target')[0]
       print target
       f.close()
  

def median_filterNaN(inputarray, size=5, *nkwarg, **kwarg):
    '''
    NaN-handling version of the scipy median filter function
    (scipy.ndimage.filters.median_filter). Any NaN values in the input array are
    simply ignored in calculating medians. Useful e.g. for filtering 'salt and pepper
    noise' (e.g. hot/dead pixels) from an image to make things clearer visually.
    (but note that quantitative applications are probably limited.)
    Works as a simple wrapper for scipy.ndimage.filters.generic-filter, to which
    calling arguments are passed.
    
    Arguments/return values are same as for scipy median_filter.
    INPUTS:
        inputarray : array-like, input array to filter (can be n-dimensional)
        size : scalar or tuple, optional, size of edge(s) of n-dimensional moving box. If 
                scalar, then same value is used for all dimensions.
    OUTPUTS:
        NaN-resistant median filtered version of inputarray.
    
    For other parameters see documentation for scipy.ndimage.filters.median_filter.

    e.g.:
        
        filteredImage = median_filterNaN(imageArray,size=3)
    
    -- returns median boxcar filtered image with a moving box size 3x3 pixels.
    
    JvE 12/28/12
    '''     
    return scipy.ndimage.filters.generic_filter(inputarray, lambda x:numpy.median(x[~numpy.isnan(x)]), size,
                                                 *nkwarg, **kwarg)
    
    
    
def mean_filterNaN(inputarray, size=3, *nkwarg, **kwarg):
    '''
    Basically a box-car smoothing filter. Same as median_filterNaN, but calculates a mean instead. 
    Any NaN values in the input array are ignored in calculating means.
    See median_filterNaN for details.
    JvE 1/4/13
    '''
         
    return scipy.ndimage.filters.generic_filter(inputarray, lambda x:numpy.mean(x[~numpy.isnan(x)]), size,
                                                 *nkwarg, **kwarg)
    

def replaceNaN(inputarray, mode='mean', boxsize=3, iterate=True):
    '''
    Replace all NaN values in an array with the mean (or median)
    of the surrounding pixels. Should work for any number of dimensions, 
    but only fully tested for 2D arrays at the moment.
    
    INPUTS:
        inputarray - input array
        mode - 'mean' or 'median' to replace with the mean or median of the neighbouring pixels.
        boxsize - scalar integer, length of edge of box surrounding bad pixels from which to
                  calculate the mean or median.
        iterate - If iterate is set to True then iterate until there are no NaN values left.
                  (To deal with cases where there are many adjacent NaN's, where some NaN
                  elements may not have any valid neighbours to calculate a mean/median. 
                  Such elements will remain NaN if only a single pass is done.)
    OUTPUTS:
        Returns 'inputarray' with NaN values replaced.
        
    TO DO: currently spits out multiple 'invalid value encoutered' warnings if 
           NaNs are not all removed on the first pass. These can safely be ignored.
           Will implement some warning catching to suppress them.
    JvE 1/4/2013    
    '''
    
    outputarray = numpy.copy(inputarray)
    print numpy.sum(numpy.isnan(outputarray))
    while numpy.sum(numpy.isnan(outputarray)) > 0:
        
        #Calculate interpolates at *all* locations (because it's easier...)
        if mode=='mean':
            interpolates=mean_filterNaN(outputarray,size=boxsize,mode='mirror')
        elif mode=='median':
            interpolates=median_filterNaN(outputarray,size=boxsize,mode='mirror')
        else:
            raise ValueError('Invalid mode selection - should be one of "mean" or "median"')
        
        #Then substitute those values in wherever there are NaN values.
        outputarray[numpy.isnan(outputarray)] = interpolates[numpy.isnan(outputarray)]
        print numpy.sum(numpy.isnan(outputarray))
        if not iterate: break 

    return outputarray
    
    
    
def rebin2D(a, ysize, xsize):
    '''
    Rebin an array to a SMALLER array. Rescales the values such that each element
    in the output array is the mean of the elememts which it encloses in the input
    array (i.e., not the total). Similar to the IDL rebin function.
    Dimensions of binned array must be an integer factor of the input array.
    Adapted from SciPy cookbook - see http://www.scipy.org/Cookbook/Rebinning
    JvE 12/28/12

    INPUTS:
        a - array to be rebinned
        ysize - new ysize (must be integer factor of y-size of input array)
        xsize - new xsize (ditto for x-size of input array)

    OUTPUTS:
        Returns the original array rebinned to the new dimensions requested.        
    '''
    
    yfactor, xfactor = numpy.asarray(a.shape) / numpy.array([ysize, xsize])
    return a.reshape(ysize, yfactor, xsize, xfactor,).mean(1).mean(2)

  
  
def gaussian_psf(fwhm, boxsize, oversample=50):
    
    '''
    Returns a simulated Gaussian PSF: an array containing a 2D Gaussian function
    of width fwhm (in pixels), binned down to the requested box size. 
    JvE 12/28/12
    
    INPUTS:
        fwhm - full-width half-max of the Gaussian in pixels
        boxsize - size of (square) output array
        oversample (optional) - factor by which the raw (unbinned) model Gaussian
                                oversamples the final requested boxsize.
    
    OUTPUTS:
        2D boxsize x boxsize array containing the binned Gaussian PSF
    
    (Verified against IDL astro library daoerf routine)
        
    '''
  
    fineboxsize = boxsize * oversample
    
    xcoord = ycoord = numpy.arange(-(fineboxsize - 1.) / 2., (fineboxsize - 1.) / 2. + 1.)
    xx, yy = numpy.meshgrid(xcoord, ycoord)
    xsigma = ysigma = fwhm / (2.*math.sqrt(2.*math.log(2.))) * oversample
    zx = (xx ** 2 / (2 * xsigma ** 2))
    zy = (yy ** 2 / (2 * ysigma ** 2))
    fineSampledGaussian = numpy.exp(-(zx + zy))

    #Bin down to the required output boxsize:
    binnedGaussian = rebin2D(fineSampledGaussian, boxsize, boxsize)

    return binnedGaussian
        
