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
from scipy.signal import convolve


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
    
