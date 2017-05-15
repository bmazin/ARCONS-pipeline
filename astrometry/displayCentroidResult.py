import numpy as np
import tables
import os
from util import ObsFile, FileName, utils


def displayCentroidResult(obsFileIn, time):
    '''
    To show an image with the location of the centroid
    measured by CentroidCalc marked on top. Use for diagnostic purposes.
    
    INPUTS:
        obsFileIn - either an ObsFile instance or the filename of an 
                  obsFile to load. If a filename, the file will be closed
                  on completion; if an ObsFile instance, it'll be left 
                  alone.
        time - time since beginning of file (in seconds) to display the 
                centroid for.
    
    OUTPUTS:
        A reconstructed image with the calculated centroid plotted on top.
        The image will be integrated over whatever time slice was used
        by CentroidCalc for calculating the centroid at the given input time.
    '''
    
    if type(obsFileIn)=='str':
        obsFile = ObsFile.ObsFile(obsFileIn)
    else:
        obsFile = obsFileIn
    
    ctrdFileName = FileName.FileName(obsFile=obsFile).centroidList()
    ctrdFile = tables.openFile(ctrdFileName, mode='r')

    #Get the boundaries of the time slices used for centroiding
    #(in seconds from start of array.)    
    sliceTimes = ctrdFile.root.centroidlist.times.read()
    xPositions = ctrdFile.root.centroidlist.xPositions.read()
    yPositions = ctrdFile.root.centroidlist.yPositions.read()    
    iSliceEnd = np.searchsorted(sliceTimes, time)
    iSliceStart = iSliceEnd-1
    sliceStartTime = sliceTimes[iSliceStart]
    sliceEndTime = sliceTimes[iSliceEnd]
    sliceXpos,sliceYpos = xPositions[iSliceStart],yPositions[iSliceStart]
    
    #Integrate to get the corresponding image
    im = obsFile.getPixelCountImage(sliceStartTime, sliceEndTime-sliceStartTime,
                                    getRawCount=True)
    
    #And display the result....
    utils.plotArray(im, sigma=3, cbar=True, plotTitle=os.path.basename(ctrdFileName)+', '
                    +str(sliceStartTime)+' - '+ str(sliceEndTime)+'sec', 
                    pixelsToMark=[sliceXpos,sliceYpos], fignum=None)
    
    
