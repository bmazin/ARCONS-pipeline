'''
Author: Julian van Eyken                    Date: Aug 16 2013
Tools for remapping pixels to correct for errors in the
original beam mapping.
'''

import numpy as np
import tables
from util.FileName import FileName

class PixelMap(object): 
    '''
    An object, basically just an array, that is used to remap
    pixels from one location to another in the event of beammapping
    errors.
    '''    
    def __init__(self,fileName=None,nRow=None,nCol=None,remapPix=None):
        '''
        Load up an .h5 file that contains a pixel 'remap'; or
        if a filename is not given, then call makeMap() to 
        create the remap instead. .h5 file format is as saved by
        function makeMap() below.
        
        INPUTS:
            fileName - if set, use load remap information from .h5
                       file of this name.
            nRow, nCol, remapPix - if all set (and not fileName), then
                        pass these all onto makeMap to automatically 
                        create a new pixel map instance without creating a
                        remap file.
        '''
        
        self.pixMap = None
        if fileName is None:
            if nRow is None or nCol is None or remapPix is None:
                raise ValueError, 'Either fileName or all of nRow,nCol,remapPix must be supplied'
            self.pixMap = makeMap(outputFileName=None,nRow=nRow,nCol=nCol,remapPix=remapPix)
        else:
            print 'Reading from file '+fileName
            h5file = tables.openFile(fileName, 'r')
            try:
                self.pixMap = h5file.root.pixMap.read()
            finally:
                h5file.close()


    def remapPix(self,yPix,xPix):
        '''
        Return a two-element array [yNew,xNew] representing the pixel to
        which a given input pixel at yPix,xPix should be mapped.
        '''
        assert self.pixMap is not None
        return (self.pixMap[:,yPix,xPix])


    def remapArray(self,inputArr,missingVal=None):
        '''
        Remap the elements of a 2D array representing data from the detector,
        according to the pixel remap object. Doubtless this function
        can be made faster, but this way should be fine for now....
        
        INPUTS:
            inputArr - array which needs pixels remapping. Size/shape must
                       match the beam-remap instance.
            missingVal - pixels which are mapped *to* somewhere else but which
                       have no source mapped to them *from* another place will
                       take this value (if None, then -999999 is used for integer
                       arrays, NaN for float arrays).
        OUTPUTS:
            Returns corrected version of the input array.
        '''
        
        if np.shape(inputArr) != np.shape(self.pixMap[0,:,:]):
            raise ValueError, 'Shape of array does not match shape of pixel remap instance.'
        if missingVal is None:
            if np.issubdtype(inputArr.dtype,int):
                mVal = -999999
            else:
                mVal = np.nan
        else:
            mVal = missingVal
        
        newArr = np.zeros_like(inputArr)
        newArr.fill(mVal)   #Set all to missing value to begin with.
        nRow,nCol = np.shape(self.pixMap[0,:,:])
        for oldRow in range(nRow):
            for oldCol in range(nCol):
                newRow,newCol = self.remapPix(oldRow,oldCol)
                if newRow>=0 and newCol>=0:
                    newArr[newRow,newCol] = inputArr[oldRow,oldCol]
        return newArr
        
        
        
def makeMap(outputFileName=FileName(run='PAL2012').pixRemap(),
             nRow=46, nCol=44,
             remapPix = [(31,29,3,17), (3,17,31,29),
                         (29,33,42,7), (42,7,32,35),
                         (32,35,40,33), (40,33,29,33)]):
    '''
    Function to make an array mapping from current pixel coordinates
    to desired pixel coordinates, and save it in an HDF file for use by the 
    remapping function.
    
    Defaults are current best guess at corrections for the Palomar 2012 run.
    
    INPUTS:
        outputFileName: name of the HDF file to save out to.
        nRow, nCol: integers, number of pixel rows/columns
        remapPix: list of 4-element tuples of the form (y1,x1,y2,x2), 
                 representing pairs of pixel coordinates. In the 
                 resulting map, pixel x1,y1 will be remapped to pixel
                 x2,y2. Any number of such pairs can be entered.
                 
    OUTPUTS:
        Returns a 3D array, and writes the same to an HDF file if outputFileName
        is not None. Array gives pixel mapping, such that
        
        makeMap()[:,y1,x1] == [y2,x2],
        
        where makeMap() is a 2 x nRow x nCol array, y1,x1 is the current location
        of a given pixel, and y2,x2 is where that pixel *should* be located. Exactly
        the same array format is saved in the HDF file.
        
        If a destination pixel is mapped from somewhere, but that pixel itself is not
        mapped *to* somewhere else, then it becomes unassigned, so that it does not
        end up with a 2-to-1 mapping to itself (the second being itself mapping to itself).
        In this case, coordinates returned are as below in 'unassignedCoord': a large
        enough -ve number that it should hopefully throw an error if you try to 
        index an image array with it.) 
    '''
    
    unassignedCoord = [-99999999,-99999999]
    
    #Make 3D array (2 x nRow x nCol) with values filled out as if all pixels
    #were mapped directly to themselves (i.e., no changes).
    mapArray = np.mgrid[0:nRow,0:nCol]
    sourcePixList = np.array(remapPix)[:,0:2]    #Source coordinates of all remap pixel pairs
    destPixList = np.array(remapPix)[:,2:4]      #Destination coordinates of all remap pixel pairs
    
    #Update pixel map for all the requested remap pairs.
    for sourcePix,destPix in zip(sourcePixList,destPixList):
        mapArray[:,sourcePix[0],sourcePix[1]] = destPix[0:2]
        if destPix not in sourcePixList:
            #any pixels that were mapped *to* but not *from* so that
            #they point nowhere (to avoid a 2-to-1 mapping to that pixel both from its
            #source pixel and itself).
            mapArray[:,destPix[0],destPix[1]] = unassignedCoord     #-ve value means unassigned.
            
    #... and that should do it.
    if outputFileName is not None:
        print 'Writing pixel map to file '+outputFileName
        h5file = tables.openFile(outputFileName, mode='w', title='Pixel remapping map')
        try:
            h5file.createArray('/','pixMap',mapArray,'Mapping from current pixel locations to desired pixel locations')
        finally:
            h5file.close()
    
    return mapArray


if __name__ == '__main__':
    mapArr = makeMap(outputFileName=None, remapPix=[(1,2,3,4),(10,11,12,13),(15,20,25,30)])
    print mapArr[:,1,2]
    print mapArr[:,3,4]
    print