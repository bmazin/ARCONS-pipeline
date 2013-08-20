import unittest
import numpy as np
from beammap import remapPixels

def testRemap():
    '''
    Run a quick test on the beam-remapping code
    '''
    
    #Number of detector rows, columns...
    nRow,nCol = 46,44

    #Arbitrary set of pixels to remap
    #Pixel at row 1,2 goes to 3,4; 10,11 goes to 12,13, etc.:
    remapPix =  [(1,2,3,4),(10,11,12,13),(15,20,25,30),
                 (25,30,1,2),(12,13,10,11)]
    
    #Value to give to unassigned pixels in the remapped result
    missingVal = -1
    
    #Make a list of the pixels involved in the remapping
    oldPixList = np.array(remapPix)[:,0:2]  #Source pixel coordinates
    newPixList = np.array(remapPix)[:,2:4]  #Destination pixel coordinates
    
    #Make the pixel map
    rmp = remapPixels.PixelMap(nRow=nRow,nCol=nCol,remapPix = remapPix)
    
    #Make a fake original array
    oldArr = np.arange(nRow*nCol).reshape((nRow,nCol))
    
    #And remap the original array
    newArr = rmp.remapArray(oldArr,missingVal=missingVal)

    #Check all those pixels which should not have been remapped:
    for iRow in range(nRow):
        for iCol in range(nCol):
            if (iRow,iCol) not in oldPixList and (iRow,iCol) not in newPixList:
                assert oldArr[iRow,iCol]==newArr[iRow,iCol]
    
    #Now check all those pairs that *should* have been remapped:
    assert newArr[3,4]==oldArr[1,2]
    assert newArr[12,13]==oldArr[10,11]
    assert newArr[25,30]==oldArr[15,20]
    assert newArr[1,2]==oldArr[25,30]
    assert newArr[10,11]==oldArr[12,13]
    assert newArr[15,20]==missingVal    
    
    #All good if we got this far
    print 'Yep, all good....'

 
if __name__ == '__main__':
    testRemap()
