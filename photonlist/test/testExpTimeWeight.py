'''
Author: Julian van Eyken        Date: Jul 9 2013
Test code for per-pixel effective exposure time weighting.
'''

import os.path
import photonlist.photlist as pl
import photonlist.RADecImage as rdi
from util.FileName import FileName

def getSimpleImage(fileName=FileName(run='PAL2012',date='20121211',tstamp='20121212-033323').photonList(),
                   firstSec=0, integrationTime=5):
    '''
    Get a simple short-exposure time RA/dec-mapped image, for
    the purposes of looking at the per-pixel effective integration
    time weighting.
    '''
    
    virtualImage = rdi.RADecImage(vPlateScale=0.1)
    print 'Loading: ',os.path.basename(fileName)
    phList = pl.PhotList(fileName)
    baseSaveName,ext=os.path.splitext(os.path.basename(fileName))
    imSaveName=baseSaveName+'.tif'
    virtualImage.loadImage(phList,stack=True,savePreStackImage=imSaveName,
                           firstSec=firstSec,integrationTime=integrationTime)
    virtualImage.display(pclip=True)

    return virtualImage



if __name__ == "__main__":
    getSimpleImage()