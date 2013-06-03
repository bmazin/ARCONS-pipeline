'''
Author: Julian van Eyken            Date: May 31 2013
A bit of photon-list image stacking testing....
'''


import photonlist.photlist as pl
import photonlist.RADecImage as rdi
import os.path
import glob

def makeCrabStack():
    dir = '/Users/vaneyken/Data/UCSB/ARCONS/turkDataCopy/Intermediate/photonLists/20121211'
    #files = ['photons_20121206-051516.h5','photons_20121206-052520.h5',
    #         'photons_20121206-052018.h5 ']
    files = glob.glob(os.path.join(dir, '*.h5'))
    
    virtualImage = rdi.RADecImage()
    for eachFile in files:
        print 'Loading: ',os.path.basename(eachFile)
        #fullFileName=os.path.join(dir,eachFile)
        phList = pl.PhotList(eachFile)
        virtualImage.loadImage(phList,stack=True)
        virtualImage.display()
    
    return virtualImage

