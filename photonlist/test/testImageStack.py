'''
Author: Julian van Eyken            Date: May 31 2013
A bit of photon-list image stacking testing....
'''

import os.path
import glob
import matplotlib.pyplot as mpl
import photonlist.photlist as pl
import photonlist.RADecImage as rdi
from util.FileName import FileName
from util import utils


def makeImageStack(fileNames='photons_*.h5',dir='/Scratch/photonLists/20121211',
                   detImage=False):
    '''
    Create an image stack
    INPUTS:
        filenames - string, list of photon-list .h5 files. Can either
                    use wildcards (e.g. 'mydirectory/*.h5') or if string
                    starts with an @, supply a text file which contains
                    a list of file names to stack. (e.g.,
                    'mydirectory/@myfilelist.txt', where myfilelist.txt 
                    is a simple text file with one file name per line.)
        dir - to provide name of a directory in which to find the files
        detImage - if True, show the images in detector x,y coordinates instead
                    of transforming to RA/dec space.
    '''
    
    #Get the list of filenames
    if fileNames[0]=='@':
        #(Note, actually untested, but should be more or less right...)'
        files=[]
        with open(fileNames[1:]) as f:
            for line in f:
                files.append(os.path.join(dir,line.strip()))
    else:
        files = glob.glob(os.path.join(dir, fileNames))

    virtualImage = rdi.RADecImage()
    for eachFile in files:
        if os.path.exists(eachFile):
            print 'Loading: ',os.path.basename(eachFile)
            #fullFileName=os.path.join(dir,eachFile)
            phList = pl.PhotList(eachFile)
            baseSaveName,ext=os.path.splitext(os.path.basename(eachFile))
            
            if detImage is True:
                imSaveName=baseSaveName+'det.tif'
                im = phList.getImageDet()
                utils.plotArray(im)
                mpl.imsave(fname=imSaveName,arr=im,colormap=mpl.cm.gnuplot2,origin='lower')
                if eachFile==files[0]:
                    virtualImage=im
                else:
                    virtualImage+=im
            else:
                imSaveName=baseSaveName+'.tif'
                virtualImage.loadImage(phList,stack=True,savePreStackImage=imSaveName)
                virtualImage.display()
        
        else:
            print 'File doesn''t exist: ',eachFile
            assert 1==0
    
    return virtualImage



def checkRotationDirection():
    '''
    Just create two Crab images which should be rotated by a significant
    amount with respect to each other on the detector, and see if they're 
    properly de-rotated
    '''
    
    plFN1 = FileName(run='PAL2012',date='20121211',tstamp='20121212-033323').photonList()
    plFN2 = FileName(run='PAL2012',date='20121211',tstamp='20121212-045902').photonList()
    
    vIm1 = rdi.RADecImage(pl.PhotList(plFN1))
    vIm2 = rdi.RADecImage(pl.PhotList(plFN2))
    vIm1.display()
    vIm2.display()
    return vIm1,vIm2

