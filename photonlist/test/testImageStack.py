'''
Author: Julian van Eyken            Date: May 31 2013
A bit of photon-list image stacking testing....
'''

import warnings
import pickle
import os.path
import glob
#from astropy import coordinates as coord
import matplotlib.pyplot as mpl
import photonlist.photlist as pl
import photonlist.RADecImage as rdi
from util.FileName import FileName
from util import utils


def makeImageStack(fileNames='photons_*.h5', dir=os.getenv('INTERM_DIR', default="/Scratch")+'/photonLists/20121211',
                   detImage=False, saveFileName='stackedImage.pkl', wvlMin=None,
                   wvlMax=None, doWeighted=True):
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
        saveFileName - name of output pickle file for saving final resulting object.
        doWeighted - boolean, if True, do the image flatfield weighting.
    
    OUTPUTS:
        Returns a stacked image object, saves the same out to a pickle file, and
        (depending whether it's still set to or not) saves out the individual non-
        stacked images as it goes. 
    '''
    
    #Get the list of filenames
    if fileNames[0]=='@':
        #(Note, actually untested, but should be more or less right...)
        files=[]
        with open(fileNames[1:]) as f:
            for line in f:
                files.append(os.path.join(dir,line.strip()))
    else:
        files = glob.glob(os.path.join(dir, fileNames))

    #Initialise empty image centered on Crab Pulsar
    virtualImage = rdi.RADecImage(nPixRA=500,nPixDec=500,vPlateScale=0.1,
                                  cenRA=1.4596725441339724, cenDec=0.38422539085925933)
                                  
    for eachFile in files:
        if os.path.exists(eachFile):
            print 'Loading: ',os.path.basename(eachFile)
            #fullFileName=os.path.join(dir,eachFile)
            phList = pl.PhotList(eachFile)
            baseSaveName,ext=os.path.splitext(os.path.basename(eachFile))
            
            if detImage is True:
                imSaveName=baseSaveName+'det.tif'
                im = phList.getImageDet(wvlMin=wvlMin,wvlMax=wvlMax)
                utils.plotArray(im)
                mpl.imsave(fname=imSaveName,arr=im,colormap=mpl.cm.gnuplot2,origin='lower')
                if eachFile==files[0]:
                    virtualImage=im
                else:
                    virtualImage+=im
            else:
                imSaveName=baseSaveName+'.tif'
                virtualImage.loadImage(phList,doStack=True,savePreStackImage=imSaveName,
                                       wvlMin=wvlMin, wvlMax=wvlMax, doWeighted=doWeighted)
                virtualImage.display(pclip=0.1)
        
        else:
            print 'File doesn''t exist: ',eachFile
    
    #Save the results
    try:
        output = open(saveFileName,'wb')
        pickle.dump(virtualImage,output,-1)
        output.close()
    except:
        warnings.warn('Unable to save results for some reason...')
    
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

