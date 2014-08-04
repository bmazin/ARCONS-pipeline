'''
Author: Julian van Eyken            Date: May 31 2013
A bit of photon-list image stacking testing....
NB - really designed for Crab pulsar image, but should
be somewhat flexible.
'''

import warnings
import pickle
import os.path
import glob
import scipy.stats
import numpy as np
#from astropy import coordinates as coord
import matplotlib.pyplot as mpl
import photonlist.photlist as pl
import photonlist.RADecImage as rdi
from util.FileName import FileName
from util import utils


def makeImageStack(fileNames='photons_*.h5', dir=os.getenv('MKID_PROC_PATH', default="/Scratch")+'/photonLists/20121211',
                   detImage=False, saveFileName='stackedImage.pkl', wvlMin=None,
                   wvlMax=None, doWeighted=True, medCombine=False, vPlateScale=0.2,
                   nPixRA=250,nPixDec=250):
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
        medCombine - experimental, if True, do a median combine of the image stack
                     instead of just adding them all.... Prob. should be implemented
                     properly at some point, just a fudge for now.
        vPlateScale - (arcsec/virtual pixel) - to set the plate scale of the virtual
                     pixels in the outputs image.
        nPixRA,nPixDec - size of virtual pixel grid in output image.
    
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
    virtualImage = rdi.RADecImage(nPixRA=nPixRA,nPixDec=nPixDec,vPlateScale=vPlateScale,
                                  cenRA=1.4596725441339724, cenDec=0.38422539085925933)
    imageStack = []
                                  
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
                virtualImage.loadImage(phList,doStack=not medCombine,savePreStackImage=imSaveName,
                                       wvlMin=wvlMin, wvlMax=wvlMax, doWeighted=doWeighted)
                imageStack.append(virtualImage.image*virtualImage.expTimeWeights)       #Only makes sense if medCombine==True, otherwise will be ignored
                if medCombine==True:
                    medComImage = scipy.stats.nanmedian(np.array(imageStack), axis=0)
                    normMin = np.percentile(medComImage[np.isfinite(medComImage)],q=0.1)
                    normMax = np.percentile(medComImage[np.isfinite(medComImage)],q=99.9)
                    toDisplay = np.copy(medComImage)
                    toDisplay[~np.isfinite(toDisplay)] = 0
                    utils.plotArray(toDisplay,normMin=normMin,normMax=normMax,colormap=mpl.cm.gray,
                                    cbar=True)
                else:
                    virtualImage.display(pclip=0.1)
                    medComImage = None
                        
        else:
            print 'File doesn''t exist: ',eachFile
    
    #Save the results.
    #Note, if median combining, 'vim' will only contain one frame. If not, medComImage will be None.
    results = {'vim':virtualImage,'imstack':imageStack,'medim':medComImage}

    try:
        output = open(saveFileName,'wb')
        pickle.dump(results,output,-1)
        output.close()
            
    except:
        warnings.warn('Unable to save results for some reason...')
    
    return results



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


if __name__ == '__main__':
    makeImageStack()
