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
import pyfits
import matplotlib.pylab as mpl
import photonlist.photlist as pl
import photonlist.RADecImage as rdi
from util.FileName import FileName
from util import utils


def makeRGBset(inputFiles='*.pkl',medIm=True):
    '''
    For now, just turns output virtual image pickle files from makeImageStack into fits files.
    Also, masks wherever there is not valid data in all three images (so that you don't get
    weird coloured pixels popping out).
    All happens in the current working directory.
    '''

    filenames = glob.glob(inputFiles)
    if len(filenames) != 3:
        print 'There should be three input files (R,G,B)...'
        return
    
    nrow,ncol = np.shape(np.load(filenames[0])['vim'].image)
    imstack = np.zeros((nrow,ncol,3))
    for i,eachFilename in enumerate(filenames):
        vim = np.load(eachFilename)
        if medIm is True:
            imstack[:,:,i] = vim['medim']
        else:
            imstack[:,:,i] = vim['vim'].image * vim['vim'].expTimeWeights
    
    #Set any non-finite values to zero
    imstack[~np.isfinite(imstack)]=0    
    mask3d = np.where(imstack>0, 1, 0)    #1 where there's flux, 0 where nothing.
    mask2d = np.prod(mask3d,axis=2)    #Multiply the masks across all three images.
    for i,eachFilename in enumerate(filenames):
        basename,ext = os.path.splitext(os.path.basename(eachFilename))
        outputfile = basename+'.fits'
        print eachFilename + ' -> '+outputfile
        hdu = pyfits.PrimaryHDU(imstack[:,:,i] * mask2d)    #Should zero out anything that doesn't exist in all 3 bands.
        hdu.writeto(outputfile)



def makeImageStack(fileNames='photons_*.h5', dir=os.getenv('MKID_PROC_PATH', 
                   default="/Scratch")+'/photonLists/20121211',
                   detImage=False, saveFileName='stackedImage.pkl', wvlMin=3500,
                   wvlMax=12000, doWeighted=True, medCombine=False, vPlateScale=0.2,
                   nPixRA=250,nPixDec=250,maxBadPixTimeFrac=0.2,integrationTime=-1,
                   outputdir=''):
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
        maxBadPixTimeFrac - Maximum fraction of time which a pixel is allowed to be 
                     flagged as bad (e.g., hot) for before it is written off as
                     permanently bad for the duration of a given image load (i.e., a
                     given obs file).
        integrationTime - the integration time to use from each input obs file (from 
                     start of file).
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
                imSaveName=os.path.join(outputdir,baseSaveName+'det.tif')
                im = phList.getImageDet(wvlMin=wvlMin,wvlMax=wvlMax)
                utils.plotArray(im)
                mpl.imsave(fname=imSaveName,arr=im,colormap=mpl.cm.gnuplot2,origin='lower')
                if eachFile==files[0]:
                    virtualImage=im
                else:
                    virtualImage+=im
            else:
                imSaveName=os.path.join(outputdir,baseSaveName+'.tif')
                virtualImage.loadImage(phList,doStack=not medCombine,savePreStackImage=imSaveName,
                                       wvlMin=wvlMin, wvlMax=wvlMax, doWeighted=doWeighted,
                                       maxBadPixTimeFrac=maxBadPixTimeFrac, integrationTime=integrationTime)
                imageStack.append(virtualImage.image*virtualImage.expTimeWeights)       #Only makes sense if medCombine==True, otherwise will be ignored
                if medCombine==True:
                    medComImage = scipy.stats.nanmedian(np.array(imageStack), axis=0)
                    toDisplay = np.copy(medComImage)
                    toDisplay[~np.isfinite(toDisplay)] = 0
                    utils.plotArray(toDisplay,pclip=0.1,cbar=True,colormap=mpl.cm.gray)
                else:
                    virtualImage.display(pclip=0.5,colormap=mpl.cm.gray)
                    medComImage = None

            mpl.show() 


        else:
            print 'File doesn''t exist: ',eachFile
    
    #Save the results.
    #Note, if median combining, 'vim' will only contain one frame. If not, medComImage will be None.
    results = {'vim':virtualImage,'imstack':imageStack,'medim':medComImage}

    try:
        output = open(os.path(outputdir,saveFileName),'wb')
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
