
#This code makes a light curve using the photometry modules

import os
import warnings
import numpy as np

from util.FileName import FileName
from util.ObsFile import ObsFile
from util.readDict import readDict
from util.getImages import *
from util.popup import *
from PSFphotometry import PSFphotometry



class LightCurve():

    def __init__(self,path='/Scratch/DisplayStack/RUN_TEMPLATE/TARGET_TEMPLATE',PSF=True,aper=False,**kwargs):
        '''
        Constructs a list of obs FileName objects from the dictionary in the path.
        'run' is assumed to be the second to last directory in the path
        
        Inputs:
            path - path to the display stack target info
            PSF - set True to do PSF fitting
            aper - for aperture photometry
            kwargs - key words for the photometry module

        '''
        
        self.path = path
        self.PSF = PSF
        self.aper = aper
        self.kwargs = kwargs

        
        for f in os.listdir(path):
            if f.endswith(".dict"):
                #if self.verbose: print 'Loading params from ',path+os.sep+f
                self.params = readDict(path+os.sep+f)
                self.params.readFromFile(path+os.sep+f)
                

    def makeLightCurve(self,images,centroids,expTimes=None):
        if expTimes==None:
            expTimes = [None]*len(images)
        flux_list = []
        for i in range(len(images)):
            flux=self.performPhotometry(images[i],centroids[i],expTimes[i])
            flux_list.append(flux)
        return flux_list
                
    def performPhotometry(self,image,centroid,expTime=None):
        '''
        Perform the photometry on an image. 
        
        Input:
            image ----- 2D image of data (0 for dead pixel, shouldn't be any nan's or infs)
                        Should be fully calibrated, dead time corrected, and scaled up to the effective integration time
            centroid -- list of (col,row) tuples. The first tuple is the target location. The next are reference stars in the field
            expTime --- 2D array of pixel exposure times (0 for dead pixels)
                        optional. But can be used for distinguishing 0 count pixels from bad pixels
        '''
        flux = np.zeros(len(centroid))
    
    
        if self.PSF:
            PSFphoto = PSFphotometry(image,centroid,expTime,**self.kwargs)
            flux = PSFphoto.PSFfit(aper_radius=5.)
            
        elif self.aper:
            pass
        else:
            print "Choose PSF fitting or aperture photometry"
            
        return flux
    
    
                
    def getFNs(self,fromObs = True, fromImageStack=False):
        '''
        Should loop through obs files. Return list of file names
        '''
        if fromObs:
            obsFNs = []
            run = os.path.basename(os.path.dirname(self.path))
            for i in range(len(self.params['sunsetDates'])):
                for j in range(len(self.params['obsTimes'][i])):
                    obsFN = FileName(run=run, date=self.params['sunsetDates'][i], tstamp=self.params['utcDates'][i]+'-'+self.params['obsTimes'][i][j])
                    obsFNs.append(obsFN)
                    
            return obsFNs
        elif fromImageStack:
            pass
        
        return None

        
    def getCentroids(self):
        '''
        Should return a list of (col,row) tuples indicating the location of the stars in the field.
        The first location should be the target star. The rest are reference stars.
        
        Needs to be implemented to grab info from centroid file in /Scratch/DisplayStack/
        '''

        return self.testCentroid()
        
    def testCentroid(self):
        target_fn = self.path+os.sep+'CentroidLists'+os.sep+'target'+os.sep+'Centroid_0.h5'
        centroidFile = tables.openFile(target_fn, mode='r')
        xPos = np.array(centroidFile.root.centroidlist.xPositions.read())
        yPos = np.array(centroidFile.root.centroidlist.yPositions.read())
        target_centroids = zip(xPos,yPos)
        centroidFile.close()
    
        ref0_fn = self.path+os.sep+'CentroidLists'+os.sep+'ref0'+os.sep+'Centroid_0.h5'
        centroidFile = tables.openFile(ref0_fn, mode='r')
        xPos = np.array(centroidFile.root.centroidlist.xPositions.read())
        yPos = np.array(centroidFile.root.centroidlist.yPositions.read())
        ref0_centroids = zip(xPos,yPos)
        centroidFile.close()

        centroids = np.asarray(zip(target_centroids,ref0_centroids))
        return centroids




if __name__ == '__main__':
    path = '/Scratch/DisplayStack/PAL2014/HAT_P1'
    verbose=True
    showPlot=True
    
    LC = LightCurve(path,verbose=verbose,showPlot=showPlot)
    #obsFNs = LC.getFNs(fromObs = True, fromImageStack=False)
    ##print 'obs[10]: ',obsFNs[10:11][0].obs()
    #obsFiles = generateObsObjectList(obsFNs,wvlLowerLimit=3500, wvlUpperLimit=5000,loadHotPix=True,loadWvlCal=True,loadFlatCal=True,loadSpectralCal=False)
    #print 'numObs: ',len(obsFiles)
    #im_dict = getImages(fromObsFile=True,fromPhotonList=False,fromImageStack=False,obsFiles=obsFiles, integrationTime=10,weighted=True, fluxWeighted=False, getRawCount=False)
    #images=im_dict['images']
    #pixIntTimes=im_dict['pixIntTimes']
    #print 'numImages: ',len(images)
    
    #writeImageStack(images=images, pixIntTimes=pixIntTimes, startTimes=im_dict['startTimes'], intTimes=im_dict['intTimes'], path=path, outputFilename='ImageStack_1.h5')
    im_dict = getImages(fromObsFile=False,fromPhotonList=False,fromImageStack=True,fullFilename=path+os.sep+'ImageStacks'+os.sep+'ImageStack_0.h5')
    images = im_dict['images']
    pixIntTimes=im_dict['pixIntTimes']
    print 'numLoadedImages: ',len(images)
    centroids = LC.getCentroids()
    print 'numLoadedCentroids: ',len(centroids)
    flux_list=LC.makeLightCurve(images,centroids,pixIntTimes)
    print flux_list
    
    pop(plotFunc=lambda fig,axes: axes.plot(flux_list),title="Flux")
    
    #for i in range(len(images)):
    #    im = images[i]
    #    centroid=LC.getCentroid()
    #    expTime = pixIntTimes[i]
    #    flux=LC.performPhotometry(im,centroid,expTime)
    #    print 'flux ',flux
 





