
#This code makes a light curve using the photometry modules

import os
import warnings
import numpy as np

from util.FileName import FileName
from util.ObsFile import ObsFile
from util.readDict import readDict
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
        run = os.path.basename(os.path.dirname(path))
        
        for f in os.listdir(path):
            if f.endswith(".dict"):
                #if self.verbose: print 'Loading params from ',path+os.sep+f
                self.params = readDict(path+os.sep+f)
                self.params.readFromFile(path+os.sep+f)
                
        self.obsFNs = []
        for i in range(len(self.params['sunsetDates'])):
            for j in range(len(self.params['obsTimes'][i])):
                obsFN = FileName(run=run, date=self.params['sunsetDates'][i], tstamp=self.params['utcDates'][i]+'-'+self.params['obsTimes'][i][j])
                self.obsFNs.append(obsFN)
                
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
    
    
                
    def getImages(self):
        '''
        Should loop through obs files. Return list of images and exposure times
        '''
        image,expTime = self.getImage()
        return [image],[expTime]


    def getImage(self):
        '''
        This code needs to be implemented
        '''
        return self.testImage()
        
    def testImage(self):
        n = 10
        deadTime = deadTime=100.e-6
        firstSec = 0
        integrationTime = 30
        
        print self.obsFNs[10].obs()
        obs = ObsFile(self.obsFNs[10].obs())
        obs.loadBestWvlCalFile()
        obs.setWvlCutoffs(wvlLowerLimit=3500, wvlUpperLimit=8000)
        obs.loadFlatCalFile(FileName(obsFile=obs).flatSoln())
        obs.loadHotPixCalFile(self.obsFNs[10].timeMask(),reasons=['hot pixel','dead pixel'])
        
        im_dict = obs.getPixelCountImage(firstSec=firstSec, integrationTime=integrationTime,weighted=True, fluxWeighted=False, getRawCount=False)
        im = im_dict['image']
        #Correct for dead time
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",'invalid value encountered in divide',RuntimeWarning)
            #warnings.simplefilter("ignore",RuntimeWarning)
            #warnings.simplefilter("ignore",FutureWarning)
            w_deadTime = 1.0-im_dict['rawCounts']*deadTime/im_dict['effIntTimes']
        im = im/w_deadTime
        #Correct for exposure time
        im = im*integrationTime/im_dict['effIntTimes']
        #Remove any funny values
        im[np.invert(np.isfinite(im))]=0.

        
        return im,im_dict['effIntTimes']
        
    def getCentroid(self):
        '''
        Should return a list of (col,row) tuples indicating the location of the stars in the field.
        The first location should be the target star. The rest are reference stars.
        
        Needs to be implemented to grab info from centroid file in /Scratch/DisplayStack/
        '''

        return self.testCentroid()
        
    def testCentroid(self):
        target = (6.,30.)
        ref = (31.,26.)
        
        #return [target]
        return [target,ref]




if __name__ == '__main__':
    path = '/Scratch/DisplayStack/PAL2014/HAT_P1'
    verbose=True
    showPlot=True
    
    LC = LightCurve(path,verbose=verbose,showPlot=showPlot)
    images,expTimes = LC.getImages()
    
    for i in range(len(images)):
        im = images[i]
        centroid=LC.getCentroid()
        expTime = expTimes[i]
        flux=LC.performPhotometry(im,centroid,expTime)
        print 'flux ',flux
 





