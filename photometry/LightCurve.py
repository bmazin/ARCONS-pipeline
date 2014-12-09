'''
Author: Alex Walter
Date: Dec 3, 2014

This code makes a light curve using the photometry modules
'''

import os
import warnings
import numpy as np

from util.FileName import FileName
from util.ObsFile import ObsFile
from util.readDict import readDict
from util.getImages import *
from util.popup import *
from photometry.PSFphotometry import PSFphotometry, writePhotometryFile, readPhotometryFile



class LightCurve():

    def __init__(self,path='/Scratch/DisplayStack/RUN_TEMPLATE/TARGET_TEMPLATE',fileID='test',PSF=True,aper=False,**kwargs):
        '''
        Constructs a list of obs FileName objects from the dictionary in the path.
        'run' is assumed to be the second to last directory in the path
        
        Inputs:
            path - path to the display stack target info
            fileID - identifier in filename. eg centroid_fileID.h5 or ImageStack_fileID.h5
            PSF - set True to do PSF fitting
            aper - for aperture photometry
            kwargs - key words for the photometry module

        '''
        
        self.path = path
        self.fileID=fileID
        self.PSF = PSF
        self.aper = aper
        self.kwargs = kwargs
        
        #if self.PSF:
        #    from photometry.PSFphotometry import PSFphotometry, writePhotometryFile
        #    print PSFphotometry
        #else:
        #    #import aperture photometry module
        #    pass

        
        for f in os.listdir(path):
            if f.endswith(".dict"):
                #if self.verbose: print 'Loading params from ',path+os.sep+f
                self.params = readDict(path+os.sep+f)
                self.params.readFromFile(path+os.sep+f)
                
    def getImages(self,fromObsFile=False,fromPhotonList=False,fromImageStack=False,**kwargs):
        assert 1.0*fromObsFile+fromPhotonList+fromImageStack==1, "Choose whether to get images from a list of ObsFiles, from a list of PhotonLists, or from an .h5 file made by imagestack"
        
        if fromImageStack:
            imStack_fn=self.path+os.sep+'ImageStacks'+os.sep+'ImageStack_'+self.fileID+'.h5'
            im_dict=getImages(fromImageStack=fromImageStack,fullFilename = imStack_fn)
        elif fromObsFile:
            print 'Not implemented yet!'
            raise IOError
            obsFNs = self.getFNs()
            obsFiles = generateObsObjectList(obsFNs,**kwargs)
            im_dict = getImages(fromObsFile=fromObsFile,obsFiles=obsFiles, integrationTime=10,weighted=True, fluxWeighted=False, getRawCount=False)
            
        return im_dict
    
    def getLightCurve(self,fromLightCurveFile=True,fromObsFile=False,fromPhotonList=False,fromImageStack=False,save=False,**kwargs):
        assert 1.0*fromLightCurveFile+fromObsFile+fromPhotonList+fromImageStack==1, "Choose how to get data."
        
        if fromLightCurveFile:
            photometryFileName = self.path+os.sep+'FittedStacks'+os.sep+'FittedStack_'+self.fileID+'.h5'
            photometryDict = readPhotometryFile(photometryFileName)
        else:
            #Need to make a lightcurve
            #imStack_fn=self.path+os.sep+'ImageStacks'+os.sep+'ImageStack_1.h5'
            #target_fn=self.path+os.sep+'CentroidLists'+os.sep+'target'+os.sep+'Centroid_1.h5'
            #ref0_fn = self.path+os.sep+'CentroidLists'+os.sep+'ref0'+os.sep+'Centroid_1.h5'

            im_dict = self.getImages(fromObsFile=fromObsFile,fromPhotonList=fromPhotonList,fromImageStack=fromImageStack,**kwargs)
            images=im_dict['images']
            pixIntTimes=im_dict['pixIntTimes']
            centroids,flags = self.getCentroids()
            
            fluxDict_list = self.makeLightCurve(images=images,centroids=centroids,flags=flags,expTimes=pixIntTimes,save=save)
            
            if save:
                photometryFileName = self.path+os.sep+'FittedStacks'+os.sep+'FittedStack_'+self.fileID+'.h5'
                writePhotometryFile(fluxDict_list=fluxDict_list, im_dict = im_dict, filename = photometryFileName,verbose=self.kwargs['verbose'])

            startTimes = im_dict['startTimes']
            intTimes = im_dict['intTimes']
            flux = np.asarray([fluxDict['flux'] for fluxDict in fluxDict_list])
            parameters = np.asarray([fluxDict['parameters'] for fluxDict in fluxDict_list])
            perrors = np.asarray([fluxDict['mpperr'] for fluxDict in fluxDict_list])
            redChi2 = np.asarray([fluxDict['redChi2'] for fluxDict in fluxDict_list])
            flags = np.asarray([fluxDict['flag'] for fluxDict in fluxDict_list])
            photometryDict = {'startTimes': startTimes, 'intTimes': intTimes, 'flux': flux, 'parameters': parameters, 'perrors': perrors, 'redChi2': redChi2, 'flag': flags}
                
        return photometryDict

    def makeLightCurve(self,images,centroids,flags,expTimes=None,save=False):
        if expTimes==None:
            expTimes = [None]*len(images)
        #flux_list = []
        fluxDict_list = []
        for i in range(len(images)):
            fluxDict=self.performPhotometry(images[i],centroids[i],expTimes[i])
            if flags[i]>0:
                fluxDict['flux']=np.asarray(fluxDict['flux'])*0.
                fluxDict['flag']=flags[i]
            fluxDict_list.append(fluxDict)
            #flux_list.append(fluxDict['flux'])
        return fluxDict_list
                
    def performPhotometry(self,image,centroid,expTime=None):
        '''
        Perform the photometry on an image. 
        
        Input:
            image ----- 2D image of data (0 for dead pixel, shouldn't be any nan's or infs)
                        Should be fully calibrated, dead time corrected, and scaled up to the effective integration time
            centroid -- list of (col,row) tuples. The first tuple is the target location. The next are reference stars in the field
            expTime --- 2D array of pixel exposure times (0 for dead pixels)
                        optional. But can be used for distinguishing 0 count pixels from bad pixels
                        
        Output:
            fluxDict - dictionary output from photometry module. Should contain keyword 'flux'
                     - 'flux:' array of flux values. [target_flux, ref0_flux, ...]
        '''

        if self.PSF:
            PSFphoto = PSFphotometry(image,centroid,expTime,**self.kwargs)
            fluxDict = PSFphoto.PSFfit(aper_radius=5.)
            del PSFphoto
        elif self.aper:
            pass
        else:
            print "Choose PSF fitting or aperture photometry"
            
        return fluxDict
    
    
                
    def getFNs(self):
        '''
        Should loop through obs files. Return list of FileName objects
        '''

        obsFNs = []
        run = os.path.basename(os.path.dirname(self.path))
        for i in range(len(self.params['sunsetDates'])):
            for j in range(len(self.params['obsTimes'][i])):
                obsFN = FileName(run=run, date=self.params['sunsetDates'][i], tstamp=self.params['utcDates'][i]+'-'+self.params['obsTimes'][i][j])
                obsFNs.append(obsFN)
                
        return obsFNs


        
    def getCentroids(self,fromFile=True,**kwargs):
        '''
        Returns:
            centroids - a list of (col,row) tuples indicating the location of the stars in the field.
                        The first location should be the target star. The rest are reference stars.
            flags - flag from centroiding algorithm. 0 means good. 1 means failed
                    Same length as number of images. If one star fails the centroiding then whole image is flagged
        '''
        
        if fromFile:
            names = os.listdir(self.path+os.sep+'CentroidLists')
            ref_fns = []
            for name in names:
                if os.path.isdir(self.path+os.sep+'CentroidLists'+os.sep+name):
                    if name=='target':
                        target_fn = self.path+os.sep+'CentroidLists'+os.sep+name+os.sep+'Centroid_'+self.fileID+'.h5'
                    else:
                        ref_fns.append(name)
                        
            target_centroid, flags = self.getCentroidFromFile(target_fn)
            if len(ref_fns) > 0:
                centroid_list = [target_centroid]
                for ref_name in sorted(ref_fns):
                    ref_fn = self.path+os.sep+'CentroidLists'+os.sep+ref_name+os.sep+'Centroid_'+self.fileID+'.h5'
                    ref_centroid, ref_flags = self.getCentroidFromFile(ref_fn)
                    centroid_list.append(ref_centroid)
                    flags+=ref_flags
                
                centroids = np.asarray(zip(*centroid_list))
            else: centroids = target_centroid
            flags=1.0*(np.asarray(flags)>0.)
        else:
            #Open up image dict, use kwargs
            print 'Not implemented!'
            raise IOError
            
            
        return centroids, flags
        
    def getCentroidFromFile(self,star_fn):
        centroidFile = tables.openFile(star_fn, mode='r')
        xPos = np.array(centroidFile.root.centroidlist.xPositions.read())
        yPos = np.array(centroidFile.root.centroidlist.yPositions.read())
        flags = np.array(centroidFile.root.centroidlist.flags.read())
        centroids = zip(xPos,yPos)
        centroidFile.close()
        return centroids, flags
        


if __name__ == '__main__':
    path = '/Scratch/DisplayStack/PAL2014/HAT_P1'
    identifier = '1'
    verbose=True
    showPlot=False
    
    LC = LightCurve(path,fileID = identifier, PSF=True,verbose=verbose,showPlot=showPlot)
    #photometryFileName=path+os.sep+'FittedStacks'+os.sep+'photometry_1.h5'
    #imageStackFilename=path+os.sep+'ImageStacks'+os.sep+'ImageStack_1.h5'
    photometryDict = LC.getLightCurve(fromLightCurveFile=False,fromObsFile=False,fromPhotonList=False,fromImageStack=True,save=True)
    #photometryDict = LC.getLightCurve(fromLightCurveFile=True)
    
    def f(fig,axes):
        time = photometryDict['startTimes']
        flux=photometryDict['flux']
        tar_flux = flux[:,0]
        ref_flux = flux[:,1]
        flags = photometryDict['flag']
        
        axes.plot(time,tar_flux,'b.-',label='target')
        axes.plot(time,ref_flux,'g.-',label='ref')
        axes.plot(time[np.where(flags==1.)],ref_flux[np.where(flags==1.)],'ro',label='Failed Centroid')
        axes.plot(time[np.where(flags>1.1)],ref_flux[np.where(flags>1.1)],'ro',label='Failed Fit')
        axes.legend()
    pop(plotFunc=f,title='Flux')
    
    

    #LC = LightCurve(path,verbose=verbose,showPlot=showPlot)
    #obsFNs = LC.getFNs(fromObs = True, fromImageStack=False)
    ##print 'obs[10]: ',obsFNs[10:11][0].obs()
    #obsFiles = generateObsObjectList(obsFNs,wvlLowerLimit=3500, wvlUpperLimit=5000,loadHotPix=True,loadWvlCal=True,loadFlatCal=True,loadSpectralCal=False)
    #print 'numObs: ',len(obsFiles)
    #im_dict = getImages(fromObsFile=True,fromPhotonList=False,fromImageStack=False,obsFiles=obsFiles, integrationTime=10,weighted=True, fluxWeighted=False, getRawCount=False)
    #images=im_dict['images']
    #pixIntTimes=im_dict['pixIntTimes']
    #print 'numImages: ',len(images)

    #writeImageStack(images=images, pixIntTimes=pixIntTimes, startTimes=im_dict['startTimes'], intTimes=im_dict['intTimes'], path=path, outputFilename='ImageStack_1.h5')
    #im_dict = getImages(fromObsFile=False,fromPhotonList=False,fromImageStack=True,fullFilename=path+os.sep+'ImageStacks'+os.sep+'ImageStack_1.h5')
    #print im_dict['startTimes']
    
    #images = im_dict['images']
    #pixIntTimes=im_dict['pixIntTimes']
    #print 'numLoadedImages: ',len(images)
    #centroids = LC.getCentroids()
    #print 'numLoadedCentroids: ',len(centroids)
    #flux_list=LC.makeLightCurve(images,centroids,pixIntTimes)
    #print flux_list
    
    #pop(plotFunc=lambda fig,axes: axes.plot(flux_list),title="Flux")
    
    #for i in range(len(images)):
    #    im = images[i]
    #    centroid=LC.getCentroid()
    #    expTime = pixIntTimes[i]
    #    flux=LC.performPhotometry(im,centroid,expTime)
    #    print 'flux ',flux
 





