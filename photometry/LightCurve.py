'''
Author: Alex Walter
Date: Dec 3, 2014

This code makes a light curve using the photometry modules
'''

import os
import warnings
import numpy as np
import inspect

from util.FileName import FileName
from util.ObsFile import ObsFile
from util.readDict import readDict
from util.getImages import *
from util.popup import *
from photometry.PSFphotometry import PSFphotometry, writePsfPhotometryFile, readPsfPhotometryFile
from photometry.AperPhotometry import AperPhotometry, writeAperPhotometryFile, readAperPhotometryFile
from astrometry.CentroidCalc import quickCentroid,saveTable



class LightCurve():

    def __init__(self,fileID='test',path='/Scratch/DisplayStack/RUN_TEMPLATE/TARGET_TEMPLATE',targetName=None,run=None,verbose=False,showPlot=False):
        '''
        Constructs a list of obs FileName objects from the dictionary in the path.
        
        Inputs:
            fileID - identifier in filename. eg centroid_fileID.h5 or ImageStack_fileID.h5
            path - path to the display stack target info
            targetName - Name of target. Assumes there's a dictionary in path called targetName.dict
                         If None, just searches for any .dict file in the path
            run - name of run. eg PAL2014
                  if None, assumes run is the second to last directory in the path
            kwargs - key words for the photometry module

        '''
        
        self.path = path
        self.fileID=fileID
        self.targetName=targetName
        self.run=run
        self.verbose=verbose
        self.showPlot = showPlot

        if self.targetName is None or not os.path.isfile(self.path+os.sep+self.targetName+'.dict'):
            #Assumes there's only one '.dict' file in the directory
            for f in os.listdir(path):
                if f.endswith(".dict"):
                    if self.verbose: print 'Loading params from ',path+os.sep+f
                    self.targetName = f.split('.dict')[0]
        try:
            self.params = readDict(self.path+os.sep+self.targetName+'.dict')
            self.params.readFromFile(self.path+os.sep+self.targetName+'.dict')
        except:
            print "Provide a target name that leads to a dictionary in the path!"
            print path
            print "target: ",self.targetName
            raise ValueError
            
        if self.run is None:
            #run is assumed to be the second to last directory in the path
            self.run = os.path.basename(os.path.dirname(self.path))
                
        #Get images
        try:
            self.im_dict = self.getImages(fromImageStack=True)
        except:
            print "Need to make image stack for this object"

        #Get Centroids
        try:
            self.centroids, self.flags = self.getCentroids(fromCentroidFile=True)
        except:
            print "Need to make centroids for this object"
    
    
    
    def getLightCurve(self,photometryType,fromPhotometryFile=True,save=False,**kwargs):
        '''
        This code grabs the in lightcurve from the photometry file. (Makes a photometry file if it doesn't exist)
        
        Inputs:
            photometryType - Can be 'PSF' or 'aperture'
            fromPhotometryFile - True if you want to grab the lightcurve from the photometry file
                                 If False, it uses self.im_dict and self.centroids to generate a photometry file
            save - if True and not fromPhotometryFile then save the newly generated photometry file
            kwargs - For photometry module
        
        Returs:
            photometryDict - dictionary of results from photometry module
                                'flux' - list of flux [counts/sec]
                                'startTimes' - list of start times [julian date]
                                'flags' - list of flags
                                - Other diagnositcs (like fit parameters for PSF fitting)
            
        '''
        
        if fromPhotometryFile:
            if photometryType is 'PSF':
                photometryFileName = self.path+os.sep+'FittedStacks'+os.sep+'FittedStack_'+self.fileID+'.h5'
                photometryDict = readPsfPhotometryFile(photometryFileName)
            elif photometryType is 'aperture':
                photometryFileName = self.path+os.sep+'ApertureStacks'+os.sep+'ApertureStack_'+self.fileID+'.h5'
                photometryDict = readAperPhotometryFile(photometryFileName)
            else:
                print 'Choose a valid type of photometry to perform!'
                raise ValueError
        else:
            #Need to make a lightcurve
            images=self.im_dict['images']
            pixIntTimes=self.im_dict['pixIntTimes']
            centroids = self.centroids
            flags=self.flags
            fluxDict_list = self.makeLightCurve(photometryType = photometryType, images=images,centroids=centroids,flags=flags,expTimes=pixIntTimes,save=save,**kwargs)
            startTimes = self.im_dict['startTimes']
            intTimes = self.im_dict['intTimes']
            
            if photometryType is 'PSF':
                if save:
                    photometryFileName = self.path+os.sep+'FittedStacks'+os.sep+'FittedStack_'+self.fileID+'.h5'
                    writePsfPhotometryFile(fluxDict_list=fluxDict_list, im_dict = self.im_dict, filename = photometryFileName,verbose=self.verbose)

                flux = np.asarray([fluxDict['flux'] for fluxDict in fluxDict_list])
                parameters = np.asarray([fluxDict['parameters'] for fluxDict in fluxDict_list])
                perrors = np.asarray([fluxDict['mpperr'] for fluxDict in fluxDict_list])
                redChi2 = np.asarray([fluxDict['redChi2'] for fluxDict in fluxDict_list])
                flags = np.asarray([fluxDict['flag'] for fluxDict in fluxDict_list])
                photometryDict = {'startTimes': startTimes, 'intTimes': intTimes, 'flux': flux, 'parameters': parameters, 'perrors': perrors, 'redChi2': redChi2, 'flags': flags}
                
            elif photometryType is 'aperture':
                if save:
                    photometryFileName = self.path+os.sep+'ApertureStacks'+os.sep+'ApertureStack_'+self.fileID+'.h5'
                    writeAperPhotometryFile(fluxDict_list=fluxDict_list, im_dict = self.im_dict, filename = photometryFileName,verbose=self.verbose)
                    
                flux = np.asarray([fluxDict['flux'] for fluxDict in fluxDict_list])
                sky = np.asarray([fluxDict['sky'] for fluxDict in fluxDict_list])
                apertureRad = np.asarray([fluxDict['apertureRad'] for fluxDict in fluxDict_list])
                annulusInnerRad = np.asarray([fluxDict['annulusInnerRad'] for fluxDict in fluxDict_list])
                annulusOuterRad = np.asarray([fluxDict['annulusOuterRad'] for fluxDict in fluxDict_list])
                flags = np.asarray([fluxDict['flag'] for fluxDict in fluxDict_list])
                photometryDict = {'startTimes': startTimes, 'intTimes': intTimes, 'flux': flux, 'sky': sky, 'apertureRad': apertureRad, 'annulusInnerRad': annulusInnerRad,'annulusOuterRad': annulusOuterRad, 'flag': flags}
            else:
                print 'Choose a valid type of photometry to perform!'
                raise ValueError
                
        return photometryDict

    def makeLightCurve(self,photometryType, images,centroids,flags,expTimes=None,save=False,**kwargs):
        if expTimes==None:
            expTimes = [None]*len(images)
        #flux_list = []
        fluxDict_list = []
        for i in range(len(images)):
            fluxDict=self.performPhotometry(photometryType, images[i],centroids[i],expTimes[i],**kwargs)
            if flags[i]>0:
                fluxDict['flux']=np.asarray(fluxDict['flux'])*0.
                fluxDict['flag']=flags[i]
            fluxDict_list.append(fluxDict)
            #flux_list.append(fluxDict['flux'])
        return fluxDict_list

    def performPhotometry(self,photometryType,image,centroid,expTime=None,**kwargs):
        '''
        Perform the photometry on an image. 
        
        Input:
            image ----- 2D image of data (0 for dead pixel, shouldn't be any nan's or infs)
                        Should be fully calibrated, dead time corrected, and scaled up to the effective integration time
            centroid -- list of (col,row) tuples. The first tuple is the target location. The next are reference stars in the field
            expTime --- 2D array of pixel exposure times (0 for dead pixels)
                        optional. But can be used for distinguishing 0 count pixels from bad pixels
            kwargs ---- For photometry module
                        
        Output:
            fluxDict - dictionary output from photometry module. Should contain keyword 'flux'
                     - 'flux:' array of flux values. [target_flux, ref0_flux, ...]
        '''

        if photometryType is 'PSF':
            PSFphoto = PSFphotometry(image,centroid,expTime,verbose=self.verbose, showPlot=self.showPlot)
            fluxDict = PSFphoto.PSFfit(**kwargs)
            del PSFphoto
        elif photometryType is 'aperture':
            aperPhoto = AperPhotometry(image,centroid,expTime,verbose=self.verbose, showPlot=self.showPlot)
            fluxDict = aperPhoto.AperPhotometry(**kwargs)
            del aperPhoto
        else:
            print 'Choose a valid type of photometry to perform!'
            raise ValueError
            
        return fluxDict
    
    
    
                  
                  
    def getImages(self,fromImageStack=False,fromObsFile=False,fromPhotonList=False,**kwargs):
        '''
        Gets a list of images
        
        Inputs:
        if fromObsFile:
            You can pass keywords. These can include
             - keywords for generateObsObjectList() in getImages
             - keywords for getImages() in getImages
             - 'save' keyword (default is False)
             
         Returns an im_dict:
            images - [Total Photon Counts] list of calibrated images. 
            pixIntTimes - [Seconds] list of pixel exposure times for each image. Same shape as images
            startTimes - [Julian Date] list of startTimes for each image. 
            intTimes - [Seconds] list of image integration times for each image. 
        '''
        assert 1.0*fromObsFile+fromPhotonList+fromImageStack==1, "Choose whether to get images from a list of ObsFiles, from a list of PhotonLists, or from an .h5 file made by imagestack"
        
        if fromImageStack:
            imStack_fn=self.path+os.sep+'ImageStacks'+os.sep+'ImageStack_'+self.fileID+'.h5'
            im_dict=getImages(fromImageStack=fromImageStack,fullFilename = imStack_fn)
        elif fromObsFile:

            obsFNs = self.getFNs()
            #get keywords for generateObsObjectList()
            arg_names = set(inspect.getargspec(generateObsObjectList)[0])
            func_kwargs = {}
            for arg_name in arg_names:
                try:
                    func_kwargs[arg_name]=kwargs.pop(arg_name)
                except:
                    pass
            print 'Generating Obs objects, calibrating...'
            obsFiles = generateObsObjectList(obsFNs,**func_kwargs)
            save = kwargs.pop('save',False)
            print 'Creating images...'
            im_dict = getImages(fromObsFile=fromObsFile,obsFiles=obsFiles, **kwargs)
            if save:
                imStack_fn='ImageStack_'+self.fileID+'.h5'
                writeImageStack(im_dict['images'], im_dict['startTimes'], intTimes=im_dict['intTimes'], pixIntTimes=im_dict['pixIntTimes'], path=self.path,outputFilename=imStack_fn)
            
        else:
            print 'not implemented yet!'
            raise IOError
            
        self.im_dict = im_dict
        return im_dict
    
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
        
    def getCentroids(self,fromCentroidFile=True,**kwargs):
        '''
        This function returns a tuple of centroid lists and flags
        
        Inputs:
            fromCentroidFile - True if you want to grab centroids from the CentroidList directory of your object
            if not fromCentroidFile:
                You can pass keywords. These can include
                 - keywords for self.getImages()
                 - 'save' keyword (default is False)
        
        Returns:
            centroids - a list of (col,row) tuples indicating the location of the stars in the field.
                        The first location should be the target star. The rest are reference stars.
            flags - flag from centroiding algorithm. 0 means good. 1 means failed
                    Same length as number of images. If one star fails the centroiding then whole image is flagged
        '''
        centroidDir = 'CentroidLists'
        targetDir = 'Target'
        refDir = 'Reference'
        
        
        if fromCentroidFile:
            names = os.listdir(self.path+os.sep+centroidDir)
            ref_fns = []
            for name in names:
                if os.path.isdir(self.path+os.sep+centroidDir+os.sep+name):
                    if name==targetDir:
                        target_fn = self.path+os.sep+centroidDir+os.sep+name+os.sep+'Centroid_'+self.fileID+'.h5'
                    else:
                        ref_fns.append(name)
                        
            target_centroid, flags = self.getCentroidFromFile(target_fn)
            if len(ref_fns) > 0:
                centroid_list = [target_centroid]
                for ref_name in sorted(ref_fns):
                    ref_fn = self.path+os.sep+centroidDir+os.sep+ref_name+os.sep+'Centroid_'+self.fileID+'.h5'
                    ref_centroid, ref_flags = self.getCentroidFromFile(ref_fn)
                    centroid_list.append(ref_centroid)
                    flags+=ref_flags
                
                centroids = np.asarray(zip(*centroid_list))
            else: centroids = target_centroid
            flags=1.0*(np.asarray(flags)>0.)
            
        #Make centroids
        else:
            save = kwargs.pop('save',False)
            im_dict = self.getImages(**kwargs)
            print 'Choose target Star.'
            xPositionList,yPositionList,flagList=quickCentroid(im_dict['images'],radiusOfSearch=10,maxMove = 4,usePsfFit=False)
            if save:
                target_fn = self.path+os.sep+centroidDir+os.sep+targetDir+os.sep+'Centroid_'+self.fileID+'.h5'
                try:
                    os.mkdir(self.path+os.sep+centroidDir+os.sep+targetDir)
                except:
                    pass
                paramsList = [-1]*4
                timeList = [-1]*len(xPositionList)
                hourAngleList = timeList
                saveTable(target_fn,paramsList,timeList,xPositionList,yPositionList,hourAngleList,flagList)
            
            target_centroid = np.asarray(zip(xPositionList,yPositionList))
            flags = np.asarray(flagList)
            
            ref_num = 0
            #raw_input() for python 2.6. just use input() for python 3+
            chooseRef = raw_input('Choose Reference Star '+str(ref_num)+'? ')
            centroid_list = [target_centroid]
            while chooseRef is 'yes' or chooseRef is 'Yes' or chooseRef is 'y' or chooseRef is 'Y' or chooseRef is 'true' or chooseRef is 'True' or chooseRef is '1':
                xPositionList,yPositionList,flagList=quickCentroid(im_dict['images'],radiusOfSearch=10,maxMove = 4,usePsfFit=False)
                if save:
                    ref_num_str = "%02.d" % (ref_num,)  #Format as 2 digit integer with leading zeros
                    ref_fn = self.path+os.sep+centroidDir+os.sep+refDir+ref_num_str+os.sep+'Centroid_'+self.fileID+'.h5'
                    try:
                        os.mkdir(self.path+os.sep+centroidDir+os.sep+refDir+ref_num_str)
                    except:
                        pass
                    paramsList = [-1]*4
                    timeList = [-1]*len(xPositionList)
                    hourAngleList = timeList
                    saveTable(ref_fn,paramsList,timeList,xPositionList,yPositionList,hourAngleList,flagList)
                ref_centroid=np.asarray(zip(xPositionList,yPositionList))
                ref_flags=np.asarray(flagList)
                centroid_list.append(ref_centroid)
                flags+=ref_flags
                ref_num+=1
                chooseRef = raw_input('Choose Reference Star '+str(ref_num)+'? ')
                
            if len(centroid_list)>1: centroids = np.asarray(zip(*centroid_list))
            else: centroids = target_centroid
            flags=1.0*(np.asarray(flags)>0.)
        
        self.centroids = centroids
        self.flags=flags
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
    path = '/Scratch/DisplayStack/PAL2014/1SWASP_J2210'
    identifier = '0'
    
    LC=LightCurve(fileID=identifier,path=path,targetName=None,run=None,verbose=True,showPlot=True)

    #beammapFileName='/ScienceData/PAL2014/beammap_SCI6_B140731-Boba_20141118flip.h5'
    #im_dict = LC.getImages(fromObsFile=True,save=True,
    #             wvlLowerLimit=3000, wvlUpperLimit=8000, beammapFileName = beammapFileName, loadHotPix=False,loadWvlCal=True,loadFlatCal=True,loadSpectralCal=False,
    #             integrationTime=10, weighted=True, fluxWeighted=False, getRawCount=False, scaleByEffInt=False)
    #LC.getCentroids(fromCentroidFile=False,fromImageStack=True,save=True)
    
    photometryDict = LC.getLightCurve(photometryType='aperture',fromPhotometryFile=True,save=False)
    
    def f(self):
        time = photometryDict['startTimes']
        flux=photometryDict['flux']
        tar_flux = flux[:,0]
        ref_flux = flux[:,1]
        flags = photometryDict['flag']
        
        self.axes.plot(time,tar_flux,'b.-',label='target')
        self.axes.plot(time,ref_flux,'g.-',label='ref')
        self.axes.plot(time[np.where(flags==1.)],ref_flux[np.where(flags==1.)],'ro',label='Failed Centroid')
        self.axes.plot(time[np.where(flags>1.1)],ref_flux[np.where(flags>1.1)],'ro',label='Failed Fit')
        self.axes.legend()
    pop(plotFunc=f,title='Flux')
    
    def f(self):
        time = photometryDict['startTimes']
        flux=photometryDict['flux']
        tar_flux = flux[:,0]
        ref_flux = flux[:,1]
        flags = photometryDict['flag']
        
        self.axes.plot(time,ref_flux/tar_flux,'b.-',label='Ratio')
        #self.axes.plot(time,ref_flux,'g.-',label='ref')
        self.axes.plot(time[np.where(flags==1.)],ref_flux[np.where(flags==1.)],'ro',label='Failed Centroid')
        self.axes.plot(time[np.where(flags>1.1)],ref_flux[np.where(flags>1.1)],'ro',label='Failed Fit')
        self.axes.legend()
    pop(plotFunc=f,title='Flux Ratio')
    

