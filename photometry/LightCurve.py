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
from util.popup import *
from util.ObsFileSeq import ObsFileSeq
from astrometry.CentroidCalc import quickCentroid
from photometry.PSFphotometry import PSFphotometry
from photometry.AperPhotometry import AperPhotometry
from headers.DisplayStackHeaders import writePhotometryFile, readPhotometryFile, PSFPhotometryDataDescription, aperPhotometryDataDescription, readImageStack, writeCentroidFile, readCentroidFile


def isPSFString(str_var):
    return str_var in ['point spread function', 'PSF', 'psf']
    
def isAperString(str_var):
    return str_var in ['aperture','Aperture','aper','Aper','aperture photometry','Aperture Photometry']


class LightCurve():
    def __init__(self,fileID='test',path='/Scratch/DisplayStack/RUN_TEMPLATE/TARGET_TEMPLATE',targetName=None,run=None,verbose=False,showPlot=False):
        '''
        Constructs a list of obs FileName objects from the dictionary in the path.
        
        Inputs:
            fileID - identifier in filename. eg centroid_fileID.h5 or ImageStack_fileID.h5
            path - path to the display stack target info
            targetName - Name of target. Assumes there's a dictionary in path called targetName.dict
                         If None, just searches for any xxx.dict file in the path and assigns targetName=xxx
            run - name of run. eg PAL2014
                  if None, assumes run is the second to last directory in the path
            verbose
            showPlot

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
            self.loadImageStack()
        except:
            print "Need to make image stack for this object"
        #Get Centroids
        try:
            self.loadAllCentroidFiles()
        except:
            print "Need to make centroids for this object"
            
    def viewLightCurve(self,photometryFilename='',photometryType=None):
        if photometryFilename is '':
            if isPSFString(photometryType): photometryFilename=self.path+os.sep+'FittedStacks'+os.sep+'FittedStack_'+self.fileID+'.h5'
            elif isAperString(photometryType): photometryFilename=self.path+os.sep+'ApertureStacks'+os.sep+'ApertureStack_'+self.fileID+'.h5'
            else:
                print 'You must provide a valid photometry type if you want to automatically generate the photometry filename'
                raise ValueError
        
        headerDict, dataDict = readPhotometryFile(photometryFilename)
        photometryType = headerDict['photometryType']
        #if isPSFString(photometryType):
        #    pop = PSF_Popup()
        print "implement this!"
        
    def loadLightCurve(self,photometryFilename='',photometryType=None):
        if photometryFilename is '':
            if isPSFString(photometryType): photometryFilename=self.path+os.sep+'FittedStacks'+os.sep+'FittedStack_'+self.fileID+'.h5'
            elif isAperString(photometryType): photometryFilename=self.path+os.sep+'ApertureStacks'+os.sep+'ApertureStack_'+self.fileID+'.h5'
            else:
                print 'You must provide a valid photometry type if you want to automatically generate the photometry filename'
                raise ValueError
        
        self.photometry_params, self.photometry_dict = readPhotometryFile(photometryFilename)
        self.photometryType = self.photometry_params['photometryType']
            
    def makeLightCurve(self,photometryType,photometryFilename='',**photometryKwargs):
        '''
        Loops through each image and performs photometry on it. 
        Saves the results in the file specified (photometryFilename='' saves in default location). 
        
        Inputs:
            photometryType - Can be 'PSF' or 'aperture.' See isPSFString() and isAperString()
            photometryFilename - full path to save photometry file. 
                                 If None then don't save.
                                 If empty string then save in self.path with self.fileID
            photometryKwargs - Provide parameters for the photometry operation.
                               You can give a single value for a parameter to be used on each image. ie., aper_radius=5
                               Or you can give a list of values (same length as number of images) to have a different value for that parameter for each image. ie., aper_radius=[5,4,...,5]
            
        Returns:
            photometryData - Dictionary where each key maps to a list of values corresponding to the value for each image
                           - valid keywords are the intersection of keywords from whatever the photometry module returns and the corresponding description in headers.DisplayStackHeaders
                                ie. for 'PSF', check PSFphotometry.PSFfit() in PSFphotometry.py and PSFPhotometryDataDescription() in DisplayStackHeaders.py
            
        Warning!
            If you need to pass a photometry parameter that is an array, and it happens that you have the 
            same number of images as the length of the array, then the code won't know if it is supposed 
            to pass the whole array with each image or one element of the array with each image. To avoid 
            this, pass a list of arrays for the images.
            ie., param=np.tile(array, (NFrames,1))
        '''
        if not hasattr(self,'im_dict'):
            print "Need to make image stack for this object"
            return
        if not hasattr(self,'centroids'):
            print "Need to make centroids for this object"
            return
        
        #Parameters needed for photometry class __init__()
        images=self.im_dict['images']
        assert len(images)>0
        expTimes=self.im_dict['pixIntTimes']
        #centroids = zip(self.centroid_dict['xPositions'], self.centroid_dict['yPositions'])
        #flags = self.centroid_dict['flag']
        flags=self.flags
        #Specific parameters for photometry operation
        if isPSFString(photometryType):
            photometryParamKeys=set(inspect.getargspec(PSFphotometry.PSFfit)[0])
            photometryDataKeys=PSFPhotometryDataDescription(1,1).keys()
        elif isAperString(photometryType):
            photometryParamKeys=set(inspect.getargspec(AperPhotometry.AperPhotometry)[0])
            photometryDataKeys=aperPhotometryDataDescription(1).keys()
        photometryParamKeys = list(photometryParamKeys.intersection(photometryKwargs.keys()))         #unique keywords given in **photometryKwargs that can be passed to photometry module
        
        #First add all the keywords that stay the same with each image
        photometryParamDict={}
        for keyword in np.copy(photometryParamKeys):
            val=photometryKwargs[keyword]
            try:
                assert type(val) is not str     #make sure it's not a single string
                val=list(val)                   #make sure it can be cast to a list
                if len(val) != len(images):
                    warnings.warn("Be careful passing an array as a parameter. You may want to pass a list of arrays, one for each image.",UserWarning)
                assert len(val)==len(images)    #make sure the list is the right length
            except:
                photometryParamDict[keyword]=val
                photometryParamKeys.remove(keyword)
                
        #Grab the rest of the keywords for the first image
        for keyword in photometryParamKeys:
            photometryParamDict[keyword]=photometryKwargs[keyword][0]
        #Perform photometry on the first image. 
        fluxDict=self.performPhotometry(photometryType, images[0],self.centroids[0],expTimes[0],**photometryParamDict)
        if flags[0]>0:
            #fluxDict['flux']=np.asarray(fluxDict['flux'])*0.       # Force flux --> if centroid flag is set
            fluxDict['flag']=flags[0]                               # Force photometry flag set if centroid flag was set
        #Check for unexpected dictionary keys
        if not set(fluxDict.keys()).issubset(photometryDataKeys):
            warnings.warn("The following keys returned by the photometry module dictionary "+
                          "don't match the keys expected by writePhotometryFile() and won't "+
                          "be saved: "+str(set(fluxDict.keys()).difference(photometryDataKeys)),UserWarning)

        #initialize dictionary of listed values for return data
        photometryData={}
        for keyword in set(fluxDict.keys()).intersection(photometryDataKeys):
            if len(images)>1:
                photometryData[keyword]=[fluxDict[keyword]]
            else:
                photometryData[keyword]=fluxDict[keyword]
        #Loop through the rest of the images
        for i in range(1,len(images)):
            for keyword in photometryParamKeys:
                photometryParamDict[keyword]=photometryKwargs[keyword][i]
            fluxDict=self.performPhotometry(photometryType, images[i],self.centroids[i],expTimes[i],**photometryParamDict)
            if flags[i]>0:
                #fluxDict['flux']=np.asarray(fluxDict['flux'])*0.       # Force flux --> if centroid flag is set
                fluxDict['flag']=flags[i]                               # Force photometry flag set if centroid flag was set
            for keyword in set(fluxDict.keys()).intersection(photometryDataKeys):
                photometryData[keyword].append(fluxDict[keyword])
        
        
        if photometryFilename is not None:
            if photometryFilename is '':
                if isPSFString(photometryType): photometryFilename=self.path+os.sep+'FittedStacks'+os.sep+'FittedStack_'+self.fileID+'.h5'
                elif isAperString(photometryType): photometryFilename=self.path+os.sep+'ApertureStacks'+os.sep+'ApertureStack_'+self.fileID+'.h5'
                else:
                    print 'Choose a valid photometry type!'
                    raise ValueError
            writePhotometryFile(photometryFilename, photometryType, targetName=self.targetName, run=self.run,
                                maxExposureTime=self.im_params['maxExposureTime'], imageStackFilename=self.imageStackFilename,
                                centroidFilenames=self.centroidFilenames, startTimes=self.im_dict['startTimes'],
                                endTimes=self.im_dict['endTimes'],intTimes=self.im_dict['intTimes'],**photometryData)

        if self.showPlot:
            self.viewLightCurve(photometryType,photometryFilename)
        return photometryData
        

    def performPhotometry(self,photometryType,image,centroid,expTime=None,**photometryKwargs):
        '''
        Perform the photometry on an image. 
        
        Input:
            photometryType - Can be 'PSF' or 'aperture.' See isPSFString() and isAperString()
            image - 2D image of data (0 for dead pixel, shouldn't be any nan's or infs)
                        Should be fully calibrated, dead time corrected, and scaled up to the effective integration time
            centroid - list of (col,row) tuples. The first tuple is the target location. The next are reference stars in the field
            expTime - 2D array of pixel exposure times (0 for dead pixels)
                        optional. But can be used for distinguishing 0 count pixels from bad pixels
            photometryKwargs - For photometry module
                        
        Output:
            fluxDict - dictionary output from photometry module. Should contain keyword 'flux'
                     - 'flux:' array of flux values. [target_flux, ref0_flux, ...]
        '''

        if isPSFString(photometryType):
            PSFphoto = PSFphotometry(image,centroid,expTime,verbose=self.verbose, showPlot=self.showPlot)
            fluxDict = PSFphoto.PSFfit(**photometryKwargs)
            del PSFphoto
        elif isAperString(photometryType):
            aperPhoto = AperPhotometry(image,centroid,expTime,verbose=self.verbose, showPlot=self.showPlot)
            fluxDict = aperPhoto.AperPhotometry(**photometryKwargs)
            del aperPhoto
        else:
            print 'Choose a valid type of photometry to perform!'
            raise ValueError
            
        return fluxDict

    def makeImageStack(self,imageStackFilename='',dt=30,wvlStart=None,wvlStop=None,
                           weighted=True, fluxWeighted=False, getRawCount=False, 
                           scaleByEffInt=True, deadTime=100.e-6):
        '''
        This function makes an image stack using the ObsFileSeq class
        
        Inputs:
            imageStackFilename - full path of file.
                               - An empty string means use the default location using self.path and self.fileID
            dt - the maximum number of seconds for one frame 
            keywords for image
        '''
        tsl = []
        for day_i in range(len(self.params['utcDates'])):
            for tstamp in self.params['obsTimes'][day_i]:
                tsl.append(self.params['utcDates'][day_i]+'-'+tstamp)
        #tsl=tsl[:2]        
        ofs = ObsFileSeq(name=self.targetName,run=self.run,date=self.params['sunsetDates'][0],timeStamps=tsl,dt=dt)
        if imageStackFilename is None or imageStackFilename is '': imageStackFilename = self.path+os.sep+'ImageStacks'+os.sep+'ImageStack_'+self.fileID+'.h5'
        self.im_params, self.im_dict = ofs.loadImageStack(imageStackFilename, wvlStart=wvlStart,wvlStop=wvlStop,
                                                           weighted=weighted, fluxWeighted=fluxWeighted, getRawCount=getRawCount, 
                                                           scaleByEffInt=scaleByEffInt, deadTime=deadTime)
        self.imageStackFilename = imageStackFilename
        #return self.im_params, self.im_dict       
                            
    def loadImageStack(self,imageStackFilename=''):
        '''
        This function will load in a new image stack from the file specified. 
        
        Inputs:
            imageStackFilename - full path of file.
                               - An empty string means use the default location using self.path and self.fileID
            kwargs - keywords for makeImageStack()
        '''
        if imageStackFilename is None or imageStackFilename is '': imageStackFilename = self.path+os.sep+'ImageStacks'+os.sep+'ImageStack_'+self.fileID+'.h5'
        self.im_params, self.im_dict = readImageStack(imageStackFilename)
        self.imageStackFilename = imageStackFilename
        #return self.im_params, self.im_dict

    def makeAllCentroidFiles(self,centroidFilenames=[''],radiusOfSearch=[10],maxMove=[4],usePsfFit=[False]):
        centroidDir = 'CentroidLists'
        targetDir = 'Target'
        refDir = 'Reference'
    
        try:
            num_images=len(self.im_dict['images'])
            assert num_images>0
            self.centroidFilenames=[]
        except:
            print "Need to make image stack for this object"
            return
        
        if centroidFilenames is None or len(centroidFilenames)<1 or centroidFilenames[0] is '':
            try: centroidFilenames[0] = self.path+os.sep+centroidDir+os.sep+targetDir+os.sep+'Centroid_'+self.fileID+'.h5'
            except NameError: centroidFilenames = [self.path+os.sep+centroidDir+os.sep+targetDir+os.sep+'Centroid_'+self.fileID+'.h5']
        for i in range(1,len(centroidFilenames)):
            ref_num_str = "%02.d" % (i-1,)  #Format as 2 digit integer with leading zeros
            centroidFilenames[i] = self.path+os.sep+centroidDir+os.sep+refDir+ref_num_str+os.sep+'Centroid_'+self.fileID+'.h5'
        
        centroid_list=[]
        flags=np.zeros(num_images)
        for i in range(len(centroidFilenames)):
            if i==0:
                print "\tSelect Target Star"
            else:
                print "\tSelect Reference Star #"+str(i-1)
            try: radiusOfSearch_i = radiusOfSearch[i]
            except IndexError: radiusOfSearch_i = radiusOfSearch[0]
            except TypeError: radiusOfSearch_i = radiusOfSearch
            try: maxMove_i = maxMove[i]
            except IndexError: maxMove_i = maxMove[0]
            except TypeError: maxMove_i = maxMove
            try: usePsfFit_i=usePsfFit[i]
            except IndexError: usePsfFit_i=usePsfFit[0]
            except TypeError: usePsfFit_i=usePsfFit
            self.makeCentroidFile(centroidFilenames[i],radiusOfSearch_i,maxMove_i,usePsfFit_i)
            centroid_list.append(self.centroids)
            flags+=self.flags
        
        if len(centroid_list)>1: self.centroids = np.asarray(zip(*centroid_list))
        self.flags=1.0*(np.asarray(flags)>0.)
        
    def loadAllCentroidFiles(self,centroidFilenames=[]):
        centroidDir = 'CentroidLists'
        targetDir = 'Target'
        refDir = 'Reference'
        try:
            num_images=len(self.im_dict['images'])
            assert num_images>0
            self.centroidFilenames=[]
        except:
            print "Need to make image stack before loading centroids"
            return
            
        if centroidFilenames is None or len(centroidFilenames)<1:
            nStars=0
            for file_i in os.listdir(path+os.sep+centroidDir):
                nStars+=int(os.path.isdir(path+os.sep+centroidDir+os.sep+file_i))
            centroidFilenames=['']*nStars
        if centroidFilenames[0] is '':
            try: centroidFilenames[0] = self.path+os.sep+centroidDir+os.sep+targetDir+os.sep+'Centroid_'+self.fileID+'.h5'
            except NameError: centroidFilenames = [self.path+os.sep+centroidDir+os.sep+targetDir+os.sep+'Centroid_'+self.fileID+'.h5']
        for i in range(1,len(centroidFilenames)):
            ref_num_str = "%02.d" % (i-1,)  #Format as 2 digit integer with leading zeros
            centroidFilenames[i] = self.path+os.sep+centroidDir+os.sep+refDir+ref_num_str+os.sep+'Centroid_'+self.fileID+'.h5'
        
        centroid_list=[]
        flags=np.zeros(num_images)
        for i in range(len(centroidFilenames)):
            self.loadCentroidFile(centroidFilenames[i])
            centroid_list.append(self.centroids)
            flags+=self.flags
        
        if len(centroid_list)>1: self.centroids = np.asarray(zip(*centroid_list))
        self.flags=1.0*(np.asarray(flags)>0.)
        self.centroidFilenames = centroidFilenames

    def makeCentroidFile(self,centroidFilename='',radiusOfSearch=10,maxMove=4,usePsfFit=False):
        '''
        This function makes a centroid file using CentroidCal.quickCentroid()
        
        Inputs:
            centroidFilename - full path of file.
                             - An empty string means use the default location using self.path and self.fileID
            kwargs - 
            
        Returns:
            centroids
            flags
        '''
        centroidDir = 'CentroidLists'
        targetDir = 'Target'
        refDir = 'Reference'
        #Get images
        try:
            images=self.im_dict['images']
        except:
            print "Need to make image stack for this object"
            return
        
        xPositionList,yPositionList,flagList=quickCentroid(images,radiusOfSearch=radiusOfSearch,maxMove = maxMove,usePsfFit=usePsfFit)
        
        if centroidFilename is None or centroidFilename is '':
            centroidFilename = self.path+os.sep+centroidDir+os.sep+targetDir+os.sep+'Centroid_'+self.fileID+'.h5'
        
        centroid_params = {'targetName':self.targetName, 'run':self.run, 'nFrames':len(images),
                                'imageStackFilename':self.imageStackFilename}
        centroid_dict = {'startTimes':self.im_dict['startTimes'], 'endTimes':self.im_dict['endTimes'],
                              'intTimes':self.im_dict['intTimes'], 'xPositions':xPositionList,
                              'yPositions':yPositionList, 'flag':flagList}
        #writeCentroidFile(centroidFilename, **centroid_params,**centroid_dict)
        writeCentroidFile(centroidFilename,**dict(centroid_params.items() + centroid_dict.items()))
        try: self.centroidFilenames.append(centroidFilename)
        except NameError: self.centroidFilenames=[centroidFilename]
        
        self.centroids = zip(xPositionList,yPositionList)
        self.flags=flagList
        
        
    def loadCentroidFile(self,centroidFilename=''):
        centroidDir = 'CentroidLists'
        targetDir = 'Target'
        refDir = 'Reference'
        if centroidFilename is None or centroidFilename is '':
            centroidFilename = self.path+os.sep+centroidDir+os.sep+targetDir+os.sep+'Centroid_'+self.fileID+'.h5'
        centroid_params, centroid_dict = readCentroidFile(centroidFilename)
        self.centroids = zip(centroid_dict['xPositions'],centroid_dict['yPositions'])
        self.flags=centroid_dict['flag']
        try: self.centroidFilenames.append(centroidFilename)
        except NameError: self.centroidFilenames=[centroidFilename]
    

                  
if __name__ == '__main__':
    path = '/Scratch/DisplayStack/PAL2014/1SWASP_J2210'
    identifier = 'manHotPix'
    #path = '/Scratch/DisplayStack/PAL2014/HAT_P1'
    #identifier = '4000-5000_flat'
    
    
    LC=LightCurve(fileID=identifier,path=path,targetName=None,run=None,verbose=True,showPlot=False)
    #LC.makeImageStack(imageStackFilename='',dt=30,wvlStart=4000,wvlStop=9000,
    #                       weighted=True, fluxWeighted=False, getRawCount=False, 
    #                       scaleByEffInt=True, deadTime=100.e-6)
    #LC.makeAllCentroidFiles(centroidFilenames=['',''])
    #print LC.centroids
    #print LC.flags
    #LC.makeLightCurve(photometryType='PSF')
    LC.loadLightCurve(photometryType='aper')
    photometryDict=LC.photometry_dict
    
    
    #transit_mid = 2456953.81673
    #duration = 0.1166
    def f(self):
        time = photometryDict['startTimes']
        flux=photometryDict['flux']
        tar_flux = flux[:,0]
        ref_flux = flux[:,1]
        flags = photometryDict['flag']
        
        self.axes.plot(time,tar_flux,'b.',label='target')
        self.axes.plot(time,ref_flux,'g.',label='ref')
        #self.axes.axvline(x=transit_mid,c='k',ls='-')
        #self.axes.axvline(x=transit_mid-duration/2.,c='k',ls='--')
        #self.axes.axvline(x=transit_mid+duration/2.,c='k',ls='--')
        #self.axes.plot(time[np.where(flags==1.)],ref_flux[np.where(flags==1.)],'ro',label='Failed Centroid')
        #self.axes.plot(time[np.where(flags>1.1)],ref_flux[np.where(flags>1.1)],'ro',label='Failed Fit')
        self.axes.legend()
    pop(plotFunc=f,title='Flux')

