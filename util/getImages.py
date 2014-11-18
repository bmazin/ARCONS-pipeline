
#This code is for getting calibrated images
import warnings
import numpy as np

from util.ObsFile import ObsFile
from util.FileName import FileName




def getImages(fromObsFile=True,fromPhotonList=False,fromImageStack=False,**kwargs):
    '''
    Inputs:
        Specify where you want to get your image from
        
    Return:
        dictionary with keywords:
        images - list of images
        expTimes - list of pixel exposure times for each image
    '''
    assert 1.0*fromObsFile+fromPhotonList+fromImageStack==1, "Choose whether to get images from a list of ObsFiles, from a list of PhotonLists, or from an .h5 file made by imagestack"
    
    if fromObsFile:
        return getObsFileImages(**kwargs)
    elif fromPhotonList:
        pass
    elif fromImageStack:
        pass
        
def generateObsObjectList(obsFNs,wvlLowerLimit=3000, wvlUpperLimit=13000,loadHotPix=True,loadWvlCal=True,loadFlatCal=True,loadSpectralCal=True):

    obsFiles = []
    for obsFN in obsFNs:
        if type(obsFN) == type(''):     #Full path to obs file
            obs = ObsFile(obsFN)
        else:                           #FileName object
            obs = ObsFile(obsFN.obs())
            
        obs.setWvlCutoffs(wvlLowerLimit=wvlLowerLimit, wvlUpperLimit=wvlUpperLimit)
        if loadHotPix:
            obs.loadHotPixCalFile(FileName(obsFile=obs).timeMask(),reasons=['hot pixel','dead pixel'])
        if loadWvlCal:
            obs.loadBestWvlCalFile()
        if loadFlatCal:
            obs.loadFlatCalFile(FileName(obsFile=obs).flatSoln())
        if loadSpectralCal:
            pass
            
        obsFiles.append(obs)
        
    return obsFiles

        
def getObsFileImages(obsFiles, integrationTime=10,**kwargs):
    '''
    Inputs:
        obsFiles - Ordered list of ObsFile objects. All the calibrations should already be applied, ie flatcal loaded etc..
        integrationTime - Split obsfiles into chuncks this many seconds long
        
    Return:
        dictionary with keywords:
        images - list of images
        expTimes - list of pixel exposure times for each image
    '''
    images=[]
    expTimes=[]
    for obs in obsFiles:
        intTime = obs.getFromHeader('exptime')
        stepStarts = np.arange(0., intTime, integrationTime)  #Start time for each step (in seconds).
        stepEnds = stepStarts + integrationTime               #End time for each step
        for startTime in stepStarts:
            image,expTime = getObsFileImage(obs,startTime,integrationTime,**kwargs)
            images.append(image)
            expTimes.append(expTime)
            
    return {'images':images,'expTimes':expTimes}
            
def getObsFileImage(obs,startTime,integrationTime,deadTime=100.e-6,**kwargs):
    im_dict = obs.getPixelCountImage(firstSec=startTime, integrationTime=integrationTime,**kwargs)
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
    
    return im, im_dict['effIntTimes']









