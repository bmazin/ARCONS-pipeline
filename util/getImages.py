
#This code is for getting calibrated images
import warnings
import os
import numpy as np
import tables
import hotpix.hotPixels as hp
from util.ObsFile import ObsFile
from util.FileName import FileName
from util.readDict import readDict
from headers.DisplayStackHeaders_old import *


def writeImageStack(images, startTimes, intTimes=None, pixIntTimes=None, path='/Scratch/DisplayStack/RUN_TEMPLATE/TARGET_TEMPLATE',outputFilename='ImageStack_0.h5'):
    '''
    Saves an h5 file containing the image stacks. See headers.DisplayStackHeaders for variable definitions
    
    Inputs:
        images - list of images. images[number, row, col]
        startTimes - list of startTimes [Julian Date]
        intTimes - (optional) list of integration times for each image [Seconds]
        pixIntTimes - (optional) list of pixel integration times. [number,row,col]  (same shape as images)
        path - path to target's data directory
        outputFilename - name of imagestack file
        
    '''
    # Create and h5 output file, create header and data groups
    fullFilename = path+os.sep+'ImageStacks'+os.sep+outputFilename
    if os.path.isfile(fullFilename):
        warnings.warn("Overwriting image stack file: "+str(fullFilename),UserWarning)
    fileh = tables.openFile(fullFilename, mode='w')
    #print fullFilename
    headerGroup = fileh.createGroup("/", headerGroupName, 'Header')
    stackGroup = fileh.createGroup("/", imagesGroupName, 'Image Stack')

    # Create row for header information
    headerTable = fileh.createTable(headerGroup, headerTableName, ImageStackHeaderDescription,'Header Info')
    header = headerTable.row
    
    run=os.path.basename(os.path.dirname(path))
    target=''
    params=None
    for f in os.listdir(path):
        if f.endswith(".dict"):
            target = f.split('.dict')[0]
            params = readDict(path+os.sep+f)
            params.readFromFile(path+os.sep+f)

    # Fill in the header with possibly useful information.
    header[runColName] = run
    header[targetColName] = target
    header[RAColName] = None if params==None else params['RA']
    header[DecColName] = None if params==None else params['Dec']
    header[HA_offsetColName] = None if params==None else params['HA_offset']
    
    nRows = len(images[0])
    nCols = len(images[0][0])
    header[nColColName] = nCols
    header[nRowColName] = nRows
    
    #header[lstColName] = self.lst
    #header[expTimeColName] = self.exptime
    #header[obsFileColName] = self.obsFn
    #header[integrationTimeColName] = self.integrationTime

    #if self.useDeadPixelMasking:
    #    header[deadPixColName] = self.deadPixelFilename
    #if self.useHotPixelMasking:
    #    header[hotPixColName] = 'Default'
    #if self.useWavelengthCalibration:
    #    header[wvlCalFileColName] = 'Default'
    #    header[lowWvlColName] = self.lowerWavelengthCutoff
    #    header[highWvlColName] = self.upperWavelengthCutoff
    #if self.useFlatCalibration:
    #    header[flatCalFileColName] = 'Default'
    header.append()

    # Create an h5 array for the starttime of each frame in the image cube.
    timeTable = fileh.createCArray(stackGroup, timeTableName, tables.Float64Atom(), (1,len(startTimes)))       
    timeTable[:] = startTimes
    # Create an h5 array for the integration time of each frame in the image cube.
    if intTimes!=None and len(intTimes)==len(startTimes):
        intTimeTable = fileh.createCArray(stackGroup, intTimeTableName, tables.Float64Atom(), (1,len(intTimes)))
        intTimeTable[:] = intTimes

    # Create an h5 table for the image cube.
    #print len(images)
    stackTable = fileh.createCArray(stackGroup, imagesTableName, tables.Float64Atom(), (nRows,nCols, len(np.asarray(images))))
    stackTable[:] = np.rollaxis(np.asarray(images),0,3)     #Paul has a weird [row,col,number] format...
    # Create an h5 table for the image cube.
    if pixIntTimes!=None and np.asarray(pixIntTimes).shape == np.asarray(images).shape:
        intTimeTable = fileh.createCArray(stackGroup, pixIntTimeTableName, tables.Float64Atom(), (nRows,nCols, len(np.asarray(pixIntTimes))))        
        intTimeTable[:] = np.rollaxis(np.asarray(pixIntTimes),0,3)     #Paul has a weird [row,col,number] format...

    # Flush the h5 output file
    fileh.flush()
    fileh.close()


def getImages(fromObsFile=False,fromPhotonList=False,fromImageStack=False,**kwargs):
    '''
    Inputs:
        Specify where you want to get your image from
        
    Return:
        dictionary with keywords:
        images - [Total Photon Counts] list of images. 
        pixIntTimes - [Seconds] list of pixel exposure times for each image. Same shape as images
        startTimes - [Julian Date] list of startTimes for each image. 
        intTimes - [Seconds] list of image integration times for each image. 
    '''
    assert 1.0*fromObsFile+fromPhotonList+fromImageStack==1, "Choose whether to get images from a list of ObsFiles, from a list of PhotonLists, or from an .h5 file made by imagestack"

    if fromObsFile:
        return getObsFileImages(**kwargs)
    elif fromPhotonList:
        print "Not implemented yet!"
        raise IOError
    elif fromImageStack:
        return loadImagesFromStack(**kwargs)
        
def loadImagesFromStack(fullFilename=''):
    print "Loading images from stack..."
    #fullFilename=imageStackFilename
    if not os.path.isfile(fullFilename): return None
    stackFile = tables.openFile(fullFilename, mode='r')
    images = stackFile.getNode('/',imagesGroupName)._f_getChild(imagesTableName).read()
    images = np.rollaxis(images,2,0)
    try:
        pixIntTimes = stackFile.getNode('/',imagesGroupName)._f_getChild(pixIntTimeTableName).read()
        pixIntTimes = np.rollaxis(pixIntTimes,2,0)
    except tables.exceptions.NoSuchNodeError:
        pixIntTimes = None
    startTimes = stackFile.getNode('/',imagesGroupName)._f_getChild(timeTableName).read()
    try:
        intTimes = stackFile.getNode('/',imagesGroupName)._f_getChild(intTimeTableName).read()
    except tables.exceptions.NoSuchNodeError:
        intTimes=None
        
    stackFile.close()
        
    return {'images':images,'pixIntTimes':pixIntTimes,'startTimes':np.asarray(startTimes).flatten(),'intTimes':np.asarray(intTimes).flatten()}
        
def generateObsObjectList(obsFNs,wvlLowerLimit=3000, wvlUpperLimit=13000,beammapFileName = None,loadHotPix=True,loadWvlCal=True,loadFlatCal=True,loadSpectralCal=True):

    obsFiles = []
    for obsFN in obsFNs:
        if type(obsFN) == type(''):     #Full path to obs file
            obs = ObsFile(obsFN)
        else:                           #FileName object
            obs = ObsFile(obsFN.obs())
            
        if beammapFileName is not None and os.path.isfile(beammapFileName):
            obs.loadBeammapFile(beammapFileName)
            
        obs.setWvlCutoffs(wvlLowerLimit=wvlLowerLimit, wvlUpperLimit=wvlUpperLimit)
        if loadHotPix:
            if not os.path.isfile(FileName(obsFile=obs).timeMask()):
                print "generating hp file ", FileName(obsFile=obs).timeMask()
                hp.findHotPixels(obsFile=obs,outputFileName=FileName(obsFile=obs).timeMask())
            obs.loadHotPixCalFile(FileName(obsFile=obs).timeMask(),reasons=['hot pixel','dead pixel'])
        if loadWvlCal:
            obs.loadBestWvlCalFile()
        if loadFlatCal:
            #obs.loadFlatCalFile(FileName(obsFile=obs).flatSoln())
            obs.loadFlatCalFile('/Scratch/flatCalSolnFiles/flatsol_1s.h5')
        if loadSpectralCal:
            pass
            
        obsFiles.append(obs)
        
    return obsFiles

        
def getObsFileImages(obsFiles, integrationTime=10,**kwargs):
    '''
    Inputs:
        obsFiles - Ordered list of ObsFile objects. All the calibrations should already be applied, ie flatcal loaded etc..
        integrationTime - Split obsfiles into chuncks this many seconds long. Can also be array of length obsFiles
        kwargs - To be passed to getObsFileImage()
        
    Return:
        dictionary with keywords:
        images - [Total Photon Counts] list of images. 
        pixIntTimes - [Seconds] list of pixel exposure times for each image. Same shape as images
        startTimes - [Julian Date] list of startTimes for each image. 
        intTimes - [Seconds] list of image integration times for each image. 
    '''
    try:
        len(integrationTime)
    except TypeError:
        integrationTime = [integrationTime]*len(obsFiles)
    
    images=[]
    pixIntTimes=[]
    startTimes=[]
    intTimes=[]
    for i in range(len(obsFiles)):
        obs=obsFiles[i]
        expTime = obs.getFromHeader('exptime')
        obsTime = obs.getFromHeader('jd')
        stepStarts = np.arange(0., expTime, integrationTime[i])  #Start time for each step (in seconds).
        stepEnds = stepStarts + integrationTime[i]               #End time for each step
        for startTime in stepStarts:
            imageStartTime = obsTime+1.0*startTime/(24.*60.*60.)
            image,pixIntTime = getObsFileImage(obs,startTime,integrationTime[i],**kwargs)
            images.append(image)
            pixIntTimes.append(pixIntTime)
            startTimes.append(imageStartTime)
            intTimes.append(integrationTime[i])
            
    return {'images':images,'pixIntTimes':pixIntTimes,'startTimes':startTimes,'intTimes':intTimes}
            
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









