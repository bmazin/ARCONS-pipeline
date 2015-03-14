'''
Author: Alex Walter
Date: Feb 2, 2015

This file contains the keywords and descriptions for display stack files
including image stacks, centroid files, and photometry files. It also
has the read and write functions.

Description functions:
    imageStackHeaderDescription()
    imageStackDataDescription()
    centroidHeaderDescription()
    centroidDataDescription()
    photometryHeaderDescription()
    PSFPhotometryDataDescription()
    aperPhotometryDataDescription()

Major functions:
    writeImageStack()
    readImageStack()
    writeCentroidFile()
    readCentroidFile()
    writePhotometryFile()
    readPhotometryFile()
    
Supporting functions:
    writeH5Table()
    readH5Table()
    
'''


import os, warnings
import tables
import numpy as np

# Parameter name strings (otherwise as defined in imageStackHeaderDescription, imageStackDataDescription)
imageListStr = 'images'
timeStampsStr = 'tStamps'
numFramesStr = 'nFrames'
photometryTypeStr = 'photometryType'
fluxesStr = 'flux'
paramsStr = 'parameters'
xPositionsStr = 'xPositions'
yPositionsStr = 'yPositions'

# Data group and table names
imageDataGroupName = 'ImageData'
imageDataTableName = 'ImageData'
imageHeaderGroupName = 'Header'
imageHeaderTableName = 'Header'
photometryDataGroupName = 'PhotometryData'
photometryDataTableName = 'PhotometryData'
photometryHeaderGroupName = 'Header'
photometryHeaderTableName = 'Header'
centroidDataGroupName = 'CentroidData'
centroidDataTableName = 'CentroidData'
centroidHeaderGroupName = 'Header'
centroidHeaderTableName = 'Header'


def centroidHeaderDescription():
    '''
    Returns a dictionary for the header table for a centroid file
    '''
    description = {
        'targetName'        : tables.StringCol(100, dflt=''),               # name of target observed
        'run'               : tables.StringCol(100, dflt=''),               # telescope run. ie., 'PAL2014'
        numFramesStr        : tables.UInt32Col(dflt=0),                     # Number of frames in image stack   [int]
        'imageStackFilename': tables.StringCol(200,dflt='')}                # full file name to image stack used for making centroid file
    #probably should add star RA and DEC
    return description

def centroidDataDescription():
    '''
    Returns a dictionary for the data table description of a centroid file
    '''
    description={
        'startTimes'    : tables.Float64Col(),      # Column for start times of image [Julian Date]
        'endTimes'      : tables.Float64Col(),      # Column for end times of image [Julian Date]
        'intTimes'      : tables.Float64Col(),      # Column for total exposure of image [seconds]
        xPositionsStr   : tables.Float64Col(),      # horizontal location of centroid on array [pixels]
        yPositionsStr   : tables.Float64Col(),      # vertical location of centroid on array [pixels]
        'flag'          : tables.UInt16Col()}       # flag to indicate if centroid is good (0)
    #Probably should add hourangle list
    return description


def imageStackHeaderDescription(nObs):
    '''
    Returns a dictionary for the header table description
    
    nObs - Number of obsFiles used to make image stack. Same as length of tStamps list
    '''
    description = {
        'targetName'        : tables.StringCol(100, dflt=''),               # name of target observed
        'run'               : tables.StringCol(100, dflt=''),               # telescope run. ie., 'PAL2014'
        numFramesStr        : tables.UInt32Col(dflt=0),                     # Number of frames in image stack   [int]
        'wvlLowerLimit'     : tables.Float64Col(dflt=np.nan),               # short wavelength limit of photons used in images  [Angstroms]
        'wvlUpperLimit'     : tables.Float64Col(dflt=np.nan),               # long wavelength limit     [Angstroms]
        'weighted'          : tables.UInt8Col(dflt=False),                  # True if flat cal applied
        'fluxWeighted'      : tables.UInt8Col(dflt=False),                  # True if flux cal applied
        'hotPixelsApplied'  : tables.UInt8Col(dflt=False),                  # True if hotpix was applied
        'maxExposureTime'   : tables.Float64Col(dflt=0.),                   # dt in ObsFileSeq.__init__()   [seconds]
        timeStampsStr       : tables.StringCol(15, dflt='', shape=nObs)}    # List of time stamps of obsfiles used. 15 chars each, ie., '20140923-111213'
        #Could add object's RA and Dec. HA offset. Etc..
    return description
    
def imageStackDataDescription(nRows,nCols):
    '''
    Returns a dictionary for the data table description
    
    nRows and nCols describe the shape of the image in pixels
    '''
    description = {
        imageListStr    : tables.Float64Col(dflt=np.nan, shape=(nRows,nCols)),  # Column for 2D image [Total counts during exposure (scaled up to intTime)]
        'pixIntTimes'   : tables.Float64Col(dflt=0., shape=(nRows,nCols)),      # Column for 2D pixel integration times [seconds]
        'startTimes'    : tables.Float64Col(),                                  # Column for start times [Julian Date]
        'endTimes'      : tables.Float64Col(),                                  # Column for end times [Julian Date]
        'intTimes'      : tables.Float64Col()}                                  # Column for total exposure [seconds]
        #Could probably add RA/Dec info here too
    return description

def photometryHeaderDescription(numStars):
    '''
    Returns a dictionary for the header table description for photometry files
    '''
    description = {
        'targetName'        : tables.StringCol(100, dflt=''),                   # name of target observed
        'run'               : tables.StringCol(100, dflt=''),                   # telescope run. ie., 'PAL2014'
        numFramesStr        : tables.UInt32Col(dflt=0),                         # Number of frames in image stack   [int]
        'maxExposureTime'   : tables.Float64Col(dflt=0.),                       # dt in ObsFileSeq.__init__()   [seconds]
        photometryTypeStr   : tables.StringCol(100,dflt=''),                    # Type of photometry used. ie., 'PSF'
        'imageStackFilename': tables.StringCol(200,dflt=''),                    # full file name to image stack used for making photometry file
        'centroidFilenames' : tables.StringCol(200,dflt='',shape=numStars)}     # full file name to centroid file used for making photometry file
    return description

def PSFPhotometryDataDescription(numStars,numParams):
    '''
    Returns a dictionary for the data table for PSF photometry files
    
    numStars - Number of stars in the images
    numParams - number of parameters in PSF fitting
    '''
    description = {
        fluxesStr       : tables.Float64Col(numStars),          # array of flux values. [target_flux, ref0_flux, ...]
        'startTimes'    : tables.Float64Col(),                  # Start time of images
        'endTimes'      : tables.Float64Col(),
        'intTimes'      : tables.Float64Col(),                  # integration time of individual images
        'model'         : tables.StringCol(100,dflt=''),        # Type of model used to fit image
        'flag'          : tables.UInt16Col(),                   # flag to indicate if fit is good (0)  
        paramsStr       : tables.Float64Col(numParams),         # parameters used to fit data
        'perrors'       : tables.Float64Col(numParams),         # the errors on the fits
        'redChi2'       : tables.Float64Col()}                  # reduced chi^2 of the fits
    return description
    
def aperPhotometryDataDescription(numStars):
    '''
    Returns a dictionary for the data table for aperture photometry files
    
    numStars - Number of stars in the images
    '''
    description = {
        fluxesStr           : tables.Float64Col(numStars),      # array of flux values. [target_flux, ref0_flux, ...]
        'startTimes'        : tables.Float64Col(),              # Start time of images
        'endTimes'          : tables.Float64Col(),
        'intTimes'          : tables.Float64Col(),              # integration time of individual images
        'flag'              : tables.UInt16Col(),               # flag to indicate if fit is good (0)  
        'skyFlux'           : tables.Float64Col(numStars),     # array of sky flux values, scaled to same n_pix as object flux. [target_sky, ref0_sky, ...]
        'apertureRad'       : tables.Float64Col(numStars),     # radius of aperture used for each object [target_apertureRad, ref0_apertureRad, ...]
        'annulusInnerRad'   : tables.Float64Col(numStars),     # radius of inner annulus used for sky subtraction. 0s if sky fitting was used.
        'annulusOuterRad'   : tables.Float64Col(numStars),     # radius of outer annulus used for sky subtraction. 0s if sky fitting was used
        'interpolation'     : tables.StringCol(100,dflt='')}    # type of interpolation performed on image before doing photometry

    return description
    
def writeCentroidFile(fullFilename, **dataKwargs):
    '''
    Writes data to centroid file
    
    Input:
        fullFilename - full path to file
        **dataKwargs - any keyword from centroidHeaderDescription() and centroidDataDescription()
    '''
    assert len(dataKwargs[xPositionsStr]) > 0                                   #Make sure there are xPosition centroids
    assert len(dataKwargs[xPositionsStr]) == len(dataKwargs[yPositionsStr])     #Make sure there are the same number of yPositions
    nFrames = len(dataKwargs[xPositionsStr])
    dataKwargs[numFramesStr] = nFrames
    
    # Create an h5 output file, create header and data groups
    if os.path.isfile(fullFilename):
        warnings.warn("Overwriting file: "+str(fullFilename),UserWarning)
    fileh = tables.openFile(fullFilename, mode='w')
    
    # Create table for header information
    headerTableDescription = centroidHeaderDescription()
    writeH5Table(fileh,centroidHeaderGroupName,centroidHeaderTableName,headerTableDescription,nRows=1,tableTitle='Header',**dataKwargs)
    
    # Create a table for the image data
    dataTableDescription = centroidDataDescription()
    filt = tables.Filters(complevel=2, complib='zlib') #See doc for pyTables 2.3 http://pytables.readthedocs.org/en/v.2.3.1/usersguide/libref.html#filtersclassdescr
    writeH5Table(fileh,centroidDataGroupName,centroidDataTableName,dataTableDescription,nRows=nFrames,tableTitle='Centroid Data',tableFilter=filt,**dataKwargs)

    # Flush the h5 output file
    fileh.flush()
    fileh.close()

def readCentroidFile(fileName):
    '''
    Returns dictionary of data in centroid file
    Input:
        fileName - full path to file
    Returns two dictionaries
        headerDict - dictionary containnig header info
        dataDict - dictionary containing data
    '''
    fileh = tables.openFile(fileName, mode='r')
    headerDict=readH5Table(fileh,centroidHeaderGroupName,centroidHeaderTableName)
    dataDict=readH5Table(fileh,centroidDataGroupName,centroidDataTableName)
    fileh.close()
    return headerDict, dataDict
    
def writePhotometryFile(fullFilename, photometryType, **dataKwargs):
    '''
    Writes data to photometry file
    
    Input:
        fullFilename - full path to file
        photometryType - Type of photometry performed. ie., 'PSF'
        **dataKwargs - any keyword from photometryHeaderDescription() and the corresponding PSFPhotometryDataDescription() or aperPhotometryDataDescription()
    '''
    from photometry.LightCurve import isPSFString, isAperString     #To avoid circular imports
    assert isPSFString(photometryType) or isAperString(photometryType)      #Make sure a valid photometry type was performed
    assert len(dataKwargs[fluxesStr]) > 0       #Make sure there was at least one image analyzed
    assert len(dataKwargs[fluxesStr][0]) > 0    #Make sure there's at least one star in the images
    dataKwargs[photometryTypeStr] = photometryType
    nFrames = len(dataKwargs[fluxesStr])
    dataKwargs[numFramesStr] = nFrames
    numStars=len(dataKwargs[fluxesStr][0])

    # Create an h5 output file, create header and data groups
    if os.path.isfile(fullFilename):
        warnings.warn("Overwriting file: "+str(fullFilename),UserWarning)
    fileh = tables.openFile(fullFilename, mode='w')
    
    # Create table for header information
    headerTableDescription = photometryHeaderDescription(numStars)
    writeH5Table(fileh,photometryHeaderGroupName,photometryHeaderTableName,headerTableDescription,nRows=1,tableTitle='Header',**dataKwargs)
        
    # Create a table for the photometry data
    if isPSFString(photometryType):
        try: numParams=len(dataKwargs[paramsStr][0])
        except: numParams=1
        dataTableDescription = PSFPhotometryDataDescription(numStars,numParams)
    elif isAperString(photometryType):
        dataTableDescription = aperPhotometryDataDescription(numStars)
    filt = tables.Filters(complevel=2, complib='zlib') #See doc for pyTables 2.3 http://pytables.readthedocs.org/en/v.2.3.1/usersguide/libref.html#filtersclassdescr
    writeH5Table(fileh,photometryDataGroupName,photometryDataTableName,dataTableDescription,nRows=nFrames,tableTitle='Photometry Data',tableFilter=filt,**dataKwargs)

    # Flush the h5 output file
    fileh.flush()
    fileh.close()
    
def readPhotometryFile(fileName):
    '''
    Returns dictionary of data in photometryFile
    Input:
        fileName - full path to file
    Returns two dictionaries
        headerDict - dictionary containnig header info
        dataDict - dictionary containing data
    '''
    fileh = tables.openFile(fileName, mode='r')
    headerDict=readH5Table(fileh,photometryHeaderGroupName,photometryHeaderTableName)
    dataDict=readH5Table(fileh,photometryDataGroupName,photometryDataTableName)
    fileh.close()
    return headerDict, dataDict

    
def writeImageStack(fullFilename, images, **dataKwargs):
                    #endTimes, intTimes, pixIntTimes, 
                    #targetName,run,nFrames,wvlLowerLimit,wvlUpperLimit,weighted,fluxWeighted,hotPixelsApplied,maxExposureTime,obsTStamps):
    '''
    Saves an h5 file containing the image stacks. See headers.DisplayStackHeaders for variable definitions
    
    Inputs:
        fullFilename - name of imagestack file
        images - list of images. images[number, row, col]
        **dataKwargs --> All the keywords from imageStackHeaderDescription and imageStackDataDescription can be used
            startTimes - list of startTimes [Julian Date]
            endTimes - list of end Times [Julian Date]. Not the same as intTime since the image may span over multiple obsFiles
            intTimes - list of integration times for each image [Seconds]
            pixIntTimes - array of pixel integration times. [number,row,col]  (same shape as images)
            targetName - 
            run - telescope run. ie., 'PAL2014'
            nFrames - You can include this keyword but it will be overwritten with the number of images in the list of images
            wvlLowerLimit - 
            wvlUpperLimit - 
            weighted - True if flat cal applied
            fluxWeighted - True if flux cal applied
            hotPixelsApplied - True if hotpix was applied
            maxExposureTime - dt in ObsFileSeq.__init__()
            tStamps - list of tstamps of included obs files
    '''
    assert(len(images)>0)  #Make sure some images are given
    dataKwargs[imageListStr]=images
    nFrames = len(images)
    nObs = len(dataKwargs.get(timeStampsStr,[0]))
    (nRows,nCols) = np.asarray(images[0]).shape
    dataKwargs[numFramesStr] = nFrames
    
    # Create an h5 output file, create header and data groups
    if os.path.isfile(fullFilename):
        warnings.warn("Overwriting image stack file: "+str(fullFilename),UserWarning)
    fileh = tables.openFile(fullFilename, mode='w')
   
    # Create table for header information
    headerTableDescription = imageStackHeaderDescription(nObs)
    writeH5Table(fileh,imageHeaderGroupName,imageHeaderTableName,headerTableDescription,nRows=1,tableTitle='Header',**dataKwargs)
        
    # Create a table for the image data
    dataTableDescription = imageStackDataDescription(nRows,nCols)
    filt = tables.Filters(complevel=2, complib='zlib') #See doc for pyTables 2.3 http://pytables.readthedocs.org/en/v.2.3.1/usersguide/libref.html#filtersclassdescr
    writeH5Table(fileh,imageDataGroupName,imageDataTableName,dataTableDescription,nRows=nFrames,tableTitle='Image Data',tableFilter=filt,**dataKwargs)

    # Flush the h5 output file
    fileh.flush()
    fileh.close()
    
    
def readImageStack(fileName):
    '''
    Returns dictionary of data in image stack file
    Input:
        fileName - full path to file
    Returns dict with keywords from imageStackHeaderDescription, imageStackDataDescription
        images
        pixIntTimes
        startTimes
        endTimes
        intTimes
        headerInfo...

    '''
    fileh = tables.openFile(fileName, mode='r')
    headerDict=readH5Table(fileh,imageHeaderGroupName,imageHeaderTableName)
    dataDict=readH5Table(fileh,imageDataGroupName,imageDataTableName)
    #dataDict['headerInfo']=headerDict
    fileh.close()
    return headerDict, dataDict
    

def writeH5Table(fileh,tableGroupName,tableName,tableDescription,nRows,tableTitle='',tableFilter=None,**dataKwargs):
    '''
    This function makes adding a table to an h5 file easy
    
    Inputs:
        fileh - tables.File instance
        tableGroupName - name of group node to put table in
        tableName - name of table node
        tableDescription - Dictionary of column descriptions for table
        nRows - number of rows for the table
        tableTitle - optional title for table
        tableFilter - tables.Filters instance. See doc for pyTables 2.3 http://pytables.readthedocs.org/en/v.2.3.1/usersguide/libref.html#filtersclassdescr
        **dataKwargs - keyword parameters for data with keyword matching that in tableDescription
    '''
    tableGroup=fileh.createGroup("/",tableGroupName)
    table=fileh.createTable(where=tableGroup,name=tableName,description=tableDescription,title=tableTitle,filters=tableFilter,expectedrows=nRows)
    validDataKeys = set(dataKwargs.keys()).intersection(tableDescription.keys())
    if nRows>1:
        for i in range(nRows):
            row = table.row
            for keyword in validDataKeys:
                row[keyword] = dataKwargs[keyword][i]
            row.append()
    else:
        row = table.row
        for keyword in validDataKeys:
            row[keyword] = dataKwargs[keyword]
        row.append()
    table.flush()
    table.close()


def readH5Table(fileh, groupName,tableName):
    '''
    This function reads the data a table in the specified group node and returns the info as a dictionary
    
    Inputs:
        fileh - tables.File instance
        groupName - name of group node
        tableName - name of table
    Returns:
        dictionary of table data. Keywords are the table column descriptions.
    '''
    data=fileh.getNode('/' + groupName, tableName).read()
    returnDict = {}
    try:
        for key in data.dtype.names:
            returnDict[key]=data[key]
        return returnDict
    except:
        print 'Unable to get table column descriptions'
        return {'data':data}

    
