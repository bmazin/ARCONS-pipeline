'''
Author: Alex Walter
Date: Feb 2, 2015

This file contains the keywords and descriptions for image stack h5 files. 
Major functions:
    writeImageStack() - writes list of images (and other info) to file
    readImageStack() - reads list of images from file
'''


import os, warnings
import tables
import numpy as np

# Data group and table names
dataGroupName = 'ImageData'
dataTableName = 'ImageData'
headerGroupName = 'Header'
headerTableName = 'Header'
    
def imageStackHeaderDescription(nObs):
    '''
    Returns a dictionary for the header table description
    
    nObs - Number of obsFiles used to make image stack. Same as length of tStamps list
    '''
    description = {
        'targetName' : tables.StringCol(100, dflt=''),           # name of target observed
        'run' : tables.StringCol(100, dflt=''),                  # telescope run. ie., 'PAL2014'
        'nFrames' : tables.UInt32Col(dflt=0),                           # Number of frames in image stack
        'wvlLowerLimit' : tables.Float64Col(dflt=np.nan),               # short wavelength limit of photons used in images
        'wvlUpperLimit' : tables.Float64Col(dflt=np.nan),               # long wavelength limit
        'weighted' : tables.UInt8Col(dflt=False),                       # True if flat cal applied
        'fluxWeighted' : tables.UInt8Col(dflt=False),                   # True if flux cal applied
        'hotPixelsApplied' : tables.UInt8Col(dflt=False),               # True if hotpix was applied
        'maxExposureTime' : tables.Float64Col(dflt=0.),                 # dt in ObsFileSeq.__init__()
        'tStamps' : tables.StringCol(15, dflt='', shape=nObs)}   # List of time stamps of obsfiles used. 15 chars each, ie., '20140923-111213'
        #Could add object's RA and Dec. HA offset. Etc..
    
    return description
    
def imageStackDataDescription(nRows,nCols):
    '''
    Returns a dictionary for the data table description
    
    nRows and nCols describe the shape of the image in pixels
    '''
    description = {
        'images' : tables.Float64Col(dflt=np.nan, shape=(nRows,nCols)),          #Column for 2D image [Total counts during exposure (scaled up to intTime)]
        'pixIntTimes' : tables.Float64Col(dflt=0., shape=(nRows,nCols)),         #Column for 2D pixel integration times [seconds]
        'startTimes' : tables.Float64Col(),                                      #Column for start times [Julian Date]
        'endTimes' : tables.Float64Col(),                                        #Column for end times [Julian Date]
        'intTimes' : tables.Float64Col()}                                        #Column for total exposure [seconds]
        #Could probably add RA/Dec info here too
    
    return description
    
    
def writeImageStack(fullFilename, images, startTimes, **kwargs):
                    #endTimes, intTimes, pixIntTimes, 
                    #targetName,run,nFrames,wvlLowerLimit,wvlUpperLimit,weighted,fluxWeighted,hotPixelsApplied,maxExposureTime,obsTStamps):
    '''
    Saves an h5 file containing the image stacks. See headers.DisplayStackHeaders for variable definitions
    
    Inputs:
        fullFilename - name of imagestack file
        images - list of images. images[number, row, col]
        startTimes - list of startTimes [Julian Date]
        **kwargs --> this info is optional. Valid keywords:
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
    # Create an h5 output file, create header and data groups
    if os.path.isfile(fullFilename):
        warnings.warn("Overwriting image stack file: "+str(fullFilename),UserWarning)
    fileh = tables.openFile(fullFilename, mode='w')
    
    try:
        nFrames = len(images)
        nObs = len(kwargs.get('tStamps',[0]))
        (nRows,nCols) = np.asarray(images[0]).shape
        
        # Create table for header information
        headerGroup = fileh.createGroup("/", headerGroupName, 'Header')
        headerTableDescription = imageStackHeaderDescription(nObs)
        headerTable = fileh.createTable(headerGroup, headerTableName, headerTableDescription,'Header Info')
        header = headerTable.row
        kwargs['nFrames'] = nFrames
        for keyword in set(kwargs.keys()).intersection(headerTableDescription.keys()):
            #grab all the correct info from the **kwargs parameters
            header[keyword] = kwargs[keyword]
        header.append()
        headerTable.flush()
        headerTable.close()
            
        # Create a table for the image data
        dataGroup = fileh.createGroup("/", dataGroupName, 'Image Data')
        dataTableDescription = imageStackDataDescription(nRows,nCols)
        filt = tables.Filters(complevel=2, complib='zlib') #See doc for pyTables 2.3 http://pytables.readthedocs.org/en/v.2.3.1/usersguide/libref.html#filtersclassdescr
        dataTable = fileh.createTable(where=dataGroup, name=dataTableName, description=dataTableDescription, title='Image Data',
                                      filters=filt, expectedrows=nFrames)
        pixIntTimes = kwargs.get('pixIntTimes',np.zeros((nRows,nCols)))
        endTimes = kwargs.get('endTimes',np.zeros(nFrames))
        intTimes = kwargs.get('intTimes',np.zeros(nFrames))
        for i in range(nFrames):
            row = dataTable.row
            row['images'] = images[i]
            row['startTimes'] = startTimes[i]
            row['pixIntTimes'] = pixIntTimes[i]
            row['endTimes'] = endTimes[i]
            row['intTimes'] = intTimes[i]
            row.append()
        dataTable.flush()
        dataTable.close()
        
    finally:
        # Flush the h5 output file
        fileh.flush()
        fileh.close()
    
def readImageStack(fileName):
    '''
    Returns dictionary of data in image stack file
    Input:
        fileName - full path to file
    Returns dict with keywords:
        images
        pixIntTimes
        startTimes
        endTimes
        intTimes
        headerInfo
    '''
    print 'Loading image stack from ',fileName
    fileh = tables.openFile(fileName, mode='r')
    try:
        headerTable = fileh.getNode('/' + headerGroupName, headerTableName)
        headerInfo = headerTable.read() #The one row from the header
        
        dataTable = fileh.getNode('/' + dataGroupName, dataTableName).read()
    finally:
        fileh.close()
    
    return {'images':dataTable['images'], 'pixIntTimes':dataTable['pixIntTimes'], 'startTimes':dataTable['startTimes'],
            'endTimes':dataTable['endTimes'],'intTimes':dataTable['intTimes'],'headerInfo':headerInfo}

    
