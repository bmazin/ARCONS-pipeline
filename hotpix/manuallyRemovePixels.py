'''
Author: Matt Strader
Date: Jan 5 2015

manually remove a pixel for the duration of an ObsFile by adding it to the hotpix cal file associated with the ObsFile
'''

import numpy as np
from interval import interval
import tables
import headers.TimeMask as tm
from hotPixelMasker import hotPixelMasker
from util.popup import plotArray,pop
from hotpix.hotPixels import constructDataTableName,dataGroupName,readHotPixels

manualHotPixelReason = tm.timeMaskReason['manual hot pixel']
manualColdPixelReason = tm.timeMaskReason['manual cold pixel']

def removePixel(timeMaskPath,pixelRow,pixelCol,reason='manual hot pixel',nodePath='/'):
    timeMaskFile = tables.openFile(timeMaskPath,mode='a')

    hotPixelsOb = readHotPixels(timeMaskFile)
    expTime = hotPixelsOb.expTime
    reasonEnum = hotPixelsOb.reasonEnum
    try:
        reasonValue = reasonEnum[reason]
    except IndexError:
        outStr = 'reason given for removing pixel,{}, is not in the enum embedded in the existing hot pixel file'.format(reason)
        raise IndexError(outStr)

    tableName = constructDataTableName(x=pixelCol, y=pixelRow)
    eventListTable = timeMaskFile.getNode(nodePath + dataGroupName, name=tableName)

    print hotPixelsOb.reasons,hotPixelsOb.reasonEnum
    print eventListTable.read()
    newRow = eventListTable.row
    newRow['tBegin'] = 0.
    newRow['tEnd'] = expTime*hotPixelsOb.ticksPerSec
    newRow['reason'] = reasonValue
    newRow.append()

    timeMaskFile.flush()
    print eventListTable.read()
    timeMaskFile.close()

def unremovePixel(timeMaskPath,pixelRow,pixelCol,reasons=['manual hot pixel','manual cold pixel'],nodePath='/'):
    timeMaskFile = tables.openFile(timeMaskPath,mode='a')

    hotPixelsOb = readHotPixels(timeMaskFile)
    expTime = hotPixelsOb.expTime
    reasonEnum = hotPixelsOb.reasonEnum
    try:
        reasonValues = [reasonEnum[eachReason] for eachReason in reasons]
    except IndexError:
        outStr = 'reason given for removing pixel,{}, is not in the timeMaskReason enum'.format(reason)
        raise IndexError(outStr)

    tableName = constructDataTableName(x=pixelCol, y=pixelRow)
    eventListTable = timeMaskFile.getNode(nodePath + dataGroupName, name=tableName)
    print eventListTable.read()
    
    #find the rows where the reason flag matches one of those in the reasons list
    rowIndicesWithGivenReasons = np.where(np.in1d(eventListTable.cols.reason,reasonValues))[0]
    print rowIndicesWithGivenReasons

    #as we start iterating through and deleting rows, the row numbers will change
    #so we need to correct the row indices accordingly
    #subtract from each row the number of rows deleted before it
    rowIndicesWithGivenReasons.sort()
    correctedRowIndices = rowIndicesWithGivenReasons - np.arange(len(rowIndicesWithGivenReasons))
    #delete rows
    for eachRowIdx in correctedRowIndices:
        eventListTable.removeRows(eachRowIdx)

    timeMaskFile.flush()

    print eventListTable.read()

    timeMaskFile.close()

if __name__ == "__main__":
    timeMaskPath = '/Scratch/timeMasks/20141020/timeMask_20141021-091336.h5'
    pixelToRemove = (0,0)
    #removePixel(timeMaskPath,pixelToRemove[0],pixelToRemove[1],reason='unknown')
    #removePixel(timeMaskPath,pixelToRemove[0],pixelToRemove[1])
    unremovePixel(timeMaskPath,pixelToRemove[0],pixelToRemove[1],reasons=['unknown'])

    

