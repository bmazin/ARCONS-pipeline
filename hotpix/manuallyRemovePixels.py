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
from PyQt4 import QtCore
from PyQt4 import QtGui

manualHotPixelReason = tm.timeMaskReason['manual hot pixel']
manualColdPixelReason = tm.timeMaskReason['manual cold pixel']

def removePixel(timeMaskPath,pixelRow,pixelCol,reason='manual hot pixel',nodePath='/',timeInterval=None):
    timeMaskFile = tables.openFile(timeMaskPath,mode='a')

        

    hotPixelsOb = readHotPixels(timeMaskFile)
    expTime = hotPixelsOb.expTime
    reasonEnum = hotPixelsOb.reasonEnum

    if not timeInterval is None:
        if len(timeInterval) == 2:
            tBegin = timeInterval[0]*hotPixelsOb.ticksPerSec
            tEnd = timeInterval[1]*hotPixelsOb.ticksPerSec
        else:
            raise TypeError('timeInterval needs to be a tuple (tBegin,tEnd), instead given:{}'.format(timeInterval))
    else:
        tBegin = 0.
        tEnd = expTime*hotPixelsOb.ticksPerSec

    try:
        reasonValue = reasonEnum[reason]
    except IndexError:
        outStr = 'reason given for removing pixel,{}, is not in the enum embedded in the existing hot pixel file'.format(reason)
        raise TypeError(outStr)

    tableName = constructDataTableName(x=pixelCol, y=pixelRow)
    eventListTable = timeMaskFile.getNode(nodePath + dataGroupName, name=tableName)

    newRow = eventListTable.row
    newRow['tBegin'] = tBegin
    newRow['tEnd'] = tEnd
    newRow['reason'] = reasonValue
    newRow.append()

    timeMaskFile.flush()
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

class EditTimeMaskDialog(QtGui.QDialog):
    def __init__(self,timeMaskPath,pixelRow,pixelCol,parent=None,nodePath='/'):
        super(EditTimeMaskDialog,self).__init__(parent=parent)
        self.parent=parent

        self.timeMaskFile = tables.openFile(timeMaskPath,mode='a')

        self.hotPixelsOb = readHotPixels(self.timeMaskFile)
        self.expTime = self.hotPixelsOb.expTime
        self.reasonEnum = self.hotPixelsOb.reasonEnum

        self.reasonList = self.hotPixelsOb.reasons[pixelRow,pixelCol]
        self.intervalList = self.hotPixelsOb.intervals[pixelRow,pixelCol]
        self.nEntries = len(self.reasonList)

        tableName = constructDataTableName(x=pixelCol, y=pixelRow)
        self.eventListTable = self.timeMaskFile.getNode(nodePath + dataGroupName, name=tableName)

        self.initUI()
        self.exec_()
        
    def closeEvent(self,evt):
        try:
            self.timeMaskFile.flush()
            self.timeMaskFile.close()
        except:
            pass
        super(EditTimeMaskDialog,self).closeEvent(evt)

    def getSelectedIntervalIndices(self):
        return np.where([checkbox.isChecked() for checkbox in self.checkboxes])[0]

    def removeSelectedIntervals(self):
        #find the rows where the reason flag matches one of those in the reasons list
        rowIndicesToRemove = self.getSelectedIntervalIndices()

        #as we start iterating through and deleting rows, the row numbers will change
        #so we need to correct the row indices accordingly
        #subtract from each row the number of rows deleted before it
        correctedRowIndices = rowIndicesToRemove - np.arange(len(rowIndicesToRemove))
        #delete rows
        for eachRowIdx in correctedRowIndices:
            self.eventListTable.removeRows(eachRowIdx)
        self.close()

    def initUI(self):
        grid = QtGui.QGridLayout()
        self.checkboxes = []
        for iEntry,(reasonIndex,inter) in enumerate(zip(self.reasonList,self.intervalList)):
            reason = self.reasonEnum(reasonIndex)
            checkbox = QtGui.QCheckBox()
            start = [eachComponent[0][0] for eachComponent in inter.components][0]
            stop = [eachComponent[0][1] for eachComponent in inter.components][0]
            intervalLabel = QtGui.QLabel('({},{})'.format(start,stop))
            reasonLabel = QtGui.QLabel('{}'.format(reason))
            grid.addWidget(checkbox,iEntry,0)
            grid.addWidget(intervalLabel,iEntry,1)
            grid.addWidget(reasonLabel,iEntry,2)
            self.checkboxes.append(checkbox)
        self.button_removeSelected = QtGui.QPushButton('Remove Selected Intervals')
        #grid.addWidget(self.button_removeSelected,self.nEntries,0,1,3)

        scrollWidget = QtGui.QWidget()
        scrollWidget.setLayout(grid)
        scrollArea = QtGui.QScrollArea()
        scrollArea.setWidget(scrollWidget)
        box = QtGui.QVBoxLayout()
        box.addWidget(scrollArea)
        box.addWidget(self.button_removeSelected)
        self.setLayout(box)
        self.connect(self.button_removeSelected,QtCore.SIGNAL('clicked()'), self.removeSelectedIntervals)
            
            
if __name__ == "__main__":
    timeMaskPath = '/Scratch/timeMasks/20141020/timeMask_20141021-091336.h5'
    pixelToRemove = (0,0)
    #removePixel(timeMaskPath,pixelToRemove[0],pixelToRemove[1],reason='unknown')
    #removePixel(timeMaskPath,pixelToRemove[0],pixelToRemove[1])
    unremovePixel(timeMaskPath,pixelToRemove[0],pixelToRemove[1],reasons=['unknown'])

    

