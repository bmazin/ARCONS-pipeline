#!/bin/python

'''
Author: Paul Szypryt		Date: July 11, 2013
'''

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util.readDict import readDict
from util import utils
import tables
import matplotlib.pyplot as plt
import hotpix.hotPixels as hp
import os
from time import time
import sys
import LoadImageStack
from tables import *

from PyQt4.QtGui import *
from PyQt4.QtGui import *
from DisplayStack_gui import Ui_DisplayStack_gui

class DisplayStack(QMainWindow):

    def __init__(self):
        # Start up gui
        QWidget.__init__(self, parent=None)
        self.ui = Ui_DisplayStack_gui()
        self.ui.setupUi(self)

        # Initialize Variables
        self.initializeVariables()

        # Lists and buttons used for specifying run and target information.
        # Click item in runList to select a run. Run number corresponds to array index.  Load up corresponding target list.
        self.ui.runList.itemClicked.connect(self.selectRun)
        # Click item in targetList to select a target.
        self.ui.targetList.itemClicked.connect(self.selectTarget)
        # Click button in 
        self.ui.targetButton.clicked.connect(self.loadTarget)


        self.ui.sunsetList.itemClicked.connect(self.createObsList)


        # Use wavelength calibration checkbox
        self.ui.wavelengthCalibrationBox.clicked.connect(self.useWaveCal)
        # Use flat calibration checkbox
        self.ui.flatCalibrationBox.clicked.connect(self.useFlatCal)

        # Buttons for obs list creation
        self.ui.addButton.clicked.connect(self.addObservation)
        self.ui.removeButton.clicked.connect(self.removeObservation)
        self.ui.clearButton.clicked.connect(self.createObsList)

        # Start process button
        self.ui.stackButton.clicked.connect(self.stackProcess)


        # Load image stack button
        self.ui.loadStackButton.clicked.connect(self.chooseStack)

    # Initialize Variables
    def initializeVariables(self):
        # Define path names
        self.displayStackPath = '/Scratch/DisplayStack/'
        self.defaultWavelengthPath = '/Scratch/waveCalSolnFiles/'
        self.defaultFlatPath = '/Scratch/flatCalSolnFiles/'
       
        # Load and display list of run names from /Scratch/DisplayStack/runList.dict
        self.loadRunData()
        
        # Load list of target names from /Scratch/DisplayStack/runName/runName.dict, for all runs
        self.loadTargetData()


    # Function to load list of run names. Runs on initialization.
    def loadRunData(self):
        # Load run data from /Scratch/DisplayStack/runList.dict
        self.runData = readDict()
        self.runData.read_from_file(self.displayStackPath + '/runList.dict')
        self.runNames = np.array(self.runData['runs'])
        # Populate runList table with run names.
        for iRun in range(len(self.runNames)):
            self.ui.runList.addItem(self.runNames[iRun])

    # Function to load a table of target names.
    def loadTargetData(self):
        self.targetNames = []
        # Cycle through runs and extract target name information from various dictionaries.
        for iRun in range(len(self.runNames)):
            self.targetData = readDict()
            self.targetData.read_from_file(self.displayStackPath + self.runNames[iRun] + '/' + self.runNames[iRun] + '.dict')
            self.iTargets = np.array(self.targetData['targets'])
            self.targetNames.append(self.iTargets)

    # Function to select a run and populate the target list for that particular run.
    def selectRun(self):
        # Clear list of target information for previously selected run.
        self.ui.targetList.clear()
        # Define a run number by the index of the selected run.
        self.runNumber = self.ui.runList.row(self.ui.runList.currentItem())
        # Populate targetList table with target names for selected run.
        for iTarget in range(len(self.targetNames[self.runNumber])):
            self.ui.targetList.addItem(self.targetNames[self.runNumber][iTarget])

    def selectTarget(self):
        self.targetNumber = self.ui.targetList.row(self.ui.targetList.currentItem())

    def loadTarget(self):
        self.run = self.runNames[self.runNumber]
        self.target = self.targetNames[self.runNumber][self.targetNumber]
        try:
            self.paramName = self.displayStackPath + self.run + '/' + self.target + '/' + self.target + '.dict'
            self.paramData = readDict()
            self.paramData.read_from_file(self.paramName)
            self.obsTimes = np.array(self.paramData['obsTimes'])
            self.utcDates = self.paramData['utcDates']
            self.sunsetDates = self.paramData['sunsetDates']
            self.calTimestamps = self.paramData['calTimestamps']
            self.flatCalDates = self.paramData['flatCalDates']
            self.RA = self.paramData['RA']
            self.Dec = self.paramData['Dec']
            self.hourAngleOffset = self.paramData['HA_offset']

            print 'Loading parameter file at ' + self.displayStackPath + self.run + '/' + self.target + '/' + self.target + '.dict'
            self.createSunsetList()
            #self.createObsList()
            self.createWavelengthList()
            self.createFlatList()
            self.paramFileExists = True
        except IOError:           
            print 'No existing parameter file at ' + self.displayStackPath + self.run + '/' + self.target + '/' + self.target + '.dict'
            self.ui.sunsetList.clear()
            self.ui.obsList.clear()
            self.ui.inputList.clear()
            self.ui.wavelengthList.clear()
            self.ui.flatList.clear()
            self.paramFileExists = False
        
    # Choose Obs File
    # Create list of available sunset dates
    def createSunsetList(self):
        self.ui.sunsetList.clear()
        for iDate in range(len(self.sunsetDates)):
            self.ui.sunsetList.addItem(self.sunsetDates[iDate])
    # Create Initial Obs file list
    def createObsList(self):
        self.ui.obsList.clear()
        self.ui.inputList.clear()
        self.currentSunsetDate = self.sunsetDates[self.ui.sunsetList.row(self.ui.sunsetList.currentItem())]
        self.currentUTCDate = self.utcDates[self.ui.sunsetList.row(self.ui.sunsetList.currentItem())]
        self.singleDayObservations = self.obsTimes[self.ui.sunsetList.row(self.ui.sunsetList.currentItem())]
        for iObs in range(len(self.singleDayObservations)):
            self.ui.obsList.addItem(self.singleDayObservations[iObs])
    # Add observation to input list
    def addObservation(self):
        self.selectedObs = self.ui.obsList.currentItem()
        self.ui.obsList.takeItem(self.ui.obsList.row(self.selectedObs))
        self.ui.inputList.addItem(self.selectedObs)
        self.ui.inputList.sortItems()
    # Remove observation from input list
    def removeObservation(self):
        self.removedObs = self.ui.inputList.currentItem()
        self.ui.inputList.takeItem(self.ui.inputList.row(self.removedObs))
        self.ui.obsList.addItem(self.removedObs)
        self.ui.obsList.sortItems()

    # Load settings
    def loadSettings(self):
        
        self.validSettings = True        

        # Run and target information
        self.run = self.runNames[self.runNumber]
        self.target = self.targetNames[self.runNumber][self.targetNumber]

        # General settings
        self.integrationTime = int(self.ui.integrationTimeLine.text())
        self.useTimeAdjustment = self.ui.timeAdjustmentBox.isChecked()
        self.useHotPixelMasking = self.ui.hotPixelBox.isChecked()

        # Wavelength calibration settings
        self.useWavelengthCalibration = self.ui.wavelengthCalibrationBox.isChecked()
        self.useBestWavelengthCalibration = self.ui.bestWavelengthCalibrationBox.isChecked()
        self.lowerWavelengthCutoff = float(self.ui.lowerWavelengthCutoffLine.text())
        self.upperWavelengthCutoff = float(self.ui.upperWavelengthCutoffLine.text())
               
        # Flat calibration settings
        self.useFlatCalibration = self.ui.flatCalibrationBox.isChecked()
        self.useDeadPixelMasking = self.ui.deadPixelBox.isChecked() 

        self.fileCount = self.ui.inputList.count()  

        self.weighted = self.useFlatCalibration
        self.useRawCounts = not (self.useWavelengthCalibration and self.useFlatCalibration)
        self.scaleByEffInt = self.useHotPixelMasking
            
        if self.ui.sunsetList.currentItem() != None:

            if self.fileCount == 0:
                print 'Please select files to process...'
                self.validSettings = False

            if self.useWavelengthCalibration: 
                if self.useBestWavelengthCalibration:
                    print 'Using best wavelength calibration...'
                elif self.ui.wavelengthList.currentItem() == None:
                    print 'Please select wavelength calibration...'
                    self.validSettings = False
                else:
                    self.selectedWvlCal = self.ui.wavelengthList.currentItem().text()
                    self.wvlCalFilename = str(FileName(run=self.run,date=self.currentSunsetDate,tstamp=self.selectedWvlCal).calSoln())

            if self.useFlatCalibration:
                if self.ui.flatList.currentItem() == None:
                    print 'Please select flat calibration...'
                    self.validSettings = False
                else:
                    self.flatCalNight = self.ui.flatList.currentItem().text()
                    self.flatCalFilename = str(FileName(run=self.run,date=self.flatCalNight).flatSoln())
                    
        else:
            print 'Please select sunset night...'
            self.validSettings = False

    # Load hot pixel mask
    def loadHotMask(self):
        self.hotPixelFilename = str(self.displayStackPath + self.run + '/' + self.target + '/HotPixelMasks/hotPixelMask_' + self.obsTS + '.h5')
        if not os.path.exists(self.hotPixelFilename):
            hp.findHotPixels(obsFile=self.ob,outputFileName=self.hotPixelFilename)
            print "Hot pixel mask saved to %s"%(self.hotPixelFilename)
        self.ob.loadHotPixCalFile(self.hotPixelFilename,switchOnMask=True)

    # Create wavelength cal file list
    def createWavelengthList(self):
        self.ui.wavelengthList.clear()
        for iCal in range(len(self.calTimestamps)):
            self.ui.wavelengthList.addItem(self.calTimestamps[iCal])
    # Enable/disable wavecal options
    def useWaveCal(self):
        if self.ui.wavelengthCalibrationBox.isChecked():
            self.ui.lowerWavelengthCutoffLine.setEnabled(True)
            self.ui.lowerWavelengthCutoffLabel.setEnabled(True)
            self.ui.upperWavelengthCutoffLine.setEnabled(True)
            self.ui.upperWavelengthCutoffLabel.setEnabled(True)
            self.ui.wavelengthList.setEnabled(True)          
            self.ui.flatCalibrationBox.setEnabled(True)
            self.ui.flatCalibrationBox.setChecked(True)
            self.ui.deadPixelBox.setEnabled(True)
            self.ui.deadPixelBox.setChecked(True)       
            self.ui.flatList.setEnabled(True)
            self.ui.bestWavelengthCalibrationBox.setEnabled(True)
            self.ui.bestWavelengthCalibrationBox.setChecked(True)
        else:
            self.ui.lowerWavelengthCutoffLine.setEnabled(False)
            self.ui.lowerWavelengthCutoffLabel.setEnabled(False)
            self.ui.upperWavelengthCutoffLine.setEnabled(False)
            self.ui.upperWavelengthCutoffLabel.setEnabled(False)
            self.ui.wavelengthList.setEnabled(False)
            self.ui.flatCalibrationBox.setEnabled(False)
            self.ui.flatCalibrationBox.setChecked(False)
            self.ui.deadPixelBox.setChecked(False)
            self.ui.deadPixelBox.setEnabled(False)
            self.ui.flatList.setEnabled(False)
            self.ui.bestWavelengthCalibrationBox.setEnabled(False)
            self.ui.bestWavelengthCalibrationBox.setChecked(False)


    # Create flat cal file list
    def createFlatList(self):
        self.ui.flatList.clear()
        for iCal in range(len(self.flatCalDates)):
            self.ui.flatList.addItem(self.flatCalDates[iCal])
    # Enable/disable flatcal options
    def useFlatCal(self):
        if self.ui.flatCalibrationBox.isChecked():
            self.ui.deadPixelBox.setEnabled(True)
            self.ui.deadPixelBox.setChecked(True)       
            self.ui.flatList.setEnabled(True)
        else:
            self.ui.deadPixelBox.setChecked(False)
            self.ui.deadPixelBox.setEnabled(False)
            self.ui.flatList.setEnabled(False)

    # Load dead pixel mask
    def loadDeadMask(self):
        self.deadPixelFilename = str(self.displayStackPath + self.run + '/' + self.target + '/DeadPixelMasks/deadPixelMask_' + self.obsTS + '.npz')
        if not os.path.exists(self.deadPixelFilename):
            self.deadMask = self.ob.getDeadPixels()
            np.savez(self.deadPixelFilename,deadMask = self.deadMask)
            print "Dead pixel mask saved to %s"%(self.deadPixelFilename)
        else:
            self.deadFile = np.load(self.deadPixelFilename)
            self.deadMask = self.deadFile['deadMask']

    # Describe the structure of the header row
    class headerDescription(tables.IsDescription):
        targetName = tables.StringCol(100, dflt='')
        run = tables.StringCol(100, dflt='')
        obsFileName = tables.StringCol(100, dflt='')
        wvlCalFileName = tables.StringCol(100, dflt=np.nan)
        flatCalFileName = tables.StringCol(100, dflt='')
        deadPixFileName = tables.StringCol(100, dflt='')
        hotPixFileName = tables.StringCol(100, dflt='')
        nCol = tables.UInt32Col(dflt=-1)
        nRow = tables.UInt32Col(dflt=-1)
        lowWvlCutoff = tables.Float64Col(dflt=np.nan)
        highWvlCutoff = tables.Float64Col(dflt=np.nan)
        exptime = tables.Float64Col(dflt=np.nan)
        lst = tables.StringCol(100, dflt='')
        integrationTime = tables.Float64Col(dflt=np.nan)
        RA = tables.StringCol(100, dflt='')
        Dec = tables.StringCol(100, dflt='')
        HA_offset = tables.Float64Col(dflt=0.0)


    # Create output name
    def createOutputName(self):
        self.rawName = str(self.displayStackPath + self.run + '/' + self.target + '/ImageStacks/' + 'ImageStack_' + self.obsTS + '_' + str(self.integrationTime) + 's')
        if self.useWavelengthCalibration and self.useHotPixelMasking:
            self.outputFilename = str(self.rawName + '_' + str(int(self.lowerWavelengthCutoff)) + '-' + str(int(self.upperWavelengthCutoff)) + '_hp.h5')
        elif self.useWavelengthCalibration and not self.useHotPixelMasking:
            self.outputFilename = str(self.rawName + '_' + str(int(self.lowerWavelengthCutoff)) + '-' + str(int(self.upperWavelengthCutoff)) + '.h5')
        elif not self.useWavelengthCalibration and self.useHotPixelMasking:
            self.outputFilename = str(self.rawName + '_hp.h5')
        else:
            self.outputFilename = str(self.rawName + '.h5')

    def createH5File(self):
        # Create header and data group and table names
        headerGroupName = 'header'
        headerTableName = 'header'
        dataGroupName = 'stack'
        dataTableName = 'stack'
        timeTableName = 'time'

        # Create lookup names for header information
        runColName = 'run'
        targetColName = 'targetName'
        obsFileColName = 'obsFileName'
        wvlCalFileColName = 'wvlCalFileName'
        flatCalFileColName = 'flatCalFileName'
        nRowColName = 'nRow'
        nColColName = 'nCol'
        RAColName = 'RA'
        DecColName = 'Dec'
        deadPixColName = 'deadPixFileName'
        hotPixColName = 'hotPixFileName'
        lowWvlColName = 'lowWvlCutoff'
        highWvlColName = 'highWvlCutoff'
        expTimeColName = 'exptime'
        lstColName = 'lst'
        integrationTimeColName = 'integrationTime'
        HA_offsetColName = 'HA_offset'


        # Create and h5 output file, create header and data groups
        fileh = tables.openFile(self.outputFilename, mode='w')
        headerGroup = fileh.createGroup("/", headerGroupName, 'Header')
        stackGroup = fileh.createGroup("/", dataGroupName, 'Image Stack')

        # Create row for header information
        headerTable = fileh.createTable(headerGroup, headerTableName, self.headerDescription,
                                        'Header Info')
        header = headerTable.row

        # Fill in the header with possibly useful information.
        header[runColName] = self.run
        header[targetColName] = self.target
        header[obsFileColName] = self.obsFn       
        header[nColColName] = self.numberCols
        header[nRowColName] = self.numberRows
        header[RAColName] = self.RA
        header[DecColName] = self.Dec
        header[expTimeColName] = self.exptime
        header[lstColName] = self.lst
        header[integrationTimeColName] = self.integrationTime
        header[HA_offsetColName] = self.hourAngleOffset
        if self.useDeadPixelMasking:
            header[deadPixColName] = self.deadPixelFilename
        if self.useHotPixelMasking:
            header[hotPixColName] = self.hotPixelFilename
        if self.useWavelengthCalibration:
            header[wvlCalFileColName] = self.ob.wvlCalFileName       
            header[lowWvlColName] = self.lowerWavelengthCutoff
            header[highWvlColName] = self.upperWavelengthCutoff
        if self.useFlatCalibration:
            header[flatCalFileColName] = self.flatCalFilename
        header.append()

        # Create an h5 array for the midtime of each frame in the image cube.
        timeTable = fileh.createCArray(stackGroup, timeTableName, Float64Atom(), (1,len(self.times)))       
        timeTable[:] = self.times
        
        # Create an h5 table for the image cube.
        stackTable = fileh.createCArray(stackGroup, dataTableName, Float64Atom(), (self.numberRows,self.numberCols, self.cube.shape[2]))        
        stackTable[:] = self.cube

        # Flush the h5 output file
        fileh.flush()
        fileh.close()


    # Start process for creating image stacks
    def stackProcess(self):

        # Check for valid params file
        if self.paramFileExists:

            # Load settings choosen from gui
            self.loadSettings()
            if self.validSettings:
        
                # Loop through all files in input list
                for iFile in range(self.fileCount):
                    # Create ObsFile instance
                    self.obsTS = str(self.currentUTCDate) + '-' + self.ui.inputList.item(iFile).text()
                    self.obsFn = str(FileName(run=self.run,date=self.currentSunsetDate,tstamp=self.obsTS).obs())
                    print 'Processing file ' + self.obsFn + '...'
                    self.ob = ObsFile(self.obsFn)

                    self.numberRows = self.ob.nRow
                    self.numberCols = self.ob.nCol

                    # Load time adjustment file
                    if self.useTimeAdjustment:
                        print 'Loading time adjustment file...'
                        self.ob.loadTimeAdjustmentFile(FileName(run=self.run).timeAdjustments())
       
                    # Load hot pixel mask
                    if self.useHotPixelMasking:
                        print 'Loading hot pixel mask...'
                        self.loadHotMask()

                    # Load wave cal solution
                    if self.useWavelengthCalibration:
                        if self.useBestWavelengthCalibration:
                            print 'Loading best wavelength calibration...'
                            self.ob.loadBestWvlCalFile()
                        else:
                            print 'Loading selected wavelength calibration...'
                            self.ob.loadWvlCalFile(self.wvlCalFilename)
                        
                    # Load flatcal solution
                    if self.useFlatCalibration:
                        print 'Loading flat calibration...'
                        self.ob.loadFlatCalFile(self.flatCalFilename)

                    # Load dead pixel mask
                    if self.useDeadPixelMasking:
                        print 'Loading dead pixel mask...'
                        self.loadDeadMask()

                    # Set wavelength cutoffs
                    if self.useWavelengthCalibration:
                        print 'Setting wavelength cutoffs...'
                        self.ob.setWvlCutoffs(self.lowerWavelengthCutoff,self.upperWavelengthCutoff)
        
                    # Timing
                    self.unix = self.ob.getFromHeader('unixtime')
                    self.startJD = self.unix/86400.+2440587.5
                    self.exptime = self.ob.getFromHeader('exptime')
                    self.lst = self.ob.getFromHeader('lst')
                    
                    self.times = []
                    self.frames = []
                        
                    # Create Image Stack
                    print 'Stacking images...'
                    for iSec in np.arange(0,self.exptime,self.integrationTime):
                        #add seconds offset to julian date, move jd to center of bin
                        self.jd = self.startJD + iSec/(24.*3600.) + self.integrationTime/2./(24.*3600.)
                        self.times.append(self.jd)
                        print 'Creating frame for time ' + str(self.jd)
                        self.frameData = self.ob.getPixelCountImage(firstSec=iSec,integrationTime=self.integrationTime,weighted=self.weighted,getRawCount=self.useRawCounts,scaleByEffInt=self.scaleByEffInt)
                        self.frame = self.frameData['image']

                        if self.ui.verticalFlipBox.isChecked():
                            self.frame = np.flipud(self.frame)                   
         
                        if self.useDeadPixelMasking:
                            self.frame[self.deadMask == 0] = np.nan
                        self.frames.append(self.frame)

                    self.cube = np.dstack(self.frames)
                    self.times = np.array(self.times)
            
                    # Create output file
                    self.createOutputName()
                    print 'Saving image stack to ' + self.outputFilename
                    self.createH5File()                 

        # Invalid params file
        else:
            print 'Invalid parameter file...'

    # Choose an image stack
    
    def chooseStack(self):
        self.defaultLoadStackDirectory = str(self.displayStackPath)
        self.stackName = ''
        self.stackName = QFileDialog.getOpenFileName(parent=None, directory=self.defaultLoadStackDirectory, caption=str("Choose Image Stack"), filter=str("H5 (*.h5)")) 
        if self.stackName == '':
            print 'No file chosen'
        else:          
            loadStackApp = LoadImageStack.LoadImageStack(stackName = self.stackName)
            loadStackApp.show()
            loadStackApp.exec_()
    

# Start up main gui
if __name__ == "__main__":
	app = QApplication(sys.argv)
	myapp = DisplayStack()
	myapp.show()
	app.exec_()



