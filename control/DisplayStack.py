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
            self.paramData = readDict()
            self.paramData.read_from_file(self.displayStackPath + self.run + '/' + self.target + '/' + self.target + '.dict')
            self.obsTimes = np.array(self.paramData['obsTimes'])
            self.utcDates = self.paramData['utcDates']
            self.sunsetDates = self.paramData['sunsetDates']
            self.calTimestamps = self.paramData['calTimestamps']
            self.flatCalDates = self.paramData['flatCalDates']  

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
                if self.ui.wavelengthList.currentItem() == None:
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
            hp.findHotPixels(self.obsFn,self.hotPixelFilename)
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

    # Create output name
    def createOutputName(self):
        self.rawName = str(self.displayStackPath + self.run + '/' + self.target + '/ImageStacks/' + 'ImageStack_' + self.obsTS + '_' + str(self.integrationTime) + 's')
        if self.useWavelengthCalibration and self.useHotPixelMasking:
            self.outputFilename = str(self.rawName + '_' + str(int(self.lowerWavelengthCutoff)) + '-' + str(int(self.upperWavelengthCutoff)) + '_hp.npz')
        elif self.useWavelengthCalibration and not self.useHotPixelMasking:
            self.outputFilename = str(self.rawName + '_' + str(int(self.lowerWavelengthCutoff)) + '-' + str(int(self.upperWavelengthCutoff)) + '.npz')
        elif not self.waveLengthCalibration and self.useHotPixelMasking:
            self.outputFilename = str(self.rawName + '_hp.npz')
        else:
            self.outputFilename = str(self.rawName + '.npz')

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
                        print 'Loading wavelength calibration...'
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
                        if self.useDeadPixelMasking:
                            self.frame[self.deadMask == 0] = np.nan
                        self.frames.append(self.frame)

                    self.cube = np.dstack(self.frames)
                    self.times = np.array(self.times)
            
                    # Create output file
                    self.createOutputName()
                    print 'Saving image stack to ' + self.outputFilename
                    np.savez(self.outputFilename, stack=self.cube, jd=self.times)

        # Invalid params file
        else:
            print 'Invalid parameter file...'

    # Choose an image stack
    
    def chooseStack(self):
        self.defaultLoadStackDirectory = str(self.displayStackPath)
        #self.defaultLoadStackDirectory = str(self.displayStackPath + self.run + '/' + self.target + '/ImageStacks')
        self.stackName = ''
        self.stackName = QFileDialog.getOpenFileName(parent=None, directory=self.defaultLoadStackDirectory, caption=str("Choose Image Stack"), filter=str("NPZ (*.npz)")) 
        if self.stackName == '':
            print 'No file chosen'
        else:          
            loadStackApp = LoadImageStack.LoadImageStack(stackName = self.stackName)
            loadStackApp.show()
            loadStackApp.exec_()
            # Change to a QDialog...
    

# Start up main gui
if __name__ == "__main__":
	app = QApplication(sys.argv)
	myapp = DisplayStack()
	myapp.show()
	app.exec_()



