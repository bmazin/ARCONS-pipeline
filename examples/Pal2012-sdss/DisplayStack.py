#!/bin/python

'''
Author: Paul Szypryt		Date: July 1, 2013
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

from PyQt4.QtGui import *
from PyQt4.QtGui import *
from DisplayStack_gui import Ui_DisplayStack_gui


class DisplayStack(QMainWindow):

    def __init__(self):
        # 0) Start up gui
        QWidget.__init__(self, parent=None)
        self.ui = Ui_DisplayStack_gui()
        self.ui.setupUi(self)

        # Initialize Variables
        self.initializeVariables()
     
        # Initialize run and target, load param file
        self.chooseRun()

        '''
        # Browse wavecal button
        self.ui.wavelengthCalibrationSolutionButton.clicked.connect(self.chooseWavelengthFile)
        # Browse flatcal button
        self.ui.flatCalibrationSolutionButton.clicked.connect(self.chooseFlatFile)       
        '''

        # Use wavelength calibration checkbox
        self.ui.wavelengthCalibrationBox.clicked.connect(self.useWaveCal)
        # Use flat calibration checkbox
        self.ui.flatCalibrationBox.clicked.connect(self.useFlatCal)

        # Buttons for obs list creation
        self.ui.addButton.clicked.connect(self.addObservation)
        self.ui.removeButton.clicked.connect(self.removeObservation)
        self.ui.clearButton.clicked.connect(self.createObsList)

        # Start process button
        self.ui.startButton.clicked.connect(self.startProcess)

    # 1) Initialize Variables
    def initializeVariables(self):
        # Define path names
        self.displayStackPath = '/Scratch/DisplayStack/'
        self.defaultWavelengthPath = '/Scratch/waveCalSolnFiles/'
        self.defaultFlatPath = '/Scratch/flatCalSolnFiles/'

        # Arrays with run names and buttons
        self.runNames = ['LICK2012','PAL2012','PAL2013']
        self.runButtons = [self.ui.lick2012Button, self.ui.pal2012Button, self.ui.pal2013Button]

        # Arrays with target names and buttons
        self.lick2012TargetNames = ['LICK2012_Object']
        self.lick2012TargetButtons = [self.ui.lick2012ObjectButton]
        self.pal2012TargetNames = ['SDSS_J0651', 'SDSS_J0926']
        self.pal2012TargetButtons = [self.ui.sdss0651Button,self.ui.sdss0926Button]
        self.pal2013TargetNames = ['PAL2013_Object']
        self.pal2013TargetButtons = [self.ui.pal2013ObjectButton]
        self.targetNames = [self.lick2012TargetNames, self.pal2012TargetNames, self.pal2013TargetNames]
        self.targetButtons = [self.lick2012TargetButtons, self.pal2012TargetButtons, self.pal2013TargetButtons]

        # Activate run and target buttons
        for iRun in range(len(self.runButtons)):
            self.runButtons[iRun].clicked.connect(self.chooseRun)
            for iTarget in range(len(self.targetButtons[iRun])):
                self.targetButtons[iRun][iTarget].clicked.connect(self.chooseTarget)

        self.ui.sunsetList.itemClicked.connect(self.createObsList)

    def testFunction(self):
        print 'test'

    # 2) Choose Run
    def chooseRun(self):
        for iRun in range(len(self.runButtons)):
            for iTarget in range(len(self.targetButtons[iRun])):
                    self.targetButtons[iRun][iTarget].setVisible(False)     
            if self.runButtons[iRun].isChecked():
                self.runNumber = iRun        
                for iTarget in range(len(self.targetButtons[iRun])):
                    self.targetButtons[iRun][iTarget].setVisible(True)  
        self.run = self.runNames[self.runNumber]      
        self.targetNumber = 0
        self.targetButtons[self.runNumber][self.targetNumber].setChecked(True)
        self.target = self.targetNames[self.runNumber][self.targetNumber]
        self.loadParamFile()
           
    # 3) Choose Target
    def chooseTarget(self):
        for iTarget in range(len(self.targetButtons[self.runNumber])):
            if self.targetButtons[self.runNumber][iTarget].isChecked():
                self.targetNumber = iTarget
        self.target = self.targetNames[self.runNumber][self.targetNumber]
        self.loadParamFile()

    # 4) Load Param File
    def loadParamFile(self):
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

    # 5) Choose Obs File
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

    # Create wavelength cal file list
    def createWavelengthList(self):
        self.ui.wavelengthList.clear()
        for iCal in range(len(self.calTimestamps)):
            self.ui.wavelengthList.addItem(self.calTimestamps[iCal])

    # Create flat cal file list
    def createFlatList(self):
        self.ui.flatList.clear()
        for iCal in range(len(self.flatCalDates)):
            self.ui.flatList.addItem(self.flatCalDates[iCal])

    # Enable/disable wavecal options
    def useWaveCal(self):
        if self.ui.wavelengthCalibrationBox.isChecked():
            self.ui.lowerWavelengthCutoffLine.setEnabled(True)
            self.ui.lowerWavelengthCutoffLabel.setEnabled(True)
            self.ui.upperWavelengthCutoffLine.setEnabled(True)
            self.ui.upperWavelengthCutoffLabel.setEnabled(True)
            self.ui.wavelengthList.setEnabled(True)
        else:
            self.ui.lowerWavelengthCutoffLine.setEnabled(False)
            self.ui.lowerWavelengthCutoffLabel.setEnabled(False)
            self.ui.upperWavelengthCutoffLine.setEnabled(False)
            self.ui.upperWavelengthCutoffLabel.setEnabled(False)
            self.ui.wavelengthList.setEnabled(False)

    # Enable/disable flatcal options
    def useFlatCal(self):
        if self.ui.flatCalibrationBox.isChecked():
            self.ui.deadPixelBox.setChecked(True)
            self.ui.deadPixelBox.setEnabled(True)
            self.ui.flatList.setEnabled(True)
        else:
            self.ui.deadPixelBox.setChecked(False)
            self.ui.deadPixelBox.setEnabled(False)
            self.ui.flatList.setEnabled(False)
        #self.printInformation()

    # Load dead pixel mask
    def loadDeadMask(self):
        self.deadPixelFilename = self.displayStackPath + 'deadPixelMasks/deadPixelMask_' + ts
        if not os.path.exists(self.deadPixelFilename + '.npz'):
            self.deadMask = self.ob.getDeadPixels()
            np.savez(self.deadPixelFilename,deadMask = self.deadMask)
            print "Dead pixel mask saved to %s"%(self.deadPixelFilename)
        else:
            self.deadFile = np.load(self.deadPixelFilename + '.npz')
            self.deadMask = self.deadFile['deadmask']


    # Load hot pixel mask
    def loadHotMask(self):
        self.hotPixelFilename = str(self.displayStackPath + self.run + '/' + self.target + '/HotPixelMasks/hotPixelMask_' + self.obsTS + '.h5')
        if not os.path.exists(self.hotPixelFilename):
            hp.findHotPixels(self.obsFn,self.hotPixelFilename)
            print "Hot pixel mask saved to %s"%(self.hotPixelFilename)
        self.ob.loadHotPixCalFile(self.hotPixelFilename,switchOnMask=True)

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
        self.lowerWavelengthCutoff = self.ui.lowerWavelengthCutoffLine.text()
        self.upperWavelengthCutoff = self.ui.upperWavelengthCutoffLine.text()
               
        # Flat calibration settings
        self.useFlatCalibration = self.ui.flatCalibrationBox.isChecked()
        self.useDeadPixelMasking = self.ui.deadPixelBox.isChecked()    

        self.fileCount = self.ui.inputList.count()  
            
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

    # Start process for creating image stacks
    def startProcess(self):

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
                        self.loadHotMask()
                    
                        

        # Invalid params file
        else:
            print 'Invalid parameter file...'
            

'''              

    # Later replace this with an information box on the gui, or probably just some log file
    def printInformation(self):
        print '................................................'
        print 'Run = ' + self.run
        print 'Target = ' + self.target
        print 'Sunset Date: ' + self.sunsetDate
        print 'UTC Date: ' + self.utcDate        
        print 'Sequence List: ' + str(self.sequenceList)
        print 'Integration Time = ' + str(self.integrationTime) + ' s'

        if self.useTimeAdjustment:
            print 'Time Adjustment: On'
        else:
            print 'Time Adjustment: Off'

        if self.useWavelengthCalibration:
            print 'Wavelength Calibration: On'
            print 'Lower Wavelength Cutoff: ' + self.lowerWavelengthCutoff
            print 'Upper Wavelength Cutoff: ' + self.upperWavelengthCutoff
        else:
            print 'Wavelength Calibration: Off'

        if self.useFlatCalibration:
            print 'Flat Calibration: On'
            if self.useDeadPixelMasking:
                print 'Dead Pixel Masking: On'
            else:
                print 'Dead Pixel Masking: Off'
        else:
            print 'Flat Calibration: Off'

        if self.useHotPixelMasking:
            print 'Hot Pixel Masking: On'
        else:
            print 'Hot Pixel Masking: Off'

        print 'Output filename = ' + self.outputFilename
    
    def startProcess(self):
        # Check to make sure an image stack can be made with the current settings      
        self.validateSettings()        
        if self.validSettings:
            # Initialize variables from settings data
            self.run = self.runNames[self.runNumber]
            self.target = self.targetNames[self.runNumber][self.targetNumber]
            self.sunsetDate = self.sunsetDates[self.sequenceNumber]
            self.utcDate = self.utcDates[self.sequenceNumber]
            self.sequenceList = self.sequences[self.sequenceNumber]
            self.integrationTime = int(self.ui.integrationTimeLine.text())
            self.useTimeAdjustment = self.ui.timeAdjustmentBox.isChecked()
            self.useWavelengthCalibration = self.ui.wavelengthCalibrationBox.isChecked()
            self.lowerWavelengthCutoff = self.ui.lowerWavelengthCutoffLine.text()
            self.upperWavelengthCutoff = self.ui.upperWavelengthCutoffLine.text()
            self.useFlatCalibration = self.ui.flatCalibrationBox.isChecked()
            self.useDeadPixelMasking = self.ui.deadPixelBox.isChecked()
            self.useHotPixelMasking = self.ui.hotPixelBox.isChecked()
            self.timestampList = [self.utcDate+'-'+str(ts) for ts in self.sequenceList]
            if self.useWavelengthCalibration:
                self.wvlCalFilename = self.ui.wavelengthCalibrationSolutionLine.text()
                self.outputFilename = self.displayStackPath + self.run + '/' + self.target + '/ImageStacks/'
            else:
                self.outputFilename = self.displayStackPath + self.run + '/' + self.target + '/ImageStacks/'
            if self.useFlatCalibration:
                self.flatCalFilename = self.ui.flatCalibrationSolutionLine.text()


            # Print information from selected settings and observation data
            self.printInformation()



            # Loop through all timestamps in a sequence
            for i,ts in enumerate(self.timestampList):
                print 'loading',ts
                # Create an ObsFile instance
                self.obsFn = FileName(run=self.run,date=self.sunsetDate,tstamp=ts).obs()
                self.ob = ObsFile(self.obsFn)

                # Load time adjustment file
                if self.useTimeAdjustment:
                    self.ob.loadTimeAdjustmentFile(FileName(run=self.run).timeAdjustments())
                
                # Load hot pixel mask
                if self.useHotPixelMasking:
                    self.hotPixFn = self.displayStackPath + 'HotPixelMasks/hotPixelMask' + ts + '.h5'
                    if not os.path.exists(self.hotPixFn):
                        hp.findHotPixels(self.obsFn,self.hotPixFn)
                        print "Flux file pixel mask saved to %s"%(self.hotPixFn)
                    self.ob.loadHotPixCalFile(self.hotPixFn,switchOnMask=True)

                # Load wavecal solution
                if self.useWavelengthCalibration:
                    self.ob.loadWvlCalFile(self.wvlCalFilename)
                    self.ob.setWvlCutoffs(self.lowerWavelengthCutoff,self.upperWavelengthCutoff)
                    print 'Loading wavelength calibration...'

                # Load flatcal solution
                if self.useFlatCalibration:
                    self.ob.loadFlatCalFile(self.flatCalFilename)
                    print 'Loading flat calibration...'
            
                # Create dead pixel mask
                if self.useDeadPixelMasking:
                    self.deadMask = self.ob.getDeadPixels()
                    print 'Creating dead pixel mask...'

                # Calculate start JD from more accurate unix time
                self.unix = self.ob.getFromHeader('unixtime')
                self.startJD = self.unix/86400.+2440587.5
                self.exptime = self.ob.getFromHeader('exptime')         

                #for sec in np.arange(0,nSecInFile,integrationTime):
'''        

# Start up main gui
if __name__ == "__main__":
	app = QApplication(sys.argv)
	myapp = DisplayStack()
	myapp.show()
	app.exec_()



