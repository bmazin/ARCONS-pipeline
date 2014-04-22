#!/bin/python

from PyQt4.QtGui import *
from PyQt4.QtGui import *
from LoadImageStack_gui import Ui_LoadImageStack_gui
import DisplayStack
import numpy as np
import matplotlib.pyplot as plt

'''
Author: Paul Szypryt		Date: November 4, 2013
'''

class LoadImageStack(QDialog):

    def __init__(self, stackName, parent=None):
        QWidget.__init__(self, parent=None)
        self.ui = Ui_LoadImageStack_gui()
        self.ui.setupUi(self)
        
        self.loadStackName = str(stackName)
        self.currentFrame = 0

        # Initialize angle array for drawing circles.
        self.an = np.linspace(0,2*np.pi,100)
        
        self.loadStackData()

        self.displayPlot()

        self.ui.leftButton.clicked.connect(self.leftButtonClicked)
        self.ui.rightButton.clicked.connect(self.rightButtonClicked)
        
        self.ui.frameNumberLine.returnPressed.connect(self.numberEntered)
        
        self.ui.xLine.returnPressed.connect(self.displayPlot)
        self.ui.yLine.returnPressed.connect(self.displayPlot)
        self.ui.apertureLine.returnPressed.connect(self.displayPlot)
        self.ui.innerAnnulusLine.returnPressed.connect(self.displayPlot)
        self.ui.outerAnnulusLine.returnPressed.connect(self.displayPlot)

        self.ui.allButton.clicked.connect(self.allButtonClicked)
        self.ui.subsequentButton.clicked.connect(self.subsequentButtonClicked)
        self.ui.previousButton.clicked.connect(self.previousButtonClicked)

        self.ui.apertureCheckbox.clicked.connect(self.displayPlot)
        self.ui.annulusCheckbox.clicked.connect(self.displayPlot)

        self.ui.psfPhotometryButton.clicked.connect(self.performPSFPhotometry)
        self.ui.aperturePhotometryButton.clicked.connect(self.performAperturePhotometry)

    def loadStackData(self):
        self.stackFile = np.load(str(self.loadStackName))
        self.stackData = np.array(self.stackFile['stack'])
        self.jdData = np.array(self.stackFile['jd'])
        self.nRow = self.stackData.shape[0]
        self.nCol = self.stackData.shape[1]
        self.totalFrames = self.stackData.shape[2]
        
        self.ui.maxLabel.setText('/ ' + str(self.totalFrames-1))
        self.centerPositions = np.zeros((self.totalFrames,2))
        self.apertureRadii = np.zeros(self.totalFrames)
        self.annulusRadii = np.zeros((self.totalFrames,2))
        for iFrame in range(self.totalFrames):
            self.centerPositions[iFrame] = [float(self.nCol)/2.0 - 0.5, float(self.nRow)/2.0 - 0.5]
            self.apertureRadii[iFrame] = 5.0
            self.annulusRadii[iFrame] = [5.0,15.0]
        self.updateFrameInfo()

    def displayPlot(self):
        self.ui.frameNumberLine.setText(str(self.currentFrame))
        self.ui.jdLabel.setText('JD = ' + str(self.jdData[self.currentFrame]))
        self.selectedFrame = self.stackData[:,:,self.currentFrame]
        self.nanMask = np.isnan(self.selectedFrame)
        self.selectedFrame[self.nanMask] = 0.0
        self.ui.plotDock.canvas.ax.clear()
        self.fig = plt.figure()
        ax = self.fig.add_subplot(111)
        self.ui.plotDock.canvas.ax.matshow(self.selectedFrame, cmap='gray', origin='lower')
        self.ui.plotDock.canvas.ax.xaxis.tick_bottom()
        self.ui.plotDock.canvas.ax.set_xlim([-0.5,self.nCol-0.5])
        self.ui.plotDock.canvas.ax.set_ylim([-0.5,self.nRow-0.5])
        self.updateCircles()
        self.ui.plotDock.canvas.draw()
        
    def rightButtonClicked(self):
        self.currentFrame+=1
        self.currentFrame%=self.totalFrames
        self.updateFrameInfo()     
        self.displayPlot()

    def leftButtonClicked(self):
        self.currentFrame-=1
        self.currentFrame%=self.totalFrames
        self.updateFrameInfo() 
        self.displayPlot()

    def numberEntered(self):
        self.currentFrame = int(self.ui.frameNumberLine.text())%self.totalFrames
        self.updateFrameInfo() 
        self.displayPlot()

    def updateCircles(self):
        self.centerPositions[self.currentFrame] = [float(self.ui.xLine.text()),float(self.ui.yLine.text())]
        self.apertureRadii[self.currentFrame] = float(self.ui.apertureLine.text())
        self.annulusRadii[self.currentFrame] = [float(self.ui.innerAnnulusLine.text()),float(self.ui.outerAnnulusLine.text())]
        if self.ui.apertureCheckbox.isChecked():
        # Plot aperture circle
            self.ui.plotDock.canvas.ax.plot( self.centerPositions[self.currentFrame][0]+self.apertureRadii[self.currentFrame]*np.cos(self.an), self.centerPositions[self.currentFrame][1]+self.apertureRadii[self.currentFrame]*np.sin(self.an), 'b')
        if self.ui.annulusCheckbox.isChecked():
        # Plot annulus circles
            self.ui.plotDock.canvas.ax.plot( self.centerPositions[self.currentFrame][0]+self.annulusRadii[self.currentFrame][0]*np.cos(self.an), self.centerPositions[self.currentFrame][1]+self.annulusRadii[self.currentFrame][0]*np.sin(self.an), 'r')
            self.ui.plotDock.canvas.ax.plot( self.centerPositions[self.currentFrame][0]+self.annulusRadii[self.currentFrame][1]*np.cos(self.an), self.centerPositions[self.currentFrame][1]+self.annulusRadii[self.currentFrame][1]*np.sin(self.an), 'r')

    def updateFrameInfo(self):
        self.ui.apertureLine.setText(str(self.apertureRadii[self.currentFrame]))
        self.ui.innerAnnulusLine.setText(str(self.annulusRadii[self.currentFrame][0]))
        self.ui.outerAnnulusLine.setText(str(self.annulusRadii[self.currentFrame][1]))
        self.ui.xLine.setText(str(self.centerPositions[self.currentFrame][0]))
        self.ui.yLine.setText(str(self.centerPositions[self.currentFrame][1]))

    def allButtonClicked(self):
        for iFrame in range(self.totalFrames):
            self.centerPositions[iFrame] = self.centerPositions[self.currentFrame]
            self.apertureRadii[iFrame] = self.apertureRadii[self.currentFrame]
            self.annulusRadii[iFrame] = self.annulusRadii[self.currentFrame]

    def subsequentButtonClicked(self):
        for iFrame in range(self.currentFrame,self.totalFrames):
            self.centerPositions[iFrame] = self.centerPositions[self.currentFrame]
            self.apertureRadii[iFrame] = self.apertureRadii[self.currentFrame]
            self.annulusRadii[iFrame] = self.annulusRadii[self.currentFrame]

    def previousButtonClicked(self):
        for iFrame in range(0,self.currentFrame):
            self.centerPositions[iFrame] = self.centerPositions[self.currentFrame]
            self.apertureRadii[iFrame] = self.apertureRadii[self.currentFrame]
            self.annulusRadii[iFrame] = self.annulusRadii[self.currentFrame]

    def performAperturePhotometry(self):
        self.doSkySubtraction = self.ui.skySubtractionCheckbox.isChecked()
        if self.doSkySubtraction:
            print 'Performing aperture photometry with sky subtraction...'
        else:
            print 'Performing aperture photometry without sky subtraction...'
        self.centerPositions[self.currentFrame] = [float(self.ui.xLine.text()),float(self.ui.yLine.text())]
        self.apertureRadii[self.currentFrame] = float(self.ui.apertureLine.text())
        self.annulusRadii[self.currentFrame] = [float(self.ui.innerAnnulusLine.text()),float(self.ui.outerAnnulusLine.text())]
        print 'Done performing aperture photometry...'

    def performPSFPhotometry(self):
        print 'Performing PSF fitting photometry...'
        self.centerPositions[self.currentFrame] = [float(self.ui.xLine.text()),float(self.ui.yLine.text())]
        self.apertureRadii[self.currentFrame] = float(self.ui.apertureLine.text())
        self.annulusRadii[self.currentFrame] = [float(self.ui.innerAnnulusLine.text()),float(self.ui.outerAnnulusLine.text())]
        print 'Done performing PSF fitting photometry...'
        


