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

        self.loadStackData()

        self.displayPlot()

        self.ui.leftButton.clicked.connect(self.leftButtonClicked)
        self.ui.rightButton.clicked.connect(self.rightButtonClicked)
        
        self.ui.frameNumberLine.returnPressed.connect(self.numberEntered)

    def loadStackData(self):
        self.stackFile = np.load(str(self.loadStackName))
        self.stackData = np.array(self.stackFile['stack'])
        self.jdData = np.array(self.stackFile['jd'])
        self.totalFrames = self.stackData.shape[2]
        self.ui.maxLabel.setText('/ ' + str(self.totalFrames-1))

    def displayPlot(self):
        self.ui.frameNumberLine.setText(str(self.currentFrame))
        self.ui.jdLabel.setText('JD = ' + str(self.jdData[self.currentFrame]))
        self.selectedFrame = self.stackData[:,:,self.currentFrame]
        self.nanMask = np.isnan(self.selectedFrame)
        self.selectedFrame[self.nanMask] = 0.0
        self.ui.plotDock.canvas.ax.clear()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        self.ui.plotDock.canvas.ax.matshow(self.selectedFrame, cmap='gray', origin='lower')
        self.ui.plotDock.canvas.ax.xaxis.tick_bottom()
        self.ui.plotDock.canvas.draw()
        
    def rightButtonClicked(self):
        self.currentFrame+=1
        self.currentFrame%=self.totalFrames
        self.displayPlot()

    def leftButtonClicked(self):
        self.currentFrame-=1
        self.currentFrame%=self.totalFrames
        self.displayPlot()

    def numberEntered(self):
        self.currentFrame = int(self.ui.frameNumberLine.text())%self.totalFrames
        self.displayPlot()
