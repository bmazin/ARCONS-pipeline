#!/bin/python

from PyQt4.QtGui import *
from PyQt4.QtGui import *
from LoadImageStack_gui import Ui_LoadImageStack_gui
import DisplayStack
import numpy as np
import matplotlib.pyplot as plt
from GaussFitter import gaussfit
from util.readDict import readDict
import tables
import ephem
import PyGuide as pg

'''
Author: Paul Szypryt		Date: November 4, 2013
'''

class LoadImageStack(QDialog):

    def __init__(self, stackName, parent=None):
        QWidget.__init__(self, parent=None)
        self.ui = Ui_LoadImageStack_gui()
        self.ui.setupUi(self)
        
        self.loadStackName = str(stackName)
        print 'Loading ' + self.loadStackName
        self.currentFrame = 0

        # Initialize angle array for drawing circles.

        self.totalAngles = 100
        self.an = np.linspace(0,2*np.pi,self.totalAngles)
        
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
        
        self.ui.centroidButton.clicked.connect(self.performCentroiding)


        

    def arcsec_to_radians(self,total_arcsec):
        total_degrees = total_arcsec/3600.0
        total_radians = total_degrees*self.d2r
        return total_radians

    # Function for converting radians to arcseconds.
    def radians_to_arcsec(self,total_radians):
        total_degrees = total_radians*self.r2d
        total_arcsec = total_degrees*3600.0
        return total_arcsec

    def loadStackData(self):

        self.d2r = np.pi/180.0
        self.r2d = 180.0/np.pi

        self.stackFile = tables.openFile(self.loadStackName, mode='r')

        self.header = self.stackFile.root.header.header
        self.headerTitles = self.header.colnames
        self.headerInfo = self.header[0]
    
        self.stackNode = self.stackFile.root.stack

        self.stackData = np.array(self.stackNode.stack.read())
        self.jdData = np.array(self.stackNode.time.read()[0])

        self.nRow = self.stackData.shape[0]
        self.nCol = self.stackData.shape[1]
        self.totalFrames = self.stackData.shape[2]

        # Load data out of display stack h5 file
        self.centroid_RA = self.headerInfo[self.headerTitles.index('RA')]
        self.centroid_DEC = self.headerInfo[self.headerTitles.index('Dec')]
        self.original_lst = self.headerInfo[self.headerTitles.index('lst')]

        self.exptime = self.headerInfo[self.headerTitles.index('exptime')]
        self.integrationTime = self.headerInfo[self.headerTitles.index('integrationTime')]
        self.deadPixelFilename = self.headerInfo[self.headerTitles.index('deadPixFileName')]
        self.HA_offset = self.headerInfo[self.headerTitles.index('HA_offset')]
        self.nRow = self.headerInfo[self.headerTitles.index('nRow')]
        self.nCol = self.headerInfo[self.headerTitles.index('nCol')]


        self.centroid_RA_radians = ephem.hours(self.centroid_RA).real
        #self.centroid_RA_arcsec = self.radians_to_arcsec(self.centroid_RA_radians)/15.0
        self.centroid_RA_arcsec = self.radians_to_arcsec(self.centroid_RA_radians)
    
        self.centroid_DEC_radians = ephem.degrees(self.centroid_DEC).real
        self.centroid_DEC_arcsec = self.radians_to_arcsec(self.centroid_DEC_radians)
    
        self.original_lst_radians = ephem.hours(self.original_lst).real
        self.original_lst_seconds = self.radians_to_arcsec(self.original_lst_radians)/15.0
        #self.original_lst_seconds = self.radians_to_arcsec(self.original_lst_radians)

        self.HA_offset_radians = self.HA_offset*self.d2r
        
        self.ui.maxLabel.setText('/ ' + str(self.totalFrames-1))
        self.centerPositions = np.zeros((self.totalFrames,2))
        self.apertureRadii = np.zeros(self.totalFrames)
        self.annulusRadii = np.zeros((self.totalFrames,2))
        for iFrame in range(self.totalFrames):
            self.centerPositions[iFrame] = [float(self.nCol)/2.0 - 0.5, float(self.nRow)/2.0 - 0.5]
            self.apertureRadii[iFrame] = 5.0
            self.annulusRadii[iFrame] = [10.0, 20.0]
        self.updateFrameInfo()

    def displayPlot(self):
        self.ui.frameNumberLine.setText(str(self.currentFrame))
        self.ui.jdLabel.setText('JD = ' + str(self.jdData[self.currentFrame]))
        self.selectedFrame = np.array(self.stackData[:,:,self.currentFrame])
        self.nanMask = np.isnan(self.selectedFrame)
        self.selectedFrame[self.nanMask] = 0.0
        self.ui.plotDock.canvas.ax.clear()
        self.ui.plotDock.canvas.ax.matshow(self.selectedFrame, cmap='gray', origin='lower')
        self.ui.plotDock.canvas.ax.xaxis.tick_bottom()
        self.ui.plotDock.canvas.ax.set_xlim([-0.5,self.nCol-0.5])
        self.ui.plotDock.canvas.ax.set_ylim([-0.5,self.nRow-0.5])
        self.updateCircles()
        cid = self.ui.plotDock.canvas.mpl_connect('motion_notify_event', self.hoverCanvas)
        self.ui.plotDock.canvas.draw()
        self.drawCompass()

    def hoverCanvas(self,event):
        if event.inaxes is self.ui.plotDock.canvas.ax:
            col = int(round(event.xdata))
            row = int(round(event.ydata))
            if row < np.shape(self.selectedFrame)[0] and col < np.shape(self.selectedFrame)[1]:
                self.ui.matshowLabel.setText('({:d},{:d}) {}'.format(row,col,self.selectedFrame[row,col]))
                #print '({:d},{:d}) {}'.format(row,col,self.selectedFrame[row,col])
        
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
        if self.ui.annulusCheckbox.isChecked():
        # Plot annulus circles
            self.ui.plotDock.canvas.ax.plot( self.centerPositions[self.currentFrame][0]+self.annulusRadii[self.currentFrame][0]*np.cos(self.an), self.centerPositions[self.currentFrame][1]+self.annulusRadii[self.currentFrame][0]*np.sin(self.an), 'r')
            self.ui.plotDock.canvas.ax.plot( self.centerPositions[self.currentFrame][0]+self.annulusRadii[self.currentFrame][1]*np.cos(self.an), self.centerPositions[self.currentFrame][1]+self.annulusRadii[self.currentFrame][1]*np.sin(self.an), 'r')
        if self.ui.apertureCheckbox.isChecked():
        # Plot aperture circle
            self.ui.plotDock.canvas.ax.plot( self.centerPositions[self.currentFrame][0]+self.apertureRadii[self.currentFrame]*np.cos(self.an), self.centerPositions[self.currentFrame][1]+self.apertureRadii[self.currentFrame]*np.sin(self.an), 'b')

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

    def aperture(self,startpx,startpy,radius):
        r = radius
        length = 2*r 
        height = length
        allx = xrange(startpx-int(np.ceil(length/2.0)),startpx+int(np.floor(length/2.0))+1)
        ally = xrange(startpy-int(np.ceil(height/2.0)),startpy+int(np.floor(height/2.0))+1)
        mask=np.zeros((46,44))
        
        for x in allx:
            for y in ally:
                if (np.abs(x-startpx))**2+(np.abs(y-startpy))**2 <= (r)**2 and 0 <= y and y < 46 and 0 <= x and x < 44:
                    mask[y,x]=1.
        return mask
    
    def drawCompass(self):
        compassAngle = np.zeros(self.totalFrames)
        currentLST = (self.original_lst_seconds + self.integrationTime/2.0 + self.currentFrame*self.integrationTime)*15.0 #in arcsec
        variableHA = currentLST - self.centroid_RA_arcsec
        variableHA_radians = self.arcsec_to_radians(variableHA)
        currentHA_radians = self.HA_offset_radians + variableHA_radians
        #print currentHA_radians

        currentHA_degrees = currentHA_radians*self.r2d

        northLine = np.linspace(0.0,1.0,2)
        northLineX = northLine*np.cos(np.pi/2.0-currentHA_radians)
        northLineY = northLine*np.sin(np.pi/2.0-currentHA_radians)

        arrowLine = np.linspace(0.0,0.3,2)

        northArrowX = northLineX[-1]
        northArrowY = northLineY[-1]

        arrowLineNX1 = northArrowX - arrowLine*np.cos(np.pi/4.-currentHA_radians)
        arrowLineNY1 = northArrowY - arrowLine*np.sin(np.pi/4.-currentHA_radians)
        arrowLineNX2 = northArrowX - arrowLine*np.cos(3.*np.pi/4.-currentHA_radians)
        arrowLineNY2 = northArrowY - arrowLine*np.sin(3.*np.pi/4.-currentHA_radians)

        eastLineX = northLine*np.cos(np.pi-currentHA_radians)
        eastLineY = northLine*np.sin(np.pi-currentHA_radians)
        
        eastArrowX = eastLineX[-1]
        eastArrowY = eastLineY[-1]

        arrowLineEX1 = eastArrowX - arrowLine*np.cos(3.*np.pi/4.-currentHA_radians)
        arrowLineEY1 = eastArrowY - arrowLine*np.sin(3.*np.pi/4.-currentHA_radians)
        arrowLineEX2 = eastArrowX - arrowLine*np.cos(5.*np.pi/4.-currentHA_radians)
        arrowLineEY2 = eastArrowY - arrowLine*np.sin(5.*np.pi/4.-currentHA_radians)
        

        self.ui.hourAngleLabel.setText(str(currentHA_degrees))

        stretchFactor = 1.8

        self.ui.compassDock.canvas.ax.clear()
        self.ui.compassDock.canvas.ax.plot(northLineX,northLineY, 'k', arrowLineNX1, arrowLineNY1, 'k', arrowLineNX2, arrowLineNY2, 'k', eastLineX,eastLineY, 'k', arrowLineEX1, arrowLineEY1, 'k', arrowLineEX2, arrowLineEY2, 'k', [0.0,0.0], [-stretchFactor,stretchFactor], 'b--', [-stretchFactor,stretchFactor], [0.0,0.0], 'b--')
        self.ui.compassDock.canvas.ax.annotate('N', (northArrowX*1.4,northArrowY*1.4), horizontalalignment='center', verticalalignment='center')
        self.ui.compassDock.canvas.ax.annotate('E', (eastArrowX*1.4,eastArrowY*1.4), horizontalalignment='center', verticalalignment='center')
        self.ui.compassDock.canvas.ax.set_xlim([-stretchFactor,stretchFactor])
        self.ui.compassDock.canvas.ax.set_ylim([-stretchFactor,stretchFactor])
        self.ui.compassDock.canvas.ax.get_xaxis().set_visible(False)
        self.ui.compassDock.canvas.ax.get_yaxis().set_visible(False)
        self.ui.compassDock.canvas.draw()
        


    def performAperturePhotometry(self):
        
        self.doSkySubtraction = self.ui.skySubtractionCheckbox.isChecked()
        if self.doSkySubtraction:
            print 'Performing aperture photometry with sky subtraction...'
        else:
            print 'Performing aperture photometry without sky subtraction...'
        # Create file name
        self.objectIdentifier = self.loadStackName[0:self.loadStackName.index('ImageStacks/ImageStack_')]
        self.obsIdentifier = self.loadStackName[self.loadStackName.rindex('ImageStacks/ImageStack_') + len('ImageStacks/ImageStack_'):-3]
        self.outputFileName = self.objectIdentifier + 'ApertureStacks/ApertureStack_' + self.obsIdentifier + '.npz'
        print 'Saving to ' + self.outputFileName

        self.centerPositions[self.currentFrame] = [float(self.ui.xLine.text()),float(self.ui.yLine.text())]
        self.apertureRadii[self.currentFrame] = float(self.ui.apertureLine.text())
        self.annulusRadii[self.currentFrame] = [float(self.ui.innerAnnulusLine.text()),float(self.ui.outerAnnulusLine.text())]

        if self.doSkySubtraction:
            # Cycle through all frames, performing aperture photometry on each individually.
            frameCounts=[]
            for iFrame in range(self.totalFrames):
                
                startpx = int(np.round(self.centerPositions[iFrame][0],0))
                startpy = int(np.round(self.centerPositions[iFrame][1],0))

                apertureRadius = self.apertureRadii[iFrame]
                innerRadius = self.annulusRadii[iFrame][0]
                outerRadius = self.annulusRadii[iFrame][1]
                
                apertureMask = self.aperture(startpx, startpy, apertureRadius)

                innerMask = self.aperture(startpx, startpy, innerRadius)
                outerMask = self.aperture(startpx, startpy, outerRadius)
                annulusMask = outerMask-innerMask

                currentImage = np.array(self.stackData[:,:,iFrame])
                

                currentImage[np.isnan(currentImage)] = 0.0#set to finite value that will be ignored
                nanMask = currentImage == 0.0
                #nanMask = np.isnan(currentImage)

                aperturePixels = np.array(np.where(np.logical_and(apertureMask==1, nanMask==False)))
                aperturePix = aperturePixels.shape[1]
                apertureCountsPerPixel = np.sum(currentImage[aperturePixels[0],aperturePixels[1]])/aperturePix

                annulusPixels = np.array(np.where(np.logical_and(annulusMask==1, nanMask==False)))
                annulusPix = annulusPixels.shape[1]
                annulusCountsPerPixel = np.sum(currentImage[annulusPixels[0],annulusPixels[1]])/annulusPix

                frameCounts.append((apertureCountsPerPixel - annulusCountsPerPixel) * aperturePix)

            np.savez(self.outputFileName, counts=frameCounts, jd=self.jdData)

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.jdData, frameCounts)
            ax.set_xlabel('JD')
            ax.set_ylabel('Counts')
            plt.show()
                             

        else:
            # Cycle through all frames, performing aperture photometry on each individually.
            frameCounts=[]
            integrationTime = self.headerInfo[self.headerTitles.index('integrationTime')]
            for iFrame in range(self.totalFrames):
                apertureRadius = self.apertureRadii[iFrame]
                startpx = int(np.round(self.centerPositions[iFrame][0],0))
                startpy = int(np.round(self.centerPositions[iFrame][1],0))
                
                apertureMask = self.aperture(startpx, startpy, apertureRadius)
                
                #aperturePixels = np.where(apertureMask==1)
                
                currentImage = np.array(self.stackData[:,:,iFrame])

                #print currentImage

                #currentImage[np.where(currentImage==0.0)] = np.nan
                currentImage[np.isnan(currentImage)] = 0.0#set to finite value that will be ignored
                nanMask = currentImage == 0.0
                #nanMask = np.isnan(currentImage)
                #nanMask = currentImage == 0.0

                #currentImage[nanMask] = 0.0

                aperturePixels = np.array(np.where(np.logical_and(apertureMask==1, nanMask==False)))

                aperturePix = aperturePixels.shape[1]
                print aperturePix
                print integrationTime

                #print aperturePixels
                #print currentImage[aperturePixels]

                #print currentImage[aperturePixels]

                apertureCountsPerSecondPerPixel = np.sum(currentImage[aperturePixels[0], aperturePixels[1]]) / aperturePix

                print currentImage[aperturePixels[0],aperturePixels[1]] / aperturePix
                #print currentImagePerSecond


                #print currentImagePerSecond[aperturePixels]

                

                #apertureCountsPerSecondPerPixel = np.sum(currentImagePerSecond)/aperturePix

                #currentImage[nanMask] = 0.0
                frameCounts.append(apertureCountsPerSecondPerPixel)
            np.savez(self.outputFileName, counts=frameCounts, jd=self.jdData)

            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(self.jdData, frameCounts)
            ax.set_xlabel('JD')
            ax.set_ylabel('Counts')
            plt.show()
        
        print 'Done performing aperture photometry...'

    def performPSFPhotometry(self):
        print 'Performing PSF fitting photometry...'
        # Create file name
        self.objectIdentifier = self.loadStackName[0:self.loadStackName.index('ImageStacks/ImageStack_')]
        self.obsIdentifier = self.loadStackName[self.loadStackName.rindex('ImageStacks/ImageStack_') + len('ImageStacks/ImageStack_'):-3]
        self.outputFileName = self.objectIdentifier + 'FittedStacks/FittedStack_' + self.obsIdentifier + '.npz'
        print 'Saving to ' + self.outputFileName

        self.centerPositions[self.currentFrame] = [float(self.ui.xLine.text()),float(self.ui.yLine.text())]
        self.apertureRadii[self.currentFrame] = float(self.ui.apertureLine.text())
        self.annulusRadii[self.currentFrame] = [float(self.ui.innerAnnulusLine.text()),float(self.ui.outerAnnulusLine.text())]
        

        paramsList = []
        errorsList = []
        fitImgList = []
        chisqList = []

        for iFrame in range(self.totalFrames):
            guessX = int(np.round(self.centerPositions[iFrame][0],0))
            guessY = int(np.round(self.centerPositions[iFrame][1],0))
            apertureRadius = self.apertureRadii[iFrame]

            apertureMask = self.aperture(guessX, guessY, apertureRadius)

            currentImage = np.array(self.stackData[:,:,iFrame])
    
            #nanMask = np.isnan(currentImage)
            currentImage[np.isnan(currentImage)] = 0.0#set to finite value that will be ignored
            nanMask = currentImage == 0.0
            print nanMask

            err = np.sqrt(currentImage)
            #err = np.ones(np.shape(currentImage))
            err[apertureMask==0] = np.inf#weight points closer to the expected psf higher
            #currentImage[nanMask]=0#set to finite value that will be ignored
            err[nanMask] = np.inf#ignore these data points
            nearDeadCutoff=1#100/15 cps for 4000-6000 angstroms
            err[currentImage<nearDeadCutoff] = np.inf
            entireMask = (err==np.inf)
            maFrame = np.ma.masked_array(currentImage,entireMask)
            guessAmp = 30.0
            guessHeight = 30.0
            guessWidth=1.5
            guessParams = [guessHeight,guessAmp,guessX,guessY,guessWidth]
            limitedmin = 5*[True] 
            limitedmax = 5*[True]
            minpars = [0,0,0,0,.1]
            maxpars = [100.0,100.0,43,45,3.0]
            usemoments=[True,True,True,True,True] #doesn't use our guess values

            out = gaussfit(data=maFrame,err=err,params=guessParams,returnfitimage=True,quiet=True,limitedmin=limitedmin,limitedmax=limitedmax,minpars=minpars,maxpars=maxpars,circle=1,usemoments=usemoments,returnmp=True)
            mp = out[0]

            outparams = mp.params
            paramErrors = mp.perror
            chisq = mp.fnorm
            dof = mp.dof
            reducedChisq = chisq/dof
            print reducedChisq
            fitimg = out[1]
            chisqList.append([chisq,dof])

            paramsList.append(outparams)
            errorsList.append(paramErrors)
            print outparams,paramErrors

            fitimg[nanMask]=0           
            fitImgList.append(fitimg)
            currentImage[nanMask]=np.nan

        cube = np.array(fitImgList)
        chisqs = np.array(chisqList)
        params = np.array(paramsList)
        errors = np.array(errorsList)

        np.savez(self.outputFileName,fitImg=cube,params=params,errors=errors,chisqs=chisqs,jd=self.jdData)

        amps = params[:,1]
        widths = params[:,4]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.jdData, amps*widths**2)
        ax.set_xlabel('JD')
        ax.set_ylabel('Power')
        plt.show()    

        print 'Done performing PSF fitting photometry...'

    def performCentroiding(self):
        '''
        # Function for converting arcseconds to radians.
        def arcsec_to_radians(total_arcsec):
            total_degrees = total_arcsec/3600.0
            total_radians = total_degrees*d2r
            return total_radians

        # Function for converting radians to arcseconds.
        def radians_to_arcsec(total_radians):
            total_degrees = total_radians*r2d
            total_arcsec = total_degrees*3600.0
            return total_arcsec

        print 'Not currently implemented...'

        d2r = np.pi/180.0
        r2d = 180.0/np.pi

        # Load data out of display stack h5 file
        centroid_RA = self.headerInfo[self.headerTitles.index('RA')]
        centroid_DEC = self.headerInfo[self.headerTitles.index('Dec')]
        original_lst = self.headerInfo[self.headerTitles.index('lst')]
        exptime = self.headerInfo[self.headerTitles.index('exptime')]
        integrationTime = self.headerInfo[self.headerTitles.index('integrationTime')]
        deadPixelFilename = self.headerInfo[self.headerTitles.index('deadPixFileName')]
        HA_offset = self.headerInfo[self.headerTitles.index('HA_offset')]
        nRow = self.headerInfo[self.headerTitles.index('nRow')]
        nCol = self.headerInfo[self.headerTitles.index('nCol')]

        centroid_RA_radians = ephem.hours(centroid_RA).real
        centroid_RA_arcsec = radians_to_arcsec(centroid_RA_radians)
    
        centroid_DEC_radians = ephem.degrees(centroid_DEC).real
        centroid_DEC_arcsec = radians_to_arcsec(centroid_DEC_radians)
    
        original_lst_radians = ephem.hours(original_lst).real
        original_lst_seconds = radians_to_arcsec(original_lst_radians)/15.0

        '''

        # Load up dead pixel mask.  Invert for PyGuide format.
        #deadFile = np.load(self.deadPixelFilename)
        #deadMask = deadFile['deadMask']
        #deadMask = -1*deadMask + 1
        
        # Saturated pixels already taken care of by hot pixel code.
        satMask = np.zeros((46,44))

        # Specify CCDInfo (bias,readNoise,ccdGain,satLevel)
        ccd = pg.CCDInfo(0,0.00001,1,2500)

        # Initialize arrays that will be saved in h5 file. 1 array element per centroid frame.
        timeList=[]
        xPositionList=[]
        yPositionList=[]
        hourAngleList=[]
        flagList=[]
    
        flag=0
    
        print 'Centroiding a total of ' + str(self.totalFrames) + ' frames...'

        for iFrame in range(self.totalFrames):
            apertureRadius = self.apertureRadii[iFrame]
            image = np.array(self.stackData[:,:,iFrame])
            nanMask = np.isnan(image)
            image[nanMask] = 0.
        
            deadMask = np.zeros((self.nRow,self.nCol))
            deadMask[np.where(image==0)] = 1

            xyguess = np.array(self.centerPositions[iFrame])
            pyguide_output = pg.centroid(image,deadMask,satMask,xyguess,apertureRadius,ccd,0,False,verbosity=2, doDS9=True)  
             # Use PyGuide centroid positions, if algorithm failed, use xy guess center positions instead
            print pyguide_output
            try:
                xycenter = [float(pyguide_output.xyCtr[0]),float(pyguide_output.xyCtr[1])]
                print 'Frame ' + str(iFrame) +': Calculated [x,y] center = ' + str((xycenter)) + '.'
                flag = 0
            except TypeError:
                print 'Cannot centroid frame' + str(iFrame) + ', using guess instead'
                xycenter = xyguess
                flag = 1


            # Calculate lst for a given frame, at midpoint of frame
            current_lst_seconds = self.original_lst_seconds + (iFrame+0.5)*self.integrationTime
            current_lst_radians = self.arcsec_to_radians(current_lst_seconds*15.0)
            # Calculate hour angle for a given frame. Include a constant offset for instrumental rotation.
            HA_variable = current_lst_radians - self.centroid_RA_radians
            HA_static = self.HA_offset*self.d2r
            HA_current = HA_variable + HA_static
            # Make lists to save to h5 file
            timeList.append((iFrame+0.5)*self.integrationTime)
            xPositionList.append(xycenter[0])
            yPositionList.append(xycenter[1])
            hourAngleList.append(HA_current)
            flagList.append(flag)

        self.objectIdentifier = self.loadStackName[0:self.loadStackName.index('ImageStacks/ImageStack_')]
        self.obsIdentifier = self.loadStackName[self.loadStackName.rindex('ImageStacks/ImageStack_') + len('ImageStacks/ImageStack_'):-3]
        fullCentroidListFileName = self.objectIdentifier + 'CentroidLists/Centroid_' + self.obsIdentifier + '.h5'
        print 'Saving to ' + fullCentroidListFileName

        # Write data to new h5 file
        centroidListFile = tables.openFile(fullCentroidListFileName,mode='w')

        centroidHeaderGroupName = 'header'
        centroidHeaderTableName = 'header'
        centroidDataGroupName = 'centroidlist'


        centroidDataGroup = centroidListFile.createGroup(centroidListFile.root,centroidDataGroupName,'Table of times, x positions, y positions, hour angles, and flags')
    

        #paramstable = tables.Array(centroidgroup,'params', object=paramsList, title = 'Object and array params')
        timestable = tables.Array(centroidDataGroup,'times',object=timeList,title='Times at which centroids were calculated')
        xpostable = tables.Array(centroidDataGroup,'xPositions',object=xPositionList,title='X centroid positions')
        ypostable = tables.Array(centroidDataGroup,'yPositions',object=yPositionList,title='Y centroid positions')
        hatable = tables.Array(centroidDataGroup,'hourAngles',object=hourAngleList,title='Hour angles at specified times')
        flagtable = tables.Array(centroidDataGroup,'flags',object=flagList,title='Flags whether or not guess had to be used, 1 for guess used')


        # Should probably just copy original header information over, more info this way!!!
        centroidHeaderGroup = centroidListFile.createGroup("/", centroidHeaderGroupName, 'Header')
        centroidHeaderTable = centroidListFile.createTable(centroidHeaderGroup, centroidHeaderTableName, DisplayStack.DisplayStack.headerDescription,
                                            'Header Info')


        centroidHeader = centroidHeaderTable.row

        for iItem in range(len(self.headerInfo)):
            centroidHeader[self.headerTitles[iItem]] = self.headerInfo[iItem]

        centroidHeader.append()

        centroidListFile.flush()
        centroidListFile.close()
  
        print 'Done performing centroiding...'






