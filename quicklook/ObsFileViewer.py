'''
Author: Matt Strader        Date: January 06, 2015
'''

import sys, os
import numpy as np
from PyQt4 import QtCore
from PyQt4 import QtGui
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib
from functools import partial

from util.FileName import FileName
from util.ObsFile import ObsFile
from util.CalLookupFile import CalLookupFile
from util import utils
from hotpix.manuallyRemovePixels import removePixel

class ObsFileViewer(QtGui.QMainWindow):
    def __init__(self, obsPath=None, obsTstamp=None):
        self.app = QtGui.QApplication([])
        self.app.setStyle('plastique')
        super(ObsFileViewer,self).__init__()
        self.setWindowTitle('ObsFile Viewer')
        self.createWidgets()
        self.initWindows()
        self.createMenu()
        self.createStatusBar()
        self.arrangeMainFrame()
        self.connectControls()

        bInitWithObs = (not obsPath is None) or (not obsTstamp is None)

        if not obsPath is None:
            self.loadObsFile(obsPath)
        elif not obsTstamp is None:
            obsPath = CalLookupFile().obs(obsTstamp)
            self.loadObsFile(obsPath)

        if bInitWithObs:
            self.obs.loadAllCals()
            if not self.obs.hotPixTimeMask is None:
                self.enableTimeMaskParams()

    def show(self):
        super(ObsFileViewer,self).show()
        self.app.exec_()

    def initWindows(self):
        self.headerWindow = HeaderWindow(self)
        self.imageParamsWindow = ImageParamsWindow(self)
        self.plotWindows = []

    def newPlotWindow(self):
        newPlotId = len(self.plotWindows)
        plotWindow = PlotWindow(parent=self,plotId=newPlotId,selectedPixels=self.arrayImageWidget.selectedPixels)
        self.plotWindows.append(plotWindow)
        self.connect(self.arrayImageWidget,QtCore.SIGNAL('newPixelSelection(PyQt_PyObject)'), plotWindow.newPixelSelection)
        
        plotWindow.show()

    def createWidgets(self):
        self.mainFrame = QtGui.QWidget()
        
        self.arrayImageWidget = ArrayImageWidget(parent=self,hoverCall=self.hoverCanvas)

        # Create the navigation toolbar, tied to the canvas
        #self.canvasToolbar = NavigationToolbar(self.canvas, self.mainFrame)

        self.label_obsFilename = QtGui.QLabel('No obs loaded')

        self.textbox_startTime = QtGui.QLineEdit('0')
        self.textbox_startTime.setMaximumWidth(50)
        self.textbox_endTime = QtGui.QLineEdit('0')
        self.textbox_endTime.setMaximumWidth(50)
        self.label_unitsAfterEndTime = QtGui.QLabel('s')

        self.textbox_startWvl = QtGui.QLineEdit('3000')
        self.textbox_endWvl = QtGui.QLineEdit('11000')
        #self.textbox_endTime.setMaximumWidth(50)

        self.button_drawPlot = QtGui.QPushButton('Plot')


        self.textbox_timeIncrement = QtGui.QLineEdit('10')
        #self.textbox_timeIncrement.setMaximumWidth(50)
        self.button_jumpToBeginning = QtGui.QPushButton('|<')
        self.button_jumpToEnd = QtGui.QPushButton('>|')
        self.button_incrementBack = QtGui.QPushButton('<')
        self.button_incrementForward = QtGui.QPushButton('>')

        self.textbox_wvlIncrement = QtGui.QLineEdit('1000')
        #self.textbox_timeIncrement.setMaximumWidth(50)
        self.button_jumpToStartWvl = QtGui.QPushButton('|<')
        self.button_jumpToEndWvl = QtGui.QPushButton('>|')
        self.button_incrementWvlBack = QtGui.QPushButton('<')
        self.button_incrementWvlForward = QtGui.QPushButton('>')

        self.button_jumpToBeginning.setMaximumWidth(30)
        self.button_jumpToEnd.setMaximumWidth(30)
        self.button_incrementBack.setMaximumWidth(30)
        self.button_incrementForward.setMaximumWidth(30)

        self.button_jumpToStartWvl.setMaximumWidth(30)
        self.button_jumpToEndWvl.setMaximumWidth(30)
        self.button_incrementWvlBack.setMaximumWidth(30)
        self.button_incrementWvlForward.setMaximumWidth(30)


    def arrangeMainFrame(self):

        canvasBox = layoutBox('V',[self.arrayImageWidget])

        timeBoundsBox = layoutBox('H',[3,self.textbox_startTime,'s to',self.textbox_endTime,'s',1,self.button_drawPlot,3])
        incrementControlsBox = layoutBox('H',[1,self.button_jumpToBeginning,self.button_incrementBack,
                    self.textbox_timeIncrement,'s',self.button_incrementForward,self.button_jumpToEnd,1])
        wvlBoundsBox = layoutBox('H',[3,self.textbox_startWvl,'Angstroms to',self.textbox_endWvl,'Angstroms',1,3])
        incrementWvlControlsBox = layoutBox('H',[1,self.button_jumpToStartWvl,self.button_incrementWvlBack,
                    self.textbox_wvlIncrement,'Angstroms',self.button_incrementWvlForward,self.button_jumpToEndWvl,1])
        mainBox = layoutBox('V',[self.label_obsFilename,canvasBox,timeBoundsBox,incrementControlsBox, wvlBoundsBox, incrementWvlControlsBox])

        self.mainFrame.setLayout(mainBox)
        self.setCentralWidget(self.mainFrame)
  
    def createStatusBar(self):
        self.statusText = QtGui.QLabel("Click File->Load Obs File")
        self.statusBar().addWidget(self.statusText, 1)
        
    def createMenu(self):        
        self.fileMenu = self.menuBar().addMenu("&File")
        
        loadObsAction = createAction(self,"Load obs file...",slot=self.loadObsFileWindow)
        loadCalsAction = createAction(self,"Load cal files...",slot=self.loadCalFiles)
        
        quitAction = createAction(self,"&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        addActions(self.fileMenu, 
            (loadObsAction, loadCalsAction, None, quitAction))

        self.windowsMenu = self.menuBar().addMenu("&Windows")
        headerAction = createAction(self,"Header Info",slot=self.headerWindow.show)
        imgParamsAction = createAction(self,"Image Plot Parameters",slot=self.imageParamsWindow.show)
        newPlotWindowAction = createAction(self,'New Plot Window',slot=self.newPlotWindow)

        addActions(self.windowsMenu,(newPlotWindowAction,None,headerAction,imgParamsAction))
        
        self.helpMenu = self.menuBar().addMenu("&Help")
        aboutAction = createAction(self,"&About", 
            shortcut='F1', slot=self.aboutMessage)
        
        addActions(self.helpMenu, (aboutAction,))

    def connectControls(self):
        self.connect(self.button_drawPlot,QtCore.SIGNAL('clicked()'), self.getObsImage)
        self.connect(self.button_jumpToBeginning,QtCore.SIGNAL('clicked()'), self.jumpToBeginning)
        self.connect(self.button_jumpToEnd,QtCore.SIGNAL('clicked()'), self.jumpToEnd)
        self.connect(self.button_incrementForward,QtCore.SIGNAL('clicked()'), self.incrementForward)
        self.connect(self.button_incrementBack,QtCore.SIGNAL('clicked()'), self.incrementBack)
        self.connect(self.button_jumpToStartWvl,QtCore.SIGNAL('clicked()'), self.jumpToStartWvl)
        self.connect(self.button_jumpToEndWvl,QtCore.SIGNAL('clicked()'), self.jumpToEndWvl)
        self.connect(self.button_incrementWvlForward,QtCore.SIGNAL('clicked()'), self.incrementWvlForward)
        self.connect(self.button_incrementWvlBack,QtCore.SIGNAL('clicked()'), self.incrementWvlBack)

    def addClickFunc(self,clickFunc):
        self.arrayImageWidget.addClickFunc(clickFunc)

    def updateHeader(self):
        keys = self.obs.titles
        self.headerInfo = {}
        for eachKey in keys:
            self.headerInfo[eachKey] = self.obs.getFromHeader(eachKey)
            self.headerWindow.append('{}:\n{}\n'.format(eachKey,self.headerInfo[eachKey]))
            
    def getObsImage(self):
        self.startTime = float(self.textbox_startTime.text())
        self.stopTime = float(self.textbox_endTime.text())

        firstSec = self.startTime
        intTime = self.stopTime - firstSec
        
        startWvl = float(self.textbox_startWvl.text())
        endWvl = float(self.textbox_endWvl.text())

        if self.obs.wvlCalFile is None:
            self.imageParamsWindow.checkbox_getRawCount.setChecked(True)
            print 'setting getRawCount, since wvlcal is not loaded'
        else:
            self.obs.setWvlCutoffs(wvlLowerLimit=startWvl, wvlUpperLimit=endWvl)

        self.setObsState()
        paramsDict = self.imageParamsWindow.getParams()
        scaleToCps = paramsDict['otherParams'].pop('scaleToCps',False)
        imgDict = self.obs.getPixelCountImage(firstSec=firstSec,integrationTime=intTime,**paramsDict['obsParams'])

        self.image = imgDict['image']
        if scaleToCps:
            self.image = self.image / intTime
        self.plotArray(self.image,**paramsDict['plotParams'])

    def getObsParams(self):
        return self.imageParamsWindow.getParams()['obsParams']

    def setObsState(self):
        self.imageParamsWindow.setObsState(self.obs)

    def jumpToBeginning(self):
        firstSec = float(self.textbox_startTime.text())
        intTime = float(self.textbox_endTime.text())-firstSec
        newStartTime = 0.
        newEndTime = intTime
        self.textbox_startTime.setText(str(newStartTime))
        self.textbox_endTime.setText(str(newEndTime))
        self.getObsImage()

    def jumpToStartWvl(self):
        newStartWvl = 3000.
        wvlIncrement = float(self.textbox_wvlIncrement.text())
        newEndWvl = newStartWvl+wvlIncrement
        self.textbox_startWvl.setText(str(newStartWvl))
        self.textbox_endWvl.setText(str(newEndWvl))
        self.getObsImage()

    def jumpToEnd(self):
        firstSec = float(self.textbox_startTime.text())
        intTime = float(self.textbox_endTime.text())-firstSec
        newEndTime = self.headerInfo['exptime']
        newStartTime = newEndTime-intTime
        self.textbox_startTime.setText(str(newStartTime))
        self.textbox_endTime.setText(str(newEndTime))
        self.getObsImage()

    def jumpToEndWvl(self):
        wvlIncrement = float(self.textbox_wvlIncrement.text())
        newEndWvl = 11000.
        newStartWvl = newEndWvl-wvlIncrement
        self.textbox_startWvl.setText(str(newStartWvl))
        self.textbox_endWvl.setText(str(newEndWvl))
        self.getObsImage()

    def incrementForward(self):
        startTime = float(self.textbox_startTime.text())
        endTime = float(self.textbox_endTime.text())
        timeIncrement = float(self.textbox_timeIncrement.text())
        newStartTime = startTime + timeIncrement
        newEndTime = endTime + timeIncrement
        self.textbox_startTime.setText(str(newStartTime))
        self.textbox_endTime.setText(str(newEndTime))
        self.getObsImage()

    def incrementWvlForward(self):
        startWvl = float(self.textbox_startWvl.text())
        endWvl = float(self.textbox_endWvl.text())
        wvlIncrement = float(self.textbox_wvlIncrement.text())
        newStartWvl = startWvl + wvlIncrement
        newEndWvl = endWvl + wvlIncrement
        self.textbox_startWvl.setText(str(newStartWvl))
        self.textbox_endWvl.setText(str(newEndWvl))
        self.getObsImage()

    def incrementBack(self):
        startTime = float(self.textbox_startTime.text())
        endTime = float(self.textbox_endTime.text())
        timeIncrement = float(self.textbox_timeIncrement.text())
        newStartTime = startTime - timeIncrement
        newEndTime = endTime - timeIncrement
        self.textbox_startTime.setText(str(newStartTime))
        self.textbox_endTime.setText(str(newEndTime))
        self.getObsImage()

    def incrementWvlBack(self):
        startWvl = float(self.textbox_startWvl.text())
        endWvl = float(self.textbox_endWvl.text())
        wvlIncrement = float(self.textbox_wvlIncrement.text())
        newStartWvl = startWvl - wvlIncrement
        newEndWvl = endWvl - wvlIncrement
        self.textbox_startWvl.setText(str(newStartWvl))
        self.textbox_endWvl.setText(str(newEndWvl))
        self.getObsImage()

    def savePlot(self):
        file_choices = "PNG (*.png)|*.png"
        
        path = unicode(QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def aboutMessage(self):
        msg = """ Use to open and view obs file images
        """
        QtGui.QMessageBox.about(self, "ObsFile Viewer", msg.strip())

    def plotArray(self,*args,**kwargs):
        self.overplotTimeMask()
        self.arrayImageWidget.plotArray(*args,**kwargs)

    def hoverCanvas(self,event):
        col = int(round(event.xdata))
        row = int(round(event.ydata))
        if row < self.arrayImageWidget.nRow and col < self.arrayImageWidget.nCol:
            self.statusText.setText('(x,y,z) = ({:d},{:d},{})'.format(col,row,self.arrayImageWidget.image[row,col]))

    def loadObsFile(self,obsPath):
        print 'loading',obsPath
        self.obs = ObsFile(obsPath)
        self.fn = FileName(obsFile=self.obs)

        self.startTime = 0.
        self.stopTime = 1.
        firstImg = self.obs.getPixelCountImage(firstSec=self.startTime,integrationTime=(self.stopTime-self.startTime),getRawCount=True)
        self.image = firstImg['image']
        self.plotArray(self.image)
        self.updateHeader()
        self.label_obsFilename.setText(os.path.basename(obsPath))
        self.textbox_endTime.setText(str(self.headerInfo['exptime']))

    def obsMethod(self,method,*args,**kwargs):
        return getattr(self.obs,method)(*args,**kwargs)
    
    def loadObsFileWindow(self):
        LoadObsDialog(parent=self)

    def loadCalFiles(self):
        LoadCalsDialog(parent=self,obsTstamp=self.fn.tstamp)

    def getTimeMaskReasons(self):
        reasonEnum = self.obs.hotPixTimeMask.reasonEnum
        reasons = [reasonPair[0] for reasonPair in reasonEnum]
        selectedReasons = self.obs.hotPixTimeMask.enabledReasons
        return {'reasons':reasons,'selectedReasons':selectedReasons}

    def enableTimeMaskParams(self):
        if not self.obs.hotPixTimeMask is None:
            reasonDict = self.getTimeMaskReasons()
            self.imageParamsWindow.enableTimeMaskParams(reasons=reasonDict['reasons'],selectedReasons=reasonDict['selectedReasons'])

    def overplotTimeMask(self):
        if self.imageParamsWindow.checkbox_overplotMasked.isChecked():
            color = str(self.imageParamsWindow.combobox_overplotColor.currentText())
            timeMaskImage = self.makeTimeMaskImage()
            self.arrayImageWidget.plotOverlayImage(timeMaskImage,color=color)
        else:
            self.arrayImageWidget.removeOverlayImage()

    def makeTimeMaskImage(self):
        timeMaskImage = np.zeros((self.obs.nRow,self.obs.nCol))
        if not self.obs.hotPixTimeMask is None:
            for iRow in xrange(self.obs.nRow):
                for iCol in xrange(self.obs.nCol):
                    maskedIntervals = self.obs.getPixelBadTimes(iRow,iCol)
                    maskedTime = utils.intervalSize(maskedIntervals)
                    totalTime = self.stopTime-self.startTime
                    timeMaskImage[iRow,iCol] = 1.*maskedTime/totalTime
        return timeMaskImage
            

class ModelessWindow(QtGui.QDialog):
    def __init__(self,parent=None):
        super(ModelessWindow,self).__init__(parent=parent)
        self.parent=parent
        self.initUI()
        self._want_to_close = False

    def closeEvent(self, evt):
        if self._want_to_close:
            super(ModelessWindow, self).closeEvent(evt)
        else:
            evt.ignore()
            self.setVisible(False)

    def initUI(self):
        pass

class PlotWindow(QtGui.QDialog):
    def __init__(self,parent=None,plotId=0,selectedPixels=[]):
        super(PlotWindow,self).__init__(parent=parent)
        self.parent=parent
        self.id = plotId
        self.selectedPixels = selectedPixels
        self.initUI()

    def closeEvent(self,evt):
        super(PlotWindow,self).closeEvent(evt)

    def draw(self):
        self.fig.canvas.draw()

    def initUI(self):
        #first gui controls that apply to all modes
        self.checkbox_trackSelection = QtGui.QCheckBox('Plot selected pixel(s)',self)
        self.checkbox_trackSelection.setChecked(True)

        self.checkbox_trackTimes = QtGui.QCheckBox('Use main window times',self)
        self.checkbox_trackTimes.setChecked(True)
        self.connect(self.checkbox_trackTimes,QtCore.SIGNAL('stateChanged(int)'),self.changeTrackTimes)

        self.checkbox_clearPlot = QtGui.QCheckBox('Clear axes before plotting',self)
        self.checkbox_clearPlot.setChecked(True)

        self.button_drawPlot = QtGui.QPushButton('Plot',self)
        self.connect(self.button_drawPlot,QtCore.SIGNAL('clicked()'), self.updatePlot)
        self.dpi = 100
        self.fig = Figure((10.0, 3.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvasToolbar = NavigationToolbar(self.canvas, self)
        self.axes = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.07,right=.93,top=.93,bottom=0.15)

        self.combobox_plotType = QtGui.QComboBox(self)
        self.plotTypeStrs = ['Light Curve','Time Mask','Spectrum','Phase Histogram']
        self.combobox_plotType.addItems(self.plotTypeStrs)
        self.connect(self.combobox_plotType,QtCore.SIGNAL('activated(QString)'), self.changePlotType)


        #light curve controls
        self.textbox_intTime = QtGui.QLineEdit('1')
        self.textbox_intTime.setFixedWidth(50)

        self.button_maskEntirePixel = QtGui.QPushButton('Mask out entire pixel(s)')
        self.connect(self.button_maskEntirePixel,QtCore.SIGNAL('clicked()'), self.maskEntirePixels)

        self.combobox_timeMaskReason = QtGui.QComboBox(self)
        reasonDict = self.parent.getTimeMaskReasons()
        self.combobox_timeMaskReason.addItems(reasonDict['reasons'])

        lightCurveControlsBox = layoutBox('H',['Int Time',self.textbox_intTime,'s',1.,self.button_maskEntirePixel,'reason:',self.combobox_timeMaskReason,10.])
        self.lightCurveControlsGroup = QtGui.QGroupBox('Light Curve Controls',parent=self)
        self.lightCurveControlsGroup.setLayout(lightCurveControlsBox)


        #spectrum controls
        self.checkbox_divideWvlBinWidths = QtGui.QCheckBox('Divide by bin widths',self)
        self.checkbox_trackWvls = QtGui.QCheckBox('Use main window wavelengths',self)
        self.connect(self.checkbox_trackWvls,QtCore.SIGNAL('stateChanged(int)'),self.changeTrackWvls)
        self.textbox_startWvl = QtGui.QLineEdit(self.parent.textbox_startWvl.text())
        self.textbox_startWvl.setFixedWidth(100)
        self.textbox_endWvl = QtGui.QLineEdit(self.parent.textbox_endWvl.text())
        self.textbox_endWvl.setFixedWidth(100)

        self.wvlGroup = QtGui.QGroupBox('',parent=self)
        wvlBox = layoutBox('H',['Start Wavelength',self.textbox_startWvl,'Angstroms',1.,'End Wavelength',self.textbox_endWvl,'Angstroms',10.])
        self.wvlGroup.setLayout(wvlBox)
        spectrumControlsBox = layoutBox('H',[self.checkbox_divideWvlBinWidths,self.checkbox_trackWvls,self.wvlGroup,10.])

        self.spectrumControlsGroup = QtGui.QGroupBox('Spectrum Controls',parent=self)
        self.spectrumControlsGroup.setLayout(spectrumControlsBox)
        self.spectrumControlsGroup.setVisible(False)

        #phase hist controls
        self.checkbox_restrictPhaseRange = QtGui.QCheckBox('Restrict phase range')
        self.connect(self.checkbox_restrictPhaseRange,QtCore.SIGNAL('stateChanged(int)'),self.changeRestrictPhaseRange)
        self.checkbox_restrictPhaseRange.setChecked(False)

        self.checkbox_keepRawPhase = QtGui.QCheckBox('Raw units')

        self.textbox_phaseHistStart = QtGui.QLineEdit('')
        self.textbox_phaseHistEnd = QtGui.QLineEdit('')
        self.changeRestrictPhaseRange()

        self.textbox_phaseHistNBins = QtGui.QLineEdit('1')
        self.combobox_phaseHistType = QtGui.QComboBox(self)
        self.phaseHistTypeStrs = ['Peaks','Baselines','Peaks - Baselines']
        self.combobox_phaseHistType.addItems(self.phaseHistTypeStrs)
        self.combobox_phaseHistType.setCurrentIndex(2)

        phaseHistControlsBoxRow0 = layoutBox('H',[self.combobox_phaseHistType,1.,self.checkbox_keepRawPhase,1.,self.checkbox_restrictPhaseRange,10.])
        phaseHistControlsBoxRow1 = layoutBox('H',['Start phase',self.textbox_phaseHistStart,1.,'End phase',self.textbox_phaseHistEnd,1.,'Number of ADC units per bin',self.textbox_phaseHistNBins,10.])
        phaseHistControlsBox = layoutBox('V',[phaseHistControlsBoxRow0,phaseHistControlsBoxRow1])

        self.phaseHistControlsGroup = QtGui.QGroupBox('Phase Histogram Controls',parent=self)
        self.phaseHistControlsGroup.setLayout(phaseHistControlsBox)
        self.phaseHistControlsGroup.setVisible(False)

        #time controls
        self.textbox_startTime = QtGui.QLineEdit(self.parent.textbox_startTime.text())
        self.textbox_startTime.setFixedWidth(50)
        self.textbox_endTime = QtGui.QLineEdit(self.parent.textbox_endTime.text())
        self.textbox_endTime.setFixedWidth(50)
        self.timesGroup = QtGui.QGroupBox('',parent=self)
        timesBox = layoutBox('H',['Start Time',self.textbox_startTime,'s',1.,'End Time',self.textbox_endTime,'s',10.])
        self.timesGroup.setLayout(timesBox)
        self.timesGroup.setVisible(False)
        timesChoiceBox = layoutBox('H',[self.checkbox_trackTimes,self.timesGroup])

        checkboxBox = layoutBox('H',[self.checkbox_trackSelection,self.combobox_plotType,2.,self.button_drawPlot])
        controlsBox = layoutBox('H',[self.lightCurveControlsGroup,self.spectrumControlsGroup,self.phaseHistControlsGroup])

        mainBox = layoutBox('V',[checkboxBox,timesChoiceBox,self.checkbox_clearPlot,controlsBox,self.canvas,self.canvasToolbar])
        self.setLayout(mainBox)

    def changeTrackTimes(self):
        if self.checkbox_trackTimes.isChecked():
            self.timesGroup.setVisible(False)
        else:
            self.timesGroup.setVisible(True)

    def changeRestrictPhaseRange(self):
        if self.checkbox_restrictPhaseRange.isChecked():
            self.textbox_phaseHistStart.setEnabled(True)
            self.textbox_phaseHistEnd.setEnabled(True)
        else:
            self.textbox_phaseHistStart.setEnabled(False)
            self.textbox_phaseHistEnd.setEnabled(False)

    def changeTrackWvls(self):
        if self.checkbox_trackWvls.isChecked():
            self.wvlGroup.setVisible(False)
        else:
            self.wvlGroup.setVisible(True)

    def changePlotType(self,plotType):
        if plotType == 'Light Curve':
            self.lightCurveControlsGroup.setVisible(True)
            self.spectrumControlsGroup.setVisible(False)
            self.phaseHistControlsGroup.setVisible(False)
        elif plotType == 'Spectrum':
            self.lightCurveControlsGroup.setVisible(False)
            self.spectrumControlsGroup.setVisible(True)
            self.phaseHistControlsGroup.setVisible(False)
        elif plotType == 'Phase Histogram':
            self.lightCurveControlsGroup.setVisible(False)
            self.spectrumControlsGroup.setVisible(False)
            self.phaseHistControlsGroup.setVisible(True)
        elif plotType == 'Time Mask':
            self.lightCurveControlsGroup.setVisible(True)
            self.spectrumControlsGroup.setVisible(False)
            self.phaseHistControlsGroup.setVisible(False)

    def newPixelSelection(self,selectedPixels):
        if self.checkbox_trackSelection.isChecked():
            self.selectedPixels = selectedPixels
            self.updatePlot()
        

    def updatePlot(self):
        plotType = self.combobox_plotType.currentText()
        if self.checkbox_clearPlot.isChecked() or plotType != self.lastPlotType:
            self.axes.cla()
        self.parent.setObsState()
        if plotType == 'Light Curve':
            self.plotLightCurve()
        elif plotType == 'Spectrum':
            self.plotSpectrum()
        elif plotType == 'Phase Histogram':
            self.plotPhaseHist()
        elif plotType == 'Time Mask':
            self.plotTimeMask()
        self.lastPlotType = plotType
        self.draw()
        print 'plot updated'

    def plotTimeMask(self):
        self.parent.obsMethod('switchOffHotPixTimeMask')
        self.plotLightCurve()
        for col,row in self.selectedPixels:
            maskedIntervals = self.parent.obsMethod('getPixelBadTimes',pixelRow=row,pixelCol=col)
            print maskedIntervals
            for inter in maskedIntervals:
                self.axes.axvspan(inter[0],inter[1],facecolor='r',alpha=0.3)
        

    def plotLightCurve(self,getRaw=False):
        paramsDict = self.parent.getObsParams()
        if self.checkbox_trackTimes.isChecked():
            startTime = float(self.parent.textbox_startTime.text())
            endTime = float(self.parent.textbox_endTime.text())
        else:
            startTime = float(self.textbox_startTime.text())
            endTime = float(self.textbox_endTime.text())
        histIntTime = float(self.textbox_intTime.text())
        firstSec = startTime
        duration = endTime-startTime
        histBinEdges = np.arange(startTime,endTime,histIntTime)
        hists = []
        if histBinEdges[-1]+histIntTime == endTime:
            histBinEdges = np.append(histBinEdges,endTime)
        for col,row in self.selectedPixels:
            returnDict = self.parent.obs.getPixelWvlList(iRow=row,iCol=col,firstSec=firstSec,integrationTime=duration,excludeBad=(not paramsDict['getRawCount']))
            timestamps = returnDict['timestamps']
            hist,_ = np.histogram(timestamps,bins=histBinEdges)
            hists.append(hist)

        hists = np.array(hists)
        hist = np.sum(hists,axis=0)

        binWidths = np.diff(histBinEdges)
        self.lightCurve = 1.*hist/binWidths
        plotHist(self.axes,histBinEdges,self.lightCurve)
        self.axes.set_xlabel('time (s)')
        self.axes.set_ylabel('counts per sec')
            
    def plotSpectrum(self):
        if self.checkbox_trackTimes.isChecked():
            startTime = float(self.parent.textbox_startTime.text())
            endTime = float(self.parent.textbox_endTime.text())
        else:
            startTime = float(self.textbox_startTime.text())
            endTime = float(self.textbox_endTime.text())
        if self.checkbox_trackWvls.isChecked():
            startWvl = float(self.parent.textbox_startWvl.text())
            endWvl = float(self.parent.textbox_endWvl.text())
        else:
            startWvl = float(self.textbox_startWvl.text())
            endWvl = float(self.textbox_endWvl.text())

        firstSec = startTime
        duration = endTime-startTime
        paramsDict = self.parent.getObsParams()

        spectrums = []
        for col,row in self.selectedPixels:
            returnDict = self.parent.obs.getPixelSpectrum(pixelRow=row,pixelCol=col,firstSec=firstSec,integrationTime=duration,wvlStart=startWvl,wvlStop=endWvl,weighted=paramsDict['weighted'],fluxWeighted=paramsDict['fluxWeighted'])
            spectrum = returnDict['spectrum']
            wvlBinEdges = returnDict['wvlBinEdges']
            
            spectrums.append(spectrum)

        spectrums = np.array(spectrums)
        self.spectrum = np.sum(spectrums,axis=0)

        binWidths = np.diff(wvlBinEdges)
        if self.checkbox_divideWvlBinWidths.isChecked():
            self.spectrum = 1. * self.spectrum / binWidths
            ylabel = 'Counts / Angstrom'
        else:
            ylabel = 'Counts'
        
        plotHist(self.axes,wvlBinEdges,self.spectrum)
        self.axes.set_xlabel('Wavelength (Angstrom)')
        self.axes.set_ylabel(ylabel)

    def plotPhaseHist(self):
        phaseHistType = self.combobox_phaseHistType.currentText()
        if self.checkbox_trackTimes.isChecked():
            startTime = float(self.parent.textbox_startTime.text())
            endTime = float(self.parent.textbox_endTime.text())
        else:
            startTime = float(self.textbox_startTime.text())
            endTime = float(self.textbox_endTime.text())
        histIntTime = float(self.textbox_intTime.text())
        firstSec = startTime
        duration = endTime-startTime
        hists = []
        nBinsPerUnit = int(self.textbox_phaseHistNBins.text())
        histBinEdges = None
        scaleToDegrees = 180./np.pi*1./2**9
        zeroOffset = 2**11 # phase is +/-4 radians, where 0 radians is represented by bin value 2048
        if self.checkbox_restrictPhaseRange.isChecked():
            range = np.array([float(self.textbox_phaseHistStart.text()),float(self.textbox_phaseHistEnd.text())])
            if not self.checkbox_keepRawPhase.isChecked():
                range = range / scaleToDegrees
        else:
            range = None
        for col,row in self.selectedPixels:
            returnDict = self.parent.obs.getTimedPacketList(iRow=row,iCol=col,firstSec=firstSec,integrationTime=duration)
            timestamps = returnDict['timestamps']
            peakHeights = returnDict['peakHeights']
            baselines = returnDict['baselines']
            if phaseHistType == 'Peaks':
                listToHist = peakHeights-zeroOffset
            elif phaseHistType == 'Baselines':
                listToHist = baselines-zeroOffset
            elif phaseHistType == 'Peaks - Baselines':
                listToHist = peakHeights-baselines

            if histBinEdges == None:
                if not self.checkbox_restrictPhaseRange.isChecked():
                    range = np.array([np.min(listToHist),np.max(listToHist)])
                bins = (range[1]-range[0])//nBinsPerUnit
            else:
                bins = histBinEdges

            hist,histBinEdges = np.histogram(listToHist,bins=bins,range=range)
            hists.append(hist)

        hists = np.array(hists)
        self.hist = np.sum(hists,axis=0)

        binWidths = np.diff(histBinEdges)
            
        if self.checkbox_keepRawPhase.isChecked():
            xlabel = 'phase'
        else:
            histBinEdges = histBinEdges * scaleToDegrees
            xlabel = 'phase (${}^\circ$)'
        plotHist(self.axes,histBinEdges,self.hist)
        self.axes.set_xlabel(xlabel)
        self.axes.set_ylabel('counts')

    def maskEntirePixels(self):
        reason = str(self.combobox_timeMaskReason.currentText())
        reasonDict = self.parent.getTimeMaskReasons()
        reply = QtGui.QMessageBox.question(self, 'Confirm',
                'Mark pixels (x,y)={} with tag \'{}\'?'.format(self.selectedPixels,reason), QtGui.QMessageBox.Yes |
                QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            if not self.parent.obs.hotPixFileName is None:
                timeMaskPath = self.parent.obs.hotPixFileName
            else:
                raise AttributeError('obs file does not have a timeMask loaded')
            self.parent.obs.hotPixFile.close()
            for col,row in self.selectedPixels:
                removePixel(timeMaskPath=timeMaskPath,pixelRow=row,
                            pixelCol=col,reason=reason)
                print 'pixel (x,y)=({},{}) tagged'.format(col,row)
            self.parent.obsMethod('loadHotPixCalFile',timeMaskPath,reasons=reasonDict['selectedReasons'])
     
class HeaderWindow(ModelessWindow):
    def append(self,text):
        self.textArea.append(text)

    def setHeader(self,newHeaderText):
        self.textArea.setText(newHeaderText)

    def initUI(self):
        self.textArea = QtGui.QTextEdit('')
        self.textArea.setReadOnly(True)
        mainBox = layoutBox('V',[self.textArea])
        self.setLayout(mainBox)

class ImageParamsWindow(ModelessWindow):
    def initUI(self):
        self.checkbox_getRawCount = QtGui.QCheckBox('getRawCount',self)
        self.checkbox_weighted = QtGui.QCheckBox('weighted',self)
        self.checkbox_fluxWeighted = QtGui.QCheckBox('fluxWeighted',self)
        self.checkbox_scaleByEffInt = QtGui.QCheckBox('scaleByEffInt',self)
        self.checkbox_applyTimeMask = QtGui.QCheckBox('apply time mask',self)
        self.checkbox_applyTimeMask.setEnabled(False)
        self.checkbox_applyCosmicMask = QtGui.QCheckBox('apply cosmic mask',self)
        self.checkbox_applyCosmicMask.setEnabled(False)

        self.group_timeMaskReasons = QtGui.QGroupBox('Reasons to mask times',self)
        self.group_timeMaskReasons.setVisible(False)
        self.timeMaskReasonsLoaded = False

        obsBox = layoutBox('V',[self.checkbox_getRawCount,
                    self.checkbox_weighted,self.checkbox_fluxWeighted,self.checkbox_scaleByEffInt,
                    self.checkbox_applyTimeMask,self.group_timeMaskReasons,self.checkbox_applyCosmicMask])
        obsGroup = QtGui.QGroupBox('ObsFile.getPixelCountImage parameters',self)
        obsGroup.setLayout(obsBox)

        self.combobox_cmap = QtGui.QComboBox(self)
        self.cmapStrs = ['hot','gray','jet','gnuplot2','Paired']
        self.combobox_cmap.addItems(self.cmapStrs)

        plotArrayBox = layoutBox('V',[self.combobox_cmap])
        plotArrayGroup = QtGui.QGroupBox('plotArray parameters',self)
        plotArrayGroup.setLayout(plotArrayBox)

        
        self.checkbox_scaleToCps = QtGui.QCheckBox('scale as counts per sec',self)
        self.checkbox_overplotMasked = QtGui.QCheckBox('Show time masked pixels (over total time)',self)
        self.connect(self.checkbox_overplotMasked,QtCore.SIGNAL('stateChanged(int)'),self.parent.overplotTimeMask)
        self.combobox_overplotColor = QtGui.QComboBox(self)
        self.overplotColorStrs = ['blue','green','red','cyan','magenta','yellow','black','white','0.5']
        self.combobox_overplotColor.addItems(self.overplotColorStrs)

        otherParamsBox = layoutBox('V',[self.checkbox_scaleToCps,self.checkbox_overplotMasked,self.combobox_overplotColor])
        otherParamsGroup = QtGui.QGroupBox('Other parameters',self)
        otherParamsGroup.setLayout(otherParamsBox)

        mainBox = layoutBox('V',[plotArrayGroup,otherParamsGroup,obsGroup])

        self.setLayout(mainBox)



    def enableTimeMaskParams(self,reasons=[],selectedReasons=[]):
        self.checkboxes_timeMaskReasons = {}
        vbox = QtGui.QVBoxLayout()
        if len(reasons) > 0:
            for reason in reasons:
                checkbox = QtGui.QCheckBox(reason,self)
                self.checkboxes_timeMaskReasons[reason] = checkbox
                vbox.addWidget(checkbox)
                if reason in selectedReasons:
                    checkbox.setChecked(True)
            self.group_timeMaskReasons.setLayout(vbox)
            self.group_timeMaskReasons.setVisible(True)
        self.checkbox_applyTimeMask.setEnabled(True)
        self.checkbox_applyTimeMask.setChecked(True)
        self.timeMaskReasonsLoaded = True
        
    def getCheckedTimeMaskReasons(self):
        allReasons = self.checkboxes_timeMaskReasons.keys()
        checkedReasons = [reason for reason in allReasons if (self.checkboxes_timeMaskReasons[reason]).isChecked()]
        return checkedReasons

    def enableCosmicMaskParams(self):
        self.checkbox_applyCosmicMask.setEnabled(True)

    def setObsState(self,obs):
        if self.checkbox_applyTimeMask.isChecked():
            reasons = self.getCheckedTimeMaskReasons()
            obs.switchOnHotPixTimeMask(reasons=reasons)
        else:
            obs.switchOffHotPixTimeMask()

        if self.checkbox_applyCosmicMask.isChecked():
            obs.switchOnCosmicTimeMask()
        else:
            obs.switchOffCosmicTimeMask()
            
    def getParams(self):
        obsParamsDict = {}
        plotParamsDict = {}
        otherParamsDict = {}
        obsParamsDict['getRawCount'] = self.checkbox_getRawCount.isChecked()
        obsParamsDict['weighted'] = self.checkbox_weighted.isChecked()
        obsParamsDict['fluxWeighted'] = self.checkbox_fluxWeighted.isChecked()
        obsParamsDict['scaleByEffInt'] = self.checkbox_scaleByEffInt.isChecked()
        cmapStr = str(self.combobox_cmap.currentText())
        cmap = getattr(matplotlib.cm,cmapStr)
        if cmapStr != 'gray':
            cmap.set_bad('0.15')
        plotParamsDict['cmap']=cmap
        otherParamsDict['scaleToCps']=self.checkbox_scaleToCps.isChecked()

        outDict = {}
        outDict['obsParams'] = obsParamsDict
        outDict['plotParams'] = plotParamsDict
        outDict['otherParams'] = otherParamsDict
        return outDict
    
class ArrayImageWidget(QtGui.QWidget):
    def __init__(self,parent=None,hoverCall=None):
        super(ArrayImageWidget,self).__init__(parent=parent)
        self.parent=parent
        # Create the mpl Figure and FigCanvas objects. 
        self.hoverCall = hoverCall
        self.selectPixelsMode = 'singlePixel'
        self.selectionPatches = []
        self.selectedPixels = []
        self.overlayImage = None
        self.initUI()
        

    def initUI(self):
        self.dpi = 100
        self.fig = Figure((5.0, 5.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.axes = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=0.07,right=.93,top=.93,bottom=0.07)
        self.plotArray(np.arange(9).reshape((3,3)))

        self.clickFuncs = []
        cid = self.fig.canvas.mpl_connect('scroll_event', self.scrollColorBar)
        cid = self.fig.canvas.mpl_connect('button_press_event', self.clickColorBar)
        cid = self.fig.canvas.mpl_connect('motion_notify_event', self.hoverCanvas)
        cid = self.fig.canvas.mpl_connect('button_press_event', self.clickCanvas)
        canvasBox = layoutBox('V',[self.canvas,])
        self.setLayout(canvasBox)

    def plotOverlayImage(self,image,color='green',**kwargs):
        self.overlayImage = image
        self.overlayImageKwargs = kwargs
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap',[color,color],256)
        cmap._init() # create the _lut array, with rgba values
        alphas = np.linspace(0, 1., cmap.N+3)
        cmap._lut[:,-1] = alphas
        self.overlayCmap = cmap
        self.drawOverlayImage()

    def removeOverlayImage(self):
        self.overlayImage = None
        try:
            self.handleOverlayMatshow.remove()
        except:
            pass
        self.draw()

    def drawOverlayImage(self):
        try:
            self.handleOverlayMatshow.remove()
        except:
            pass
        self.handleOverlayMatshow = self.axes.matshow(self.overlayImage,cmap=self.overlayCmap,origin='lower',**self.overlayImageKwargs)
        self.draw()

    def plotArray(self,image,normNSigma=3,title='',**kwargs):
        self.image = image
        self.imageShape = np.shape(image)
        self.nRow = self.imageShape[0]
        self.nCol = self.imageShape[1]
        if not 'vmax' in kwargs:
            goodImage = image[np.isfinite(image)]
            kwargs['vmax'] = np.mean(goodImage)+normNSigma*np.std(goodImage)
        if not 'cmap' in kwargs:
            defaultCmap=matplotlib.cm.hot
            defaultCmap.set_bad('0.15')
            kwargs['cmap'] = defaultCmap
        if not 'origin' in kwargs:
            kwargs['origin'] = 'lower'

        self.fig.clf()
        self.selectionPatches = []
        self.axes = self.fig.add_subplot(111)
        self.axes.set_title(title)

        self.matshowKwargs = kwargs
        self.handleMatshow = self.axes.matshow(image,**kwargs)
        self.fig.cbar = self.fig.colorbar(self.handleMatshow)
        if not self.overlayImage is None:
            self.drawOverlayImage()
        else:
            self.draw()
        print 'image drawn'

    def drawSelections(self):
        for patch in self.selectionPatches:
            patch.remove()
        self.selectionPatches = []
        
        for pixelCoord in self.selectedPixels:
            lowerLeftCorner = tuple(np.subtract(pixelCoord, (0.5,0.5)))
            patch = matplotlib.patches.Rectangle(xy=lowerLeftCorner,width=1.,
                    height=1.,edgecolor='blue',facecolor='none')
            self.selectionPatches.append(patch)
            self.axes.add_patch(patch)

    def addClickFunc(self,clickFunc):
        self.clickFuncs.append(clickFunc)

    def emitNewSelection(self):
        self.emit(QtCore.SIGNAL('newPixelSelection(PyQt_PyObject)'),self.selectedPixels)

    def clickCanvas(self,event):
        if event.inaxes is self.axes:
            col = round(event.xdata)
            row = round(event.ydata)
            if self.selectPixelsMode == 'singlePixel':
                self.selectedPixels = [(col,row)]
                self.draw()
                self.emitNewSelection()
            for func in self.clickFuncs:
                func(row=row,col=col)


    def hoverCanvas(self,event):
        if event.inaxes is self.axes:
            self.hoverCall(event)

    def scrollColorBar(self, event):
        if event.inaxes is self.fig.cbar.ax:
            increment=0.05
            currentClim = self.fig.cbar.mappable.get_clim()
            currentRange = currentClim[1]-currentClim[0]
            if event.button == 'up':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]+increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]+increment*currentRange)
            if event.button == 'down':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]-increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]-increment*currentRange)
            self.fig.cbar.mappable.set_clim(newClim)
            self.fig.canvas.draw()

    def clickColorBar(self,event):
        if event.inaxes is self.fig.cbar.ax:
            self.fig.currentClim = self.fig.cbar.mappable.get_clim()
            lower = self.fig.currentClim[0]
            upper = self.fig.currentClim[1]
            fraction = event.ydata
            currentRange = upper-lower
            clickedValue = lower+fraction*currentRange
            extrapolatedValue = lower+event.ydata*currentRange
            if event.button == 1:
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (clickedValue,upper)
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (lower,clickedValue)
            if event.button == 3:
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = ((lower-fraction*upper)/(1.-fraction),upper)
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (lower,lower+currentRange/fraction)
            self.fig.cbar.mappable.set_clim(newClim)
            self.fig.canvas.draw()

    def draw(self):
        self.drawSelections()
        self.fig.canvas.draw()

    def addClickFunc(self,clickFunc):
        self.clickFuncs.append(clickFunc)


class LoadWidget(QtGui.QWidget):
    def __init__(self,parent=None,fileNameMethod='obs',promptStr='',obsTstamp=''):
        super(LoadWidget,self).__init__(parent=parent)
        self.parent=parent
        self.fileNameMethod = fileNameMethod
        self.promptStr=promptStr
        if obsTstamp != '':
            self.calLookupFile = CalLookupFile()
            self.obsTstamp = obsTstamp
            initialLoadPath = getattr(self.calLookupFile,fileNameMethod)(self.obsTstamp)
            
            print initialLoadPath
            initialParams = self.calLookupFile.getComponents(self.obsTstamp,fileNameMethod)
        else:
            initialParams = {}
            initialLoadPath = ''
        self.initUI(initialLoadPath=initialLoadPath,initialParams=initialParams)

    def initUI(self,initialLoadPath='',initialParams={}):
        self.textbox_run = QtGui.QLineEdit(initialParams.get('run'))
        self.textbox_date = QtGui.QLineEdit(initialParams.get('date'))
        self.textbox_tstamp = QtGui.QLineEdit(initialParams.get('tstamp'))
        self.button_loadFilename = QtGui.QPushButton('Lookup FileName Path')

        self.button_openDialog = QtGui.QPushButton('Open Path')
        self.textbox_filename = QtGui.QLineEdit(initialLoadPath)
        #self.setFilenameFromParams()

#        self.textbox_run.setMaximumWidth(200)
#        self.textbox_date.setMaximumWidth(200)
#        self.textbox_tstamp.setMaximumWidth(400)
#        self.textbox_filename.setMaximumWidth(1000)
        self.button_openDialog.setMaximumWidth(100)

        fnBox = layoutBox('H',['run:',self.textbox_run,'date:',self.textbox_date,
                    'timestamp:',self.textbox_tstamp,self.button_loadFilename])
        dialogBox = layoutBox('H',[self.button_openDialog,self.textbox_filename])
        mainBox = layoutBox('V',[self.promptStr,fnBox,dialogBox])

        self.setLayout(mainBox)

        self.connect(self.button_loadFilename,QtCore.SIGNAL('clicked()'), self.setFilenameFromParams)
        self.connect(self.button_openDialog,QtCore.SIGNAL('clicked()'), self.setFilenameFromDialog)


    def setFilenameFromDialog(self):
        run = str(self.textbox_run.text())
        date = str(self.textbox_date.text())
        initialPath = os.path.join(os.path.join(os.environ['MKID_RAW_PATH'],run),date)
        filename = QtGui.QFileDialog.getOpenFileName(self,'Open File',initialPath,"(*.h5)")
        self.textbox_filename.setText(filename)

    def setFilenameFromParams(self):
        run = str(self.textbox_run.text())
        date = str(self.textbox_date.text())
        tstamp = str(self.textbox_tstamp.text())
        if run != '' or date != '' or tstamp != '':
            fn = FileName(run=run,date=date,tstamp=tstamp)
            path = getattr(fn,self.fileNameMethod)()
            self.textbox_filename.setText(path)

class LoadObsDialog(QtGui.QDialog):
    def __init__(self,parent=None):
        super(LoadObsDialog,self).__init__(parent=parent)
        self.parent=parent
        self.initUI()
        self.connectButtons()
        self.exec_()

    def initUI(self):

        self.loadObsWidget = LoadWidget(self,'obs','Choose the Obs file')

        self.button_load = QtGui.QPushButton('Load Obs')
        self.button_cancel = QtGui.QPushButton('Cancel')
        buttonBox = layoutBox('H',[self.button_cancel,self.button_load])
        mainBox = layoutBox('V',[self.loadObsWidget,buttonBox])
        self.setLayout(mainBox)

    def connectButtons(self):
        self.connect(self.button_load,QtCore.SIGNAL('clicked()'),self.load)
        self.connect(self.button_cancel,QtCore.SIGNAL('clicked()'),self.close)

    def load(self):
        filename = str(self.loadObsWidget.textbox_filename.text())
        self.parent.loadObsFile(filename)
        self.close()


class LoadCalsDialog(QtGui.QDialog):
    def __init__(self,parent=None,obsTstamp=''):
        super(LoadCalsDialog,self).__init__(parent=parent)
        self.parent=parent
        self.obsTstamp = obsTstamp
        self.initUI()
        self.connectButtons()
        self.exec_()

    def initUI(self):

        self.loadWvlWidget = LoadWidget(self,'calSoln','Choose the wavelength cal',obsTstamp=self.obsTstamp)
        self.loadFlatWidget = LoadWidget(self,'flatSoln','Choose the flat cal',obsTstamp=self.obsTstamp)
        self.loadFluxWidget = LoadWidget(self,'fluxSoln','Choose the flux cal',obsTstamp=self.obsTstamp)
        self.loadTimeMaskWidget = LoadWidget(self,'timeMask','Choose the time mask file',obsTstamp=self.obsTstamp)
        self.loadTimeAdjustmentWidget = LoadWidget(self,'timeAdjustments','Choose the time adjustment file',obsTstamp=self.obsTstamp)
        self.loadCosmicWidget = LoadWidget(self,'cosmicMask','Choose the cosmic mask file',obsTstamp=self.obsTstamp)
        self.loadBeammapWidget = LoadWidget(self,'beammap','Choose the beammap file',obsTstamp=self.obsTstamp)

        self.button_loadWvl = QtGui.QPushButton('Load Wvl Cal')
        self.button_loadFlat = QtGui.QPushButton('Load Flat Cal')
        self.button_loadFlux = QtGui.QPushButton('Load Flux Cal')
        self.button_loadTimeMask = QtGui.QPushButton('Load Time Mask')
        self.button_loadTimeAdjustment = QtGui.QPushButton('Load Time Adjustments')
        self.button_loadCosmic = QtGui.QPushButton('Load Cosmic Mask')
        self.button_loadBeammap = QtGui.QPushButton('Load Beammap')

        self.button_loadAll = QtGui.QPushButton('Load All Cals')

        wvlBox = layoutBox('H',[self.loadWvlWidget,self.button_loadWvl])
        flatBox = layoutBox('H',[self.loadFlatWidget,self.button_loadFlat])
        fluxBox = layoutBox('H',[self.loadFluxWidget,self.button_loadFlux])
        timeMaskBox = layoutBox('H',[self.loadTimeMaskWidget,self.button_loadTimeMask])
        timeAdjustBox = layoutBox('H',[self.loadTimeAdjustmentWidget,self.button_loadTimeAdjustment])
        cosmicBox = layoutBox('H',[self.loadCosmicWidget,self.button_loadCosmic])
        beammapBox = layoutBox('H',[self.loadBeammapWidget,self.button_loadBeammap])

        mainBox = layoutBox('V',[wvlBox,flatBox,fluxBox,timeMaskBox,timeAdjustBox,cosmicBox,beammapBox,self.button_loadAll])
        self.setLayout(mainBox)

    def connectButtons(self):
        self.connect(self.button_loadWvl,QtCore.SIGNAL('clicked()'),partial(self.loadCal,self.loadWvlWidget,'loadWvlCalFile'))
        self.connect(self.button_loadFlat,QtCore.SIGNAL('clicked()'),partial(self.loadCal,self.loadFlatWidget,'loadFlatCalFile'))
        self.connect(self.button_loadFlux,QtCore.SIGNAL('clicked()'),partial(self.loadCal,self.loadFluxWidget,'loadFluxCalFile'))
        self.connect(self.button_loadTimeMask,QtCore.SIGNAL('clicked()'),partial(self.loadCal,self.loadTimeMaskWidget,'loadHotPixCalFile'))
        self.connect(self.button_loadTimeAdjustment,QtCore.SIGNAL('clicked()'),partial(self.loadCal,self.loadTimeAdjustmentWidget,'loadTimeAdjustmentFile'))
        self.connect(self.button_loadCosmic,QtCore.SIGNAL('clicked()'),partial(self.loadCal,self.loadCosmicWidget,'loadCosmicMaskFile'))
        self.connect(self.button_loadBeammap,QtCore.SIGNAL('clicked()'),partial(self.loadCal,self.loadBeammapWidget,'loadBeammapFile'))
        self.connect(self.button_loadAll,QtCore.SIGNAL('clicked()'),self.loadAllCals)

    def loadCal(self,calWidget,obsMethod,*args,**kwargs):
        path = str(calWidget.textbox_filename.text())
        if path != '':
            print 'loading',path
            self.parent.obsMethod(obsMethod,path,*args,**kwargs)
    
    def loadTimeMask(self):
        self.loadCal(self.loadTimeMaskWidget,'loadHotPixCalFile',reasons=['hot pixel','unknown'])
        path = str(self.loadTimeMaskWidget.textbox_filename.text())
        self.parent.enableTimeMaskParams()
        
    def loadAllCals(self):
        self.loadCal(self.loadWvlWidget,'loadWvlCalFile')
        self.loadCal(self.loadFlatWidget,'loadFlatCalFile')
        self.loadCal(self.loadFluxWidget,'loadFluxCalFile')
        self.loadTimeMask()
        self.loadCal(self.loadTimeAdjustmentWidget,'loadTimeAdjustmentFile')
        self.loadCal(self.loadCosmicWidget,'loadCosmicMaskFile')
        self.loadCal(self.loadBeammapWidget,'loadBeammapFile')
        self.parent.imageParamsWindow.checkbox_getRawCount.setChecked(False)
        print 'setting getRawCount=False'
        self.close()

#gui functions
def addActions(target, actions):
    for action in actions:
        if action is None:
            target.addSeparator()
        else:
            target.addAction(action)

def createAction(  gui, text, slot=None, shortcut=None, 
                    icon=None, tip=None, checkable=False, 
                    signal="triggered()"):
    action = QtGui.QAction(text, gui)
    if icon is not None:
        action.setIcon(QIcon(":/%s.png" % icon))
    if shortcut is not None:
        action.setShortcut(shortcut)
    if tip is not None:
        action.setToolTip(tip)
        action.setStatusTip(tip)
    if slot is not None:
        gui.connect(action, QtCore.SIGNAL(signal), slot)
    if checkable:
        action.setCheckable(True)
    return action

def layoutBox(type,elements):
    if type == 'vertical' or type == 'V':
        box = QtGui.QVBoxLayout()
    elif type == 'horizontal' or type == 'H':
        box = QtGui.QHBoxLayout()
    else:
        raise TypeError('type should be one of [\'vertical\',\'horizontal\',\'V\',\'H\']')

    for element in elements:
        try:
            box.addWidget(element)
        except:
            try:
                box.addLayout(element)
            except:
                try:
                    box.addStretch(element)
                except:
                    try:
                        label = QtGui.QLabel(element)
                        box.addWidget(label)
                        #label.adjustSize()
                    except:
                        print 'could\'t add {} to layout box'.format(element)
    return box

def plotHist(ax,histBinEdges,hist,**kwargs):
    ax.plot(histBinEdges,np.append(hist,hist[-1]),drawstyle='steps-post',**kwargs)


if __name__ == "__main__":
    kwargs = {}
    if len(sys.argv) > 1:
        if sys.argv[1] == '-h' or sys.argv[1] == '--help':
            print 'Usage: {} obsFilePath/obsTstamp'.format(sys.argv[0])
        elif os.path.exists(sys.argv[1]):
            kwargs['obsPath'] = sys.argv[1]
        else:
            kwargs['obsTstamp'] = sys.argv[1]
    else:
        obsPath = None
    form = ObsFileViewer(**kwargs)
    form.show()



