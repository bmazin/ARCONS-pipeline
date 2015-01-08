import sys, os
import numpy as np
from PyQt4 import QtCore
from PyQt4 import QtGui
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

from util.FileName import FileName
from util.ObsFile import ObsFile

class ObsFileViewer(QtGui.QMainWindow):
    def __init__(self, obsPath=None,showMe=True):
        self.app = QtGui.QApplication([])
        self.app.setStyle('plastique')
        super(self.__class__,self).__init__()
        self.setWindowTitle('ObsFile Viewer')
        self.createWidgets()
        self.createWindows()
        self.createMenu()
        self.createStatusBar()
        self.arrangeMainFrame()
        self.connectControls()

        if not obsPath is None:
            self.loadObsFile(obsPath)
            
        self.show()

    def show(self):
        super(self.__class__,self).show()
        self.app.exec_()

    def createWindows(self):
        self.headerWindow = HeaderWindow(self)
        self.imageParamsWindow = ImageParamsWindow(self)

    def createWidgets(self):
        self.mainFrame = QtGui.QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        self.dpi = 100
        self.fig = Figure((5.0, 5.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.mainFrame)
        self.axes = self.fig.add_subplot(111)

        # Create the navigation toolbar, tied to the canvas
        #self.canvasToolbar = NavigationToolbar(self.canvas, self.mainFrame)

        self.label_obsFilename = QtGui.QLabel('No obs loaded')

        self.textbox_startTime = QtGui.QLineEdit('0')
        self.textbox_startTime.setMaximumWidth(50)
        self.label_startToEnd = QtGui.QLabel('s to')
        self.textbox_endTime = QtGui.QLineEdit('0')
        self.textbox_endTime.setMaximumWidth(50)
        self.label_unitsAfterEndTime = QtGui.QLabel('s')

        self.button_drawPlot = QtGui.QPushButton('Plot')


        self.textbox_timeIncrement = QtGui.QLineEdit('10')
        self.textbox_timeIncrement.setMaximumWidth(50)
        self.button_jumpToBeginning = QtGui.QPushButton('|<')
        self.button_jumpToEnd = QtGui.QPushButton('>|')
        self.button_incrementBack = QtGui.QPushButton('<')
        self.button_incrementForward = QtGui.QPushButton('>')

        self.button_jumpToBeginning.setMaximumWidth(30)
        self.button_jumpToEnd.setMaximumWidth(30)
        self.button_incrementBack.setMaximumWidth(30)
        self.button_incrementForward.setMaximumWidth(30)


    def arrangeMainFrame(self):

        canvasBox = layoutBox('V',(self.canvas,))#self.canvasToolbar

        timeBoundsBox = layoutBox('H',(3,self.textbox_startTime,'s to',self.textbox_endTime,'s',1,self.button_drawPlot,3))
        incrementControlsBox = layoutBox('H',(1,self.button_jumpToBeginning,self.button_incrementBack,
                    self.textbox_timeIncrement,'s',self.button_incrementForward,self.button_jumpToEnd,1))
        mainBox = layoutBox('V',(self.label_obsFilename,canvasBox,timeBoundsBox,incrementControlsBox))

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
        addActions(self.windowsMenu,(headerAction,imgParamsAction))
        
        self.helpMenu = self.menuBar().addMenu("&Help")
        aboutAction = createAction(self,"&About", 
            shortcut='F1', slot=self.aboutMessage)
        
        addActions(self.helpMenu, (aboutAction,))

    def connectControls(self):
        self.connect(self.button_drawPlot,QtCore.SIGNAL('clicked()'), self.drawPlot)
        self.connect(self.button_jumpToBeginning,QtCore.SIGNAL('clicked()'), self.jumpToBeginning)
        self.connect(self.button_jumpToEnd,QtCore.SIGNAL('clicked()'), self.jumpToEnd)
        self.connect(self.button_incrementForward,QtCore.SIGNAL('clicked()'), self.incrementForward)
        self.connect(self.button_incrementBack,QtCore.SIGNAL('clicked()'), self.incrementBack)

    def updateHeader(self):
        keys = self.obs.titles
        self.headerInfo = {}
        for eachKey in keys:
            self.headerInfo[eachKey] = self.obs.getFromHeader(eachKey)
            self.headerWindow.append('{}:\n{}\n'.format(eachKey,self.headerInfo[eachKey]))
            
    def drawPlot(self):
        firstSec = float(self.textbox_startTime.text())
        intTime = float(self.textbox_endTime.text())-firstSec
        
        if self.obs.wvlCalFile is None:
            self.imageParamsWindow.checkbox_getRawCount.setCheckState(True)
            print 'setting getRawCount, since wvlcal is not loaded'
        paramsDict = self.imageParamsWindow.getParams()
        
        imgDict = self.obs.getPixelCountImage(firstSec=firstSec,integrationTime=intTime,**paramsDict['obsParams'])

        self.image = imgDict['image']
        self.plotArray(self.image,**paramsDict['plotParams'])

    def jumpToBeginning(self):
        firstSec = float(self.textbox_startTime.text())
        intTime = float(self.textbox_endTime.text())-firstSec
        newStartTime = 0.
        newEndTime = intTime
        self.textbox_startTime.setText(str(newStartTime))
        self.textbox_endTime.setText(str(newEndTime))
        self.drawPlot()

    def jumpToEnd(self):
        firstSec = float(self.textbox_startTime.text())
        intTime = float(self.textbox_endTime.text())-firstSec
        newEndTime = self.headerInfo['exptime']
        newStartTime = newEndTime-intTime
        self.textbox_startTime.setText(str(newStartTime))
        self.textbox_endTime.setText(str(newEndTime))
        self.drawPlot()

    def incrementForward(self):
        startTime = float(self.textbox_startTime.text())
        endTime = float(self.textbox_endTime.text())
        timeIncrement = float(self.textbox_timeIncrement.text())
        newStartTime = startTime + timeIncrement
        newEndTime = endTime + timeIncrement
        self.textbox_startTime.setText(str(newStartTime))
        self.textbox_endTime.setText(str(newEndTime))
        self.drawPlot()

    def incrementBack(self):
        startTime = float(self.textbox_startTime.text())
        endTime = float(self.textbox_endTime.text())
        timeIncrement = float(self.textbox_timeIncrement.text())
        newStartTime = startTime - timeIncrement
        newEndTime = endTime - timeIncrement
        self.textbox_startTime.setText(str(newStartTime))
        self.textbox_endTime.setText(str(newEndTime))
        self.drawPlot()

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

    def draw(self):
        self.fig.canvas.draw()

    def plotArray(self,image,normNSigma=3,title='',showColorBar=True,**kwargs):
        self.image = image
        if not 'vmax' in kwargs:
            goodImage = image[np.isfinite(image)]
            kwargs['vmax'] = np.mean(goodImage)+normNSigma*np.std(goodImage)
        if not 'cmap' in kwargs:
            defaultCmap=matplotlib.cm.hot
            defaultCmap.set_bad('0.15')
            kwargs['cmap'] = defaultCmap
        if not 'origin' in kwargs:
            kwargs['origin'] = 'lower'

        if 'button_press_event' in kwargs:
            cid = self.fig.canvas.mpl_connect('button_press_event',partial(kwargs.pop('button_press_event'),self))

        self.fig.clf()
        self.axes = self.fig.add_subplot(111)
        self.handleMatshow = self.axes.matshow(image,**kwargs)
        if showColorBar:
            self.fig.cbar = self.fig.colorbar(self.handleMatshow)
            cid = self.fig.canvas.mpl_connect('scroll_event', self.scrollColorBar)
            cid = self.fig.canvas.mpl_connect('button_press_event', self.clickColorBar)
        self.axes.set_title(title)
        cid = self.fig.canvas.mpl_connect('motion_notify_event', self.hoverCanvas)
        print 'draw'
        self.draw()

    def hoverCanvas(self,event):
        if event.inaxes is self.axes:
            col = int(round(event.xdata))
            row = int(round(event.ydata))
            if row < np.shape(self.image)[0] and col < np.shape(self.image)[1]:
                self.statusText.setText('({:d},{:d}) {}'.format(col,row,self.image[row,col]))


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

    def loadObsFile(self,obsPath):
        print 'loading',obsPath
        self.obs = ObsFile(obsPath)
        self.fn = FileName(obsFile=self.obs)
        firstImg = self.obs.getPixelCountImage(firstSec=0,integrationTime=1,getRawCount=True)
        self.image = firstImg['image']
        self.plotArray(self.image)
        self.updateHeader()
        self.label_obsFilename.setText(os.path.basename(obsPath))
        self.textbox_endTime.setText(str(self.headerInfo['exptime']))

    
    def loadObsFileWindow(self):
        LoadObsDialog(parent=self)

    def loadCalFiles(self):
        LoadCalsDialog(parent=self)

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

class HeaderWindow(ModelessWindow):
    def append(self,text):
        self.textArea.append(text)

    def setHeader(self,newHeaderText):
        self.textArea.setText(newHeaderText)

    def initUI(self):
        self.textArea = QtGui.QTextEdit('')
        mainBox = layoutBox('V',(self.textArea,))
        self.setLayout(mainBox)

class ImageParamsWindow(ModelessWindow):
    def initUI(self):
        self.checkbox_getRawCount = QtGui.QCheckBox('getRawCount',self)
        self.checkbox_weighted = QtGui.QCheckBox('weighted',self)
        self.checkbox_fluxWeighted = QtGui.QCheckBox('fluxWeighted',self)
        self.checkbox_scaleByEffInt = QtGui.QCheckBox('scaleByEffInt',self)

        mainBox = layoutBox('V',('ObsFile.getPixelCountImage Parameters',self.checkbox_getRawCount,
                    self.checkbox_weighted,self.checkbox_fluxWeighted,self.checkbox_scaleByEffInt))

        self.setLayout(mainBox)


    def getParams(self):
        obsParamsDict = {}
        plotParamsDict = {}
        obsParamsDict['getRawCount'] = self.checkbox_getRawCount.isChecked()
        obsParamsDict['weighted'] = self.checkbox_weighted.isChecked()
        obsParamsDict['fluxWeighted'] = self.checkbox_fluxWeighted.isChecked()
        obsParamsDict['scaleByEffInt'] = self.checkbox_scaleByEffInt.isChecked()
        outDict = {}
        outDict['obsParams'] = obsParamsDict
        outDict['plotParams'] = plotParamsDict
        return outDict
    

class LoadWidget(QtGui.QWidget):
    def __init__(self,parent=None,fileNameMethod='obs',promptStr='',initialParams={}):
        super(LoadWidget,self).__init__(parent=parent)
        self.parent=parent
        self.fileNameMethod = fileNameMethod
        self.promptStr=promptStr
        self.initialParams = initialParams
        self.initUI()

    def initUI(self):
        self.textbox_run = QtGui.QLineEdit(self.initialParams.get('run'))
        self.textbox_date = QtGui.QLineEdit(self.initialParams.get('date'))
        self.textbox_tstamp = QtGui.QLineEdit(self.initialParams.get('tstamp'))
        self.button_loadFilename = QtGui.QPushButton('Lookup FileName Path')

        self.button_openDialog = QtGui.QPushButton('Open Path')
        self.textbox_filename = QtGui.QLineEdit('')
        self.setFilenameFromParams()

        self.textbox_run.setMaximumWidth(200)
        self.textbox_date.setMaximumWidth(200)
        self.textbox_tstamp.setMaximumWidth(400)
        self.textbox_filename.setMaximumWidth(1000)
        self.button_openDialog.setMaximumWidth(100)

        fnBox = layoutBox('H',('run:',self.textbox_run,'date:',self.textbox_date,
                    'timestamp:',self.textbox_tstamp,self.button_loadFilename))
        dialogBox = layoutBox('H',(self.button_openDialog,self.textbox_filename))
        mainBox = layoutBox('V',(self.promptStr,fnBox,dialogBox))

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

        fn = FileName(run=run,date=date,tstamp=tstamp)
        path = getattr(fn,self.fileNameMethod)()
        self.textbox_filename.setText(path)

class LoadObsDialog(QtGui.QDialog):
    def __init__(self,parent=None):
        super(self.__class__,self).__init__(parent=parent)
        self.parent=parent
        self.initUI()
        self.connectButtons()
        self.exec_()

    def initUI(self):

        self.loadObsWidget = LoadWidget(self,'obs','Choose the Obs file')

        self.button_load = QtGui.QPushButton('Load Obs')
        self.button_cancel = QtGui.QPushButton('Cancel')
        buttonBox = layoutBox('H',(self.button_cancel,self.button_load))
        mainBox = layoutBox('V',(self.loadObsWidget,buttonBox))
        self.setLayout(mainBox)

    def connectButtons(self):
        self.connect(self.button_load,QtCore.SIGNAL('clicked()'),self.load)
        self.connect(self.button_cancel,QtCore.SIGNAL('clicked()'),self.close)

    def load(self):
        filename = str(self.loadObsWidget.textbox_filename.text())
        self.parent.loadObsFile(filename)
        self.close()


class LoadCalsDialog(QtGui.QDialog):
    def __init__(self,parent=None):
        super(self.__class__,self).__init__(parent=parent)
        self.parent=parent
        self.initUI()
        self.connectButtons()
        self.exec_()

    def initUI(self):

        self.loadWvlWidget = LoadWidget(self,'calSoln','Choose the wavelength cal')
        self.loadFlatWidget = LoadWidget(self,'flatSoln','Choose the flat cal')
        self.loadFluxWidget = LoadWidget(self,'fluxSoln','Choose the flux cal')
        self.loadTimeMaskWidget = LoadWidget(self,'timeMask','Choose the time mask file')
        self.loadTimeAdjustmentWidget = LoadWidget(self,'timeAdjustments','Choose the time adjustment file')
        self.loadCosmicWidget = LoadWidget(self,'cosmicMask','Choose the cosmic mask file')
        self.loadBeammapWidget = LoadWidget(self,'beammap','Choose the beammap file')

        self.button_loadWvl = QtGui.QPushButton('Load Wvl Cal')
        self.button_loadFlat = QtGui.QPushButton('Load Flat Cal')
        self.button_loadFlux = QtGui.QPushButton('Load Flux Cal')
        self.button_loadTimeMask = QtGui.QPushButton('Load Time Mask')
        self.button_loadTimeAdjustment = QtGui.QPushButton('Load Time Adjustments')
        self.button_loadCosmic = QtGui.QPushButton('Load Cosmic Mask')
        self.button_loadBeammap = QtGui.QPushButton('Load Beammap')

        self.button_load = QtGui.QPushButton('Load All Cals')

        wvlBox = layoutBox('H',(self.loadWvlWidget,self.button_loadWvl))
        flatBox = layoutBox('H',(self.loadFlatWidget,self.button_loadFlat))
        fluxBox = layoutBox('H',(self.loadFluxWidget,self.button_loadFlux))
        timeMaskBox = layoutBox('H',(self.loadTimeMaskWidget,self.button_loadTimeMask))
        timeAdjustBox = layoutBox('H',(self.loadTimeAdjustmentWidget,self.button_loadTimeAdjustment))
        cosmicBox = layoutBox('H',(self.loadCosmicWidget,self.button_loadCosmic))
        beammapBox = layoutBox('H',(self.loadBeammapWidget,self.button_loadBeammap))

        mainBox = layoutBox('V',(wvlBox,flatBox,fluxBox,timeMaskBox,timeAdjustBox,cosmicBox,beammapBox,self.button_load))
        self.setLayout(mainBox)

    def connectButtons(self):
        self.connect(self.button_load,QtCore.SIGNAL('clicked()'),self.load)

    def load(self):
        filename = str(self.loadObsWidget.textbox_filename.text())
        self.close()
#QtGui.QFileDialog.getOpenFileName(self, 'Open File',".","(*.writer)")
        

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
                    box.addSpacing(element)
                except:
                    try:
                        box.addWidget(QtGui.QLabel(element))
                    except:
                        print 'could\'t add {} to layout box'.format(element)
    return box

if __name__ == "__main__":
    form = ObsFileViewer()



