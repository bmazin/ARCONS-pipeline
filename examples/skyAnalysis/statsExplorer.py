from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import numpy as np
import sys
import os

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

from util.ObsFile import ObsFile
from util.FileName import FileName
from util.readDict import readDict
from util.popup import PopUp
from hotpix import hotPixelsMatt as hotPixels

def mad(a,axis=-1):
    if axis==-1:
        return 1.4826*np.median(np.abs(a-np.median(a)))
    elif axis == 0:
        return 1.4826*np.median(np.abs(a-np.median(a)))
        
class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Pixel Explorer')

        paramFile = sys.argv[1]
        self.params = readDict()
        self.params.read_from_file(paramFile)
        self.createMainFrame()
        self.createStatusBar()

    def createMainFrame(self):
        self.main_frame = QWidget()
      # Create the mpl Figure and FigCanvas objects. 
        self.dpi = 100
        self.fig = Figure((7.0, 7.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        self.axes0 = self.fig.add_subplot(111)
        cid=self.canvas.mpl_connect('button_press_event', self.clickCanvas)

        # Create the navigation toolbar, tied to the canvas
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def createStatusBar(self):
        self.status_text = QLabel("Awaiting orders.")
        self.statusBar().addWidget(self.status_text, 1)

    def openImage(self):

        timestampList = [self.params['obsUtcDate']+'-'+ts for ts in self.params['obsSequence']] 
        run = self.params['run']
        sunsetDate = self.params['obsSunsetDate']
        utcDate = self.params['obsUtcDate']
        self.intTime = self.params['intTime']
        wvlLowerCutoff = self.params['wvlLowerCutoff']
        wvlUpperCutoff = self.params['wvlUpperCutoff']

        calTimestamp = self.params['wvlTimestamp']
        wfn = FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln()
        calfn = FileName(run=run,date=self.params['wvlSunsetDate'],tstamp=calTimestamp).cal()
        ffn = FileName(run=run,date=self.params['flatCalSunsetDate'],tstamp='').flatSoln()

        obsFns = [FileName(run=run,date=sunsetDate,tstamp=timestamp).obs() for timestamp in timestampList]
        self.obList = [ObsFile(obsFn) for obsFn in obsFns]
        for ob in self.obList:
            print 'Loading ',ob.fullFileName
            ob.loadWvlCalFile(wfn)
            ob.loadFlatCalFile(ffn)

        self.cal = ObsFile(calfn)
        self.cal.loadWvlCalFile(wfn)
        self.cal.loadFlatCalFile(ffn)
        self.loadSpectra()

    def loadSpectra(self):
        fileName = self.params['outFileName']
        if os.path.exists(fileName):
            data = np.load(fileName)
            self.spectra = data['cube']
            self.frame = data['frame']
            self.frameValues = self.frame[~np.isnan(self.frame)]
            self.wvlBinEdges = np.array(self.obList[0].flatCalWvlBins)
            self.frameIntTime = 300*6
        else:
            self.spectra,self.wvlBinEdges = self.obList[0].getSpectralCube(weighted=True)
            self.frameIntTime = 0
            for ob in self.obList[1:]:
                print ob.fileName
                cube,wvlBinEdges = ob.getSpectralCube(weighted=True)
                self.spectra += cube
                self.frameIntTime += ob.getFromHeader('exptime')
            self.spectra = np.array(self.spectra,dtype=np.float64)
            self.frame = np.sum(self.spectra,axis=2)
            hotPixMask = hotPixels.findHotPixels(image=self.frame,nsigma=2)['badflag']
            self.frame[hotPixMask != 0] = np.nan
            self.frame[self.frame == 0] = np.nan
            self.frameValues = self.frame[~np.isnan(self.frame)]
            np.savez(fileName,cube=self.spectra,frame=self.frame)

    def plotWeightedImage(self):
        self.showFrame = np.array(self.frame)
        self.showFrame[np.isnan(self.frame)] = 0
        handleMatshow = self.axes0.matshow(self.showFrame,cmap=matplotlib.cm.gnuplot2,origin='lower',vmax=np.mean(self.showFrame)+3*np.std(self.showFrame))
        self.fig.colorbar(handleMatshow)

    def clickCanvas(self,event):
        self.showLaserSpectrum = True
        self.showPixelSpectrum = True
        self.showWvlLightCurves = True
        self.showWvlLightCurveHists = False
        self.showStdVsIntTime = True
        self.showNormStdVsIntTime = True

        col = round(event.xdata)
        row = round(event.ydata)
        
        
        if self.showPixelSpectrum:
            #next plot the integrated spectrum for this pixel in the total image
            spectrum = self.spectra[row,col]
            print sum(spectrum),' counts in broadband spectrum'
            def plotFunc(fig,axes):
                axes.plot(self.wvlBinEdges[:-1],spectrum)
                axes.set_xlabel(r'$\lambda$ ($\AA$)')
                axes.set_ylabel(r'total counts')
            popup = PopUp(parent=self,plotFunc=plotFunc,title='spectrum, pixel %d,%d (intTime=%d)'%(row,col,self.frameIntTime))

        rebinSpecBins = 5
        firstAfterConvolve = rebinSpecBins//2
        rebinnedWvlEdges = self.wvlBinEdges[::rebinSpecBins]

        if self.showLaserSpectrum:
            #First plot the laser cal spectrum for this pixel to see if it's good
            laserSpectrum,binEdges = self.cal.getPixelSpectrum(row,col,weighted=True)
            def plotFunc(fig,axes):
                axes.plot(binEdges[:-1],laserSpectrum)
                axes.set_xlabel(r'$\lambda$ ($\AA$)')
                axes.set_ylabel(r'total counts')
            popup = PopUp(parent=self,plotFunc=plotFunc,title='Laser Cal Spectrum, pixel %d,%d'%(row,col))
            
        if self.showWvlLightCurves:
            spectrumInTime = []
            for iOb,ob in enumerate(self.obList):
                for sec in range(0,ob.getFromHeader('exptime'),self.intTime):
                    spectrum,binEdges = ob.getPixelSpectrum(pixelRow=row,pixelCol=col,firstSec=sec,integrationTime=self.intTime,weighted=True)
                    spectrum = np.convolve(spectrum,np.ones(rebinSpecBins),'same')[firstAfterConvolve::rebinSpecBins]
                    spectrumInTime.append(spectrum)
            spectrumInTime = np.array(spectrumInTime)
            nBins = np.shape(spectrumInTime)[1]
            def plotFunc(fig,axes):
                #plot counts vs time for each wavelength bin
                t=np.arange(len(spectrumInTime[:,0]))*self.intTime
                for iBin in xrange(nBins):
                    axes.plot(t,1.0*spectrumInTime[:,iBin]/self.intTime,
                        c=cm.jet((iBin+1.)/nBins),
                        label=r'%d-%d $\AA$'%(rebinnedWvlEdges[iBin],
                        rebinnedWvlEdges[iBin+1]))
                    axes.set_xlabel('time (s)')
                    axes.set_ylabel('cps')
                #plot counts vs time summed over all wavelengths
                axes.plot(t,np.sum(spectrumInTime,axis=1)/self.intTime,c='k',
                    label=r'%d-%d $\AA$'%(rebinnedWvlEdges[0],rebinnedWvlEdges[-1]))
                #axes.legend(loc='center right')
            popup = PopUp(parent=self,plotFunc=plotFunc,title='Light Curve by Band, Pixel %d,%d'%(row,col))

        if self.showWvlLightCurveHists or self.showStdVsIntTime or self.showNormStdVsIntTime:
            intTimes = [1,2,3,5,10,15,30]
            spectrumVsIntTimeVsTime = []
            for intTime in intTimes:
                spectrumInTime = []
                for iOb,ob in enumerate(self.obList):
                    for sec in range(0,ob.getFromHeader('exptime'),intTime):
                        spectrum,binEdges = ob.getPixelSpectrum(pixelRow=row,pixelCol=col,firstSec=sec,integrationTime=intTime,weighted=True)
                        spectrum = np.convolve(spectrum,np.ones(rebinSpecBins),'same')[firstAfterConvolve::rebinSpecBins]
                        
                        spectrumInTime.append(spectrum)
                spectrumInTime = np.array(spectrumInTime)
                spectrumVsIntTimeVsTime.append(spectrumInTime)
            #resulting array indexed as
            #spectrumVsIntTimeVsTime[iIntTime][iTimeChunk][iWvlBin]

            #sum over wavelength for total counts
            countsVsIntTimeVsTime = [np.sum(spectrumInTime,axis=1) for spectrumInTime in spectrumVsIntTimeVsTime]
            #countsVsIntTimeVsTime[iIntTime][iTimeChunk]
            if self.showWvlLightCurveHists:
                for iIntTime,countsVsTime in enumerate(countsVsIntTimeVsTime):
                    def plotFunc(fig,axes):
                        axes.hist(countsVsTime,bins=20)
                        axes.set_xlabel('counts per intTime %d s'%intTimes[iIntTime])
                    popup = PopUp(parent=self,plotFunc=plotFunc,title='Int Time %d s, pixel %d,%d'%(intTimes[iIntTime],row,col))

            countStds = [np.std(countsVsTime) for countsVsTime in countsVsIntTimeVsTime]
            countStds = np.array(countStds)
            countSqrts = [np.sqrt(np.median(countsVsTime)) for countsVsTime in countsVsIntTimeVsTime]
            countSqrts = np.array(countSqrts)
            spectrumStds = [np.std(spectrumVsTime,axis=0) for spectrumVsTime in spectrumVsIntTimeVsTime]
            spectrumSqrts = [np.sqrt(np.median(spectrumVsTime,axis=0)) for spectrumVsTime in spectrumVsIntTimeVsTime]
            spectrumStds = np.array(spectrumStds)
            spectrumSqrts = np.array(spectrumSqrts)

            if self.showStdVsIntTime:
                def plotFunc(fig,axes):
                    axes.set_xlabel('integration time (s)')
                    axes.plot(intTimes,countStds,'k',label=r'total $\sigma$')
                    axes.plot(intTimes,countSqrts,'k--',label=r'$\sqrt{med(N)}$')
                    nBins = np.shape(spectrumStds)[1]
                    for iBin in xrange(nBins):
                        axes.plot(intTimes,spectrumStds[:,iBin],
                            c=cm.jet((iBin+1.)/nBins),
                            label=r'%d-%d $\AA$ $\sigma$'%(rebinnedWvlEdges[iBin],
                            rebinnedWvlEdges[iBin+1]))
                        axes.plot(intTimes,spectrumSqrts[:,iBin],
                            c=cm.jet((iBin+1.)/nBins),linestyle='--')
                    axes.legend(loc='upper left')

                popup = PopUp(parent=self,plotFunc=plotFunc,
                    title=r'$\sigma$ vs Integration Time, Pixel %d,%d'%(row,col))

                
            if self.showNormStdVsIntTime:
                def plotFunc(fig,axes):
                    axes.set_xlabel('integration time (s)')
                    axes.set_ylabel('normalized $\sigma$')
                    axes.plot(intTimes,countSqrts/np.max(countSqrts),'k--',
                        label=r'$\sqrt{N}$')
                    axes.plot(intTimes,countStds/np.max(countSqrts),'k',
                        label=r'%d-%d $\AA$'%(rebinnedWvlEdges[0],rebinnedWvlEdges[-1]))
                    nBins = np.shape(spectrumStds)[1]
                    for iBin in xrange(nBins):
                        axes.plot(intTimes,
                            spectrumStds[:,iBin]/np.max(spectrumSqrts[:,iBin]),
                            c=cm.jet((iBin+1.0)/nBins),
                            label=r'%d-%d $\AA$'%(rebinnedWvlEdges[iBin],
                                rebinnedWvlEdges[iBin+1]))
                    axes.legend(loc='upper left')
                popup = PopUp(parent=self,plotFunc=plotFunc,
                    title='Normalized Standard Deviation vs IntegrationTime, Pixel %d,%d'%(row,col))
            
        

        
def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.openImage()
    form.plotWeightedImage()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
