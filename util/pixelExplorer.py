from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4 import QtGui
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import numpy as np
import sys
import os
import tables
from scipy.stats import chi2

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
from functools import partial

from util.ObsFile import ObsFile
from util.FileName import FileName
from util.readDict import readDict
from util.popup import PopUp,onscroll_cbar,onclick_cbar
#from hotpix import hotPixelsMatt as hotPixels
from hotpix import hotPixels

class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('Pixel Explorer')

        paramFile = sys.argv[1]
        self.params = readDict()
        self.params.read_from_file(paramFile)

        self.createMainFrame()
        self.createStatusBar()

    def __del__(self):
        for stackLabel in self.stackObsFileLists:
            for obs in self.stackObsFileLists[stackLabel]:
                try:
                    obs.close()
                except:
                    pass
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

    def loadStacks(self):
        #make a dictionary to put lists of ObsFile objects into
        self.stackObsFileLists = dict()
        self.stackWvlCals = dict()

        for stackLabel in ['obs','sky','twilight']:
            self.loadObs(stackLabel)
        #Unload the lists of ObsFile objects into individual lists
        self.obList = self.stackObsFileLists['obs']
        self.skyObList = self.stackObsFileLists['sky']
        self.twilightObList = self.stackObsFileLists['twilight']
        self.cal = self.stackWvlCals['obs']

        #If there isn't already an npz holding the obs stack spectra 
        #or if the params say to remake the npz regardless, then make the stack npz file
        if os.path.exists(self.params['obsStackFileName'])==False or self.params['makeObsStackFile'] == True:
            self.createStack('obs')
        if os.path.exists(self.params['rawStackFileName'])==False or self.params['makeRawObsStackFile'] == True:
            self.createStack('raw')
        if os.path.exists(self.params['skyStackFileName'])==False or self.params['makeSkyStackFile'] == True:
            self.createStack('sky')
        if os.path.exists(self.params['twilightStackFileName'])==False or self.params['makeTwilightStackFile'] == True:
            self.createStack('twilight')

        #Now load all the info from the files, whether made previously or just now
        stackFileName = self.params['obsStackFileName']
        data = np.load(stackFileName)
        self.spectra = data['spectra']
        self.frame = data['frame']
        #self.obsTotalIntTime = data['totalIntTime']
        self.wvlBinEdges = data['wvlBinEdges']
        self.wvlBinWidths = np.diff(self.wvlBinEdges)
        self.nRow,self.nCol,self.nWvlBins = np.shape(self.spectra)
        self.intTime = self.params['obsIntTime']
        
        stackFileName = self.params['rawStackFileName']
        data = np.load(stackFileName)
        self.rawSpectra = data['spectra']
        self.rawFrame = data['frame']
        #self.rawTotalIntTime = data['totalIntTime']

        stackFileName = self.params['skyStackFileName']
        data = np.load(stackFileName)
        self.skySpectra = data['spectra']
        self.skyFrame = data['frame']
        #self.skyTotalIntTime = data['totalIntTime']
        self.skyIntTime = self.params['skyIntTime']
        
        stackFileName = self.params['twilightStackFileName']
        data = np.load(stackFileName)
        self.twilightSpectra = data['spectra']
        self.twilightFrame = data['frame']
        #self.twilightTotalIntTime = data['totalIntTime']
        self.twilightIntTime = self.params['twilightIntTime']

        
        #Load flat info for the current obs stack
        flatSolnFileName = FileName(run=self.params['run'],date=self.params['obsFlatCalSunsetDate'],tstamp=self.params['obsFlatCalTimestamp']).flatSoln()
        flatInfoFileName =  os.path.splitext(flatSolnFileName)[0]+'.npz' 
        self.flatInfo= np.load(flatInfoFileName)

        #Load flat info for all flats specified in params file
        flats = self.params['flatInfoFiles']
        self.flatInfos = []
        for flat in flats:
            flatSolnFileName = FileName(run=self.params['run'],date=flat,tstamp='').flatSoln()
            flatInfoFileName =  os.path.splitext(flatSolnFileName)[0]+'.npz' 
            self.flatInfos.append(np.load(flatInfoFileName))



    def prepareForClickPlots(self):
        matplotlib.rcParams['font.size'] = 10
        #create wavelength array more coarsely binned for use in showPixelWvlLightCurves and similar
        self.rebinSpecBins = self.params['nWvlBands']
        self.firstAfterConvolve = self.rebinSpecBins//2
        self.rebinnedWvlEdges = self.wvlBinEdges[::self.rebinSpecBins]
        self.averageTwilightSpectrum = np.zeros(self.nWvlBins)
        spectra2d = np.reshape(self.twilightSpectra,[self.nRow*self.nCol,self.nWvlBins ])
        fractionOfPixelsToTrim = .1
        for iWvl in xrange(self.nWvlBins):
            spectrum = spectra2d[:,iWvl]
            goodSpectrum = spectrum[spectrum != 0]#dead pixels need to be taken out before calculating medians
            goodSpectrum = np.sort(goodSpectrum)
            nGoodPixels = len(goodSpectrum)
            trimmedSpectrum = goodSpectrum[fractionOfPixelsToTrim*nGoodPixels:(1-fractionOfPixelsToTrim)*nGoodPixels]
            self.averageTwilightSpectrum[iWvl] = np.mean(trimmedSpectrum)

    def loadObs(self,stackLabel):
        timestampList = self.params[stackLabel+'Sequence']
        run = self.params['run']
        sunsetDate = self.params[stackLabel+'SunsetDate']
        utcDate = self.params[stackLabel+'UtcDate']
        intTime = self.params[stackLabel+'IntTime']
        wvlLowerCutoff = self.params[stackLabel+'WvlLowerCutoff']
        wvlUpperCutoff = self.params[stackLabel+'WvlUpperCutoff']

        calTimestamp = self.params[stackLabel+'WvlTimestamp']
        print stackLabel,calTimestamp
        wvlSolnFileName = FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln()
        wvlCalFileName = FileName(run=run,date=self.params[stackLabel+'WvlSunsetDate'],tstamp=calTimestamp).cal()
        flatSolnFileName = FileName(run=run,date=self.params[stackLabel+'FlatCalSunsetDate'],tstamp=self.params[stackLabel+'FlatCalTimestamp']).flatSoln()

        obsFileNames = [FileName(run=run,date=sunsetDate,tstamp=timestamp).obs() for timestamp in timestampList]
        obList = [ObsFile(obsFn) for obsFn in obsFileNames]
        for ob in obList:
            ob.loadWvlCalFile(wvlSolnFileName)
            ob.loadFlatCalFile(flatSolnFileName)

        self.stackObsFileLists[stackLabel] = obList
        
        cal = ObsFile(wvlCalFileName)
        cal.loadWvlCalFile(wvlSolnFileName)
        cal.loadFlatCalFile(flatSolnFileName)
        self.stackWvlCals[stackLabel] = cal
        
    def createStack(self,stackLabel):

        paramsLabel = stackLabel
        weighted = True
        if stackLabel == 'raw':
            paramsLabel = 'obs'
            weighted = False
            getRawCounts = True
        x = self.stackObsFileLists[paramsLabel][0].getSpectralCube(weighted=weighted)
        spectra = x['cube']
        wvlBinEdges = x['wvlBinEdges']
        totalIntTime = 0
        for ob in self.stackObsFileLists[paramsLabel][1:]:
            print ob.fileName
            x = ob.getSpectralCube(weighted=weighted)
            cube = x['cube']
            wvlBinEdges = x['wvlBinEdges']
            spectra += cube
            totalIntTime += ob.getFromHeader('exptime')
        spectra = np.array(spectra,dtype=np.float64)
        frame = np.sum(spectra,axis=2)
        hotPixMask = hotPixels.checkInterval(image=frame)['mask']
        frame[hotPixMask != 0] = np.nan
        frame[frame == 0] = np.nan
        stackFileName = self.params[stackLabel+'StackFileName']
        np.savez(stackFileName,spectra=spectra,frame=frame,wvlBinEdges=wvlBinEdges,totalIntTime=totalIntTime)

    def plotWeightedImage(self):
        self.showFrame = np.array(self.frame)
        self.showFrame[np.isnan(self.frame)] = 0
        handleMatshow = self.axes0.matshow(self.showFrame,cmap=matplotlib.cm.gnuplot2,origin='lower',vmax=np.mean(self.showFrame)+3*np.std(self.showFrame))
        self.fig.cbar = self.fig.colorbar(handleMatshow)
        cid = self.fig.canvas.mpl_connect('scroll_event', partial(onscroll_cbar, self.fig))
        cid = self.fig.canvas.mpl_connect('button_press_event', partial(onclick_cbar, self.fig))

    def arrayPlots(self):
        if self.params['showArrayRawImage']:
            self.showArrayRawImage()
        if self.params['showArrayStdVsIntTime']:
            self.showArrayStdVsIntTime()
        if self.params['showArrayWvlCalRange']:
            self.showArrayWvlCalRange()
        if self.params['showTwilightArrayImage']:
            self.showTwilightArrayImage()
        if self.params['showTwilightArrayStdVsFlux']:
            self.showTwilightArrayStdVsFlux()
        if self.params['showTwilightArrayReducedChisqImage']:
            self.showTwilightArrayReducedChisqImage()
        if self.params['showSkyArrayImage']:
            self.showSkyArrayImage()
        if self.params['showArrayLaserImage']:
            self.showArrayLaserImage()

    def clickCanvas(self,event):
        if event.inaxes is self.axes0:
            col = round(event.xdata)
            row = round(event.ydata)
            print '(%d,%d)'%(row,col)
            
            if self.params['showPixelSpectrum'] or self.params['showPixelRawSpectrum']:
                self.showPixelSpectrum(row,col)
            if self.params['showPixelLightCurve']:
                self.showPixelLightCurve(row,col)
            if self.params['showPixelWvlLightCurves']:
                self.showPixelWvlLightCurves(row,col)
            if self.params['showPixelRawPhaseHist']:
                self.showPixelRawPhaseHist(row,col)
            if self.params['showPixelRawBaselineHist']:
                self.showPixelRawBaselineHist(row,col)
            if self.params['showPixelRawPeakHist']:
                self.showPixelRawPeakHist(row,col)

            if self.params['showPixelFlatWeights']:
                self.showPixelFlatWeights(row,col)

            if self.params['showPixelLaserSpectrum']:
                self.showPixelLaserSpectrum(row,col)

            if self.params['showTwilightPixelSpectrum']:
                self.showTwilightPixelSpectrum(row,col)
            if self.params['showTwilightPixelStdVsFlux']:
                self.showTwilightPixelStdVsFlux(row,col)
            if self.params['showTwilightPixelDeviationFromMedian']:
                self.showTwilightPixelDeviationFromMedian(row,col)
            
            if self.params['showPixelStdVsIntTime']:
                self.showPixelStdVsIntTime(row,col)

    def showArrayRawImage(self):
        output = self.obList[0].getPixelCountImage(getRawCount=True,weighted=False)
        frame = output['image']
        for ob in self.obList[1:]:
            output = ob.getPixelCountImage(getRawCount=True,weighted=False)
            frame += output['image']
        PopUp(parent=self,title='showArrayRawImage').plotArray(image=frame,title='raw image')

    def showArrayStdVsIntTime(self):
        intTimes = [1,2,3,5,10,15,30]
        sdevVsIntTime = []
        madVsIntTime = []
        medVsIntTime = []
        for intTime in intTimes:
            image = np.zeros((self.nRow,self.nCol))
            for iOb,ob in enumerate(self.skyObList):
                x = ob.getPixelCountImage(firstSec=0,integrationTime=intTime)
                image+=x['image']
            hotPixMask = hotPixels.checkInterval(image=image)['mask']
            image[hotPixMask!=0]=0
            countList = image[image!=0]
            sdevVsIntTime.append(np.std(countList))
            madVsIntTime.append(np.median(np.abs(countList-np.median(countList))))
            medVsIntTime.append(np.median(countList))
            PopUp(parent=self).plotArray(image,title=r'%d std=%f mad=%f med=%f'%(intTime,sdevVsIntTime[-1],madVsIntTime[-1],medVsIntTime[-1]))
        medVsIntTime = np.array(medVsIntTime)
        sqrtNVsIntTime = np.sqrt(medVsIntTime)
        pop = PopUp(parent=self,title='showArrayStdVsIntTime')
        pop.axes.set_xlabel('integration time (s)')
        pop.axes.set_ylabel('$\sigma$')
        pop.axes.plot(intTimes,sqrtNVsIntTime,'k--',label=r'$\sqrt(med(N))$')
        pop.axes.plot(intTimes,sdevVsIntTime,'k')
        pop.axes.plot(intTimes,madVsIntTime,'r')
        pop.draw()

    def showTwilightArrayImage(self):
        image = self.twilightFrame
        image[np.isnan(image)] = 0
        #self.popUpArray(image=image,title='twilight image')
        PopUp(parent=self,title='showTwilightArrayImage').plotArray(image=image,title='Twilight Image')
        
    def showTwilightArrayStdVsFlux(self):
        pass

    def showTwilightArrayReducedChisqImage(self):
        chisqImage = np.zeros((self.nRow,self.nCol))
        nDeltaFromZero = np.zeros((self.nRow,self.nCol,self.nWvlBins))
        for iRow in range(self.nRow):
            for iCol in range(self.nCol):
                x = self.getChisq(iRow,iCol)
                chisqImage[iRow,iCol] = x['reducedChisq']
                nDeltaFromZero[iRow,iCol,:] = x['nDeltaFromZero']
        chisqImage[np.isnan(chisqImage)]=0
        chisqImage[chisqImage == np.inf]=0
        nDeltaFromZero = np.ma.array(nDeltaFromZero,mask=np.logical_and(np.isnan(nDeltaFromZero),nDeltaFromZero==np.inf))
        hotPixMask = hotPixels.checkInterval(image=chisqImage)['mask']
        chisqImage[hotPixMask != 0] = 0
        #self.popUpArray(image=chisqImage,title='Flat Cal $\chi^{2}_{red}$',normNSigma=1.)
        PopUp(parent=self,title='showTwilightArrayReducedChisqImage').plotArray(image=chisqImage,title='Flat Cal $\chi^{2}_{red}$',normNSigma=1.)
        #PopUp(parent=self,title='showTwilightArrayNDeltaFromZero').plotArray(image=nDeltaFromZero,title='Flat Cal n$\sigma$ from 0',normNSigma=3.)
#        verbose = True
#        pdfFullPath='/Scratch/flatCalSolnFiles2/nDeltaWvlSlicesTwiAppliedToSky.pdf'
#        pp = PdfPages(pdfFullPath)
#        nPlotsPerRow = 3
#        nPlotsPerCol = 4
#        nPlotsPerPage = nPlotsPerRow*nPlotsPerCol
#        iPlot = 0
#        if verbose:
#            print 'plotting weights in wavelength sliced images'
#
#        #matplotlib.rcParams['font.size'] = 4
#        wvls = self.wvlBinEdges[0:-1]
#
#        for iWvl,wvl in enumerate(wvls):
#            if verbose:
#                print 'wvl ',iWvl
#            if iPlot % nPlotsPerPage == 0:
#                fig = plt.figure(figsize=(10,10),dpi=100)
#
#            ax = fig.add_subplot(nPlotsPerCol,nPlotsPerRow,iPlot%nPlotsPerPage+1)
#            ax.set_title(r'%.0f $\AA$'%wvl)
#            image = nDeltaFromZero[:,:,iWvl]
#
#            cmap = matplotlib.cm.gnuplot2
#            cmap.set_bad('.1')
#            
#            handleMatshow = ax.matshow(image,cmap=cmap,origin='lower',vmin=-10.,vmax=10.)
#            cbar = fig.colorbar(handleMatshow)
#
#            if iPlot%nPlotsPerPage == nPlotsPerPage-1:
#                pp.savefig(fig)
#            iPlot += 1
#        pp.savefig(fig)
#        pp.close()


    def showSkyArrayImage(self):
        image = self.skyFrame
        image[np.isnan(image)] = 0
        #self.popUpArray(image,title='sky image')
        PopUp(parent=self,title='showSkyArrayImage').plotArray(image=image,title='Sky Image')

    def showArrayWvlCalRange(self):
        rangeTable = self.cal.wvlRangeTable[:,:,1]#Get upper limit from valid wvl ranges 
        #self.popUpArray(image=rangeTable,title=r'Upper Wavecal Limits ($\AA$)')
        PopUp(parent=self,title='showArrayWvlCalRange').plotArray(image=rangeTable,title=r'Upper Wavecal Limits ($\AA$)')


    def showPixelSpectrum(self,row,col):
        spectrum = self.spectra[row,col]
        binWidths = np.diff(self.wvlBinEdges)
        if self.params['showPixelRawSpectrum']:
            weights = self.flatInfo['weights'][row,col]
            rawSpectrum = self.spectra[row,col]/weights
            rawSpectrum/=binWidths

        spectrum/=binWidths
        pop = PopUp(parent=self,title='showPixelSpectrum')
        pop.axes.step(self.wvlBinEdges[:-1],spectrum,label='calibrated',color='b',where='post')
        if self.params['showPixelRawSpectrum']:
            pop.axes.step(self.wvlBinEdges[:-1],rawSpectrum,label='raw',color='r',where='post')
        pop.axes.set_xlabel(r'$\lambda$ ($\AA$)')
        pop.axes.set_ylabel(r'counts/$\AA$')
        pop.axes.legend(loc='lower right')
        pop.axes.set_title('spectrum (%d,%d)'%(row,col))
        pop.draw()

    def showPixelWvlLightCurves(self,row,col):
        spectrumInTime = []
        for iOb,ob in enumerate(self.obList):
            for sec in range(0,ob.getFromHeader('exptime'),self.intTime):
                x = ob.getPixelSpectrum(pixelRow=row,pixelCol=col,firstSec=sec,integrationTime=self.intTime,weighted=True)
                spectrum = x['spectrum']
                spectrum = np.convolve(spectrum,np.ones(self.rebinSpecBins),'same')[self.firstAfterConvolve::self.rebinSpecBins]
                spectrumInTime.append(spectrum)
        spectrumInTime = np.array(spectrumInTime)
        nBins = np.shape(spectrumInTime)[1]
        pop = PopUp(parent=self,title='showPixelWvlLightCurves')
        #plot counts vs time for each wavelength bin
        times=np.arange(len(spectrumInTime[:,0]))*self.intTime
        for iBin in xrange(nBins):
            pop.axes.plot(times,1.0*spectrumInTime[:,iBin]/self.intTime,
                c=cm.jet((iBin+1.)/nBins),
                label=r'%d-%d $\AA$'%(self.rebinnedWvlEdges[iBin],
                self.rebinnedWvlEdges[iBin+1]))
        pop.axes.set_xlabel('time (s)')
        pop.axes.set_ylabel('cps')
        #plot counts vs time summed over all wavelengths
        pop.axes.legend(loc='upper right')
        pop.axes.set_title('Light Curve by Band (%d,%d)'%(row,col))
        pop.draw()


    def showPixelLightCurve(self,row,col):
        lightCurve = []
        for iOb,ob in enumerate(self.obList):
            for sec in range(0,ob.getFromHeader('exptime'),self.intTime):
                x = ob.getPixelCount(iRow=row,iCol=col,firstSec=sec,integrationTime=self.intTime,weighted=True)
                counts = x['counts']/self.intTime
                lightCurve.append(counts)
        pop = PopUp(parent=self,title='showPixelLightCurve')
        times=np.arange(0,len(lightCurve)*self.intTime,self.intTime)
        pop.axes.set_xlabel('time (s)')
        pop.axes.set_ylabel('cps')
        pop.axes.plot(times,lightCurve,c='k')
        pop.axes.set_title('Light Curve (%d,%d)'%(row,col))
        pop.draw()

    def showArrayLaserImage(self):
        getImageOutput = self.cal.getPixelCountImage(getRawCount=True,weighted=False)
        frame = getImageOutput['image']
        #self.popUpArray(image=self.rawFrame,title='raw image')
        pop = PopUp(parent=self,title='showArrayLaserImage')
        pop.plotArray(image=frame,title='laser cal raw image')
        
    def showPixelLaserSpectrum(self,row,col):
        #First plot the laser cal spectrum for this pixel to see if it's good
        x = self.cal.getTimedPacketList(row,col)
        phases=np.array(x['peakHeights'],dtype=np.double)-np.array(x['baselines'],dtype=np.double)
        pop = PopUp(parent=self,title='showPixelLaserSpectrum')
        nBins=np.max(phases)-np.min(phases)
        histPhases,binEdges = np.histogram(phases,bins=nBins)
        lambdaBinEdges = self.cal.convertToWvl(binEdges,row,col,excludeBad=True)
        pop.axes.set_xlabel(r'$\lambda$ ($\AA$)')
        if len(lambdaBinEdges)==0: #no wavecal for this pixel, so lambdaBinEdges came back empty
            lambdaBinEdges = binEdges
            pop.axes.set_xlabel('phase (ADU)')
        pop.axes.step(lambdaBinEdges[:-1],histPhases,where='post',color='k')

        pop.axes.set_ylabel('counts')
        pop.axes.set_title('Raw Laser Cal Spectrum (%d,%d)'%(row,col))
        wvlCalSigma = self.cal.wvlErrorTable[row,col]
        xOffset = self.cal.wvlCalTable[row,col,0]
        yOffset = self.cal.wvlCalTable[row,col,1]
        amplitude = self.cal.wvlCalTable[row,col,2]
        #energy(eV) = amplitude*(pulseHeight-xOffset)**2+yOffset
       
        stackLabel = 'obs'
        run = self.params['run']
        sunsetDate = self.params[stackLabel+'SunsetDate']
        calTimestamp = self.params[stackLabel+'WvlTimestamp']
        wvlDriftFileName = FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calDriftInfo()
        wvlDriftFile = tables.openFile(wvlDriftFileName,mode='r')
        wvlDriftInfo = wvlDriftFile.root.params_drift.driftparams.read()
        wvlDriftFile.close()
        driftEntry = wvlDriftInfo[np.logical_and(wvlDriftInfo['pixelrow']==row ,wvlDriftInfo['pixelcol']==col)][0]
        #extract gaussianparams in first column of selected row
        bluePhaseSigma=driftEntry[0][2]
        bluePhaseAmp = driftEntry[0][1]
        bluePhaseOffset = driftEntry[0][0]

        redPhaseSigma=driftEntry[0][5]
        redPhaseAmp = driftEntry[0][4]
        redPhaseOffset = driftEntry[0][3]

        phases = np.linspace(np.min(phases),np.max(phases),(nBins+1)*100.)
        blueGaussFit = bluePhaseAmp*np.exp(-1/2*((phases-bluePhaseOffset)/bluePhaseSigma)**2)
        redGaussFit = redPhaseAmp*np.exp(-1/2*((phases-redPhaseOffset)/redPhaseSigma)**2)
        wavelengths = self.cal.convertToWvl(phases,row,col)
        if len(wavelengths)==0: #no wavecal for this pixel, so lambdaBinEdges came back empty
            wavelengths=phases
        pop.axes.plot(wavelengths,blueGaussFit,'b')
        pop.axes.plot(wavelengths,redGaussFit,'r')
        pop.draw()


    def showPixelStdVsIntTime(self,row,col):
        intTimes = [1,2,3,5,10,15,30]
        spectrumVsIntTimeVsTime = []
        for intTime in intTimes:
            spectrumInTime = []
            for iOb,ob in enumerate(self.skyObList):
                for sec in range(0,ob.getFromHeader('exptime'),intTime):
                    x = ob.getPixelSpectrum(pixelRow=row,pixelCol=col,firstSec=sec,integrationTime=intTime,weighted=True)
                    spectrum = x['spectrum']
                    spectrum = np.convolve(spectrum,np.ones(self.rebinSpecBins),'same')[self.firstAfterConvolve::self.rebinSpecBins]
                    
                    spectrumInTime.append(spectrum)
            spectrumInTime = np.array(spectrumInTime)
            spectrumVsIntTimeVsTime.append(spectrumInTime)
        #resulting array indexed as
        #spectrumVsIntTimeVsTime[iIntTime][iTimeChunk][iWvlBin]

        #sum over wavelength for total counts
        countsVsIntTimeVsTime = [np.sum(spectrumInTime,axis=1) for spectrumInTime in spectrumVsIntTimeVsTime]
        #countsVsIntTimeVsTime[iIntTime][iTimeChunk]

        countStds = [np.std(countsVsTime) for countsVsTime in countsVsIntTimeVsTime]
        countStds = np.array(countStds)
        countSqrts = [np.sqrt(np.median(countsVsTime)) for countsVsTime in countsVsIntTimeVsTime]
        countSqrts = np.array(countSqrts)
        spectrumStds = [np.std(spectrumVsTime,axis=0) for spectrumVsTime in spectrumVsIntTimeVsTime]
        spectrumSqrts = [np.sqrt(np.median(spectrumVsTime,axis=0)) for spectrumVsTime in spectrumVsIntTimeVsTime]
        spectrumStds = np.array(spectrumStds)
        spectrumSqrts = np.array(spectrumSqrts)

            
        pop = PopUp(parent=self,title='showPixelStdVsIntTime')
        pop.axes.set_xlabel('integration time (s)')
        pop.axes.set_ylabel('normalized $\sigma$')
        pop.axes.plot(intTimes,countSqrts/np.max(countSqrts),'k--',
            label=r'$\sqrt{N}$')
        pop.axes.plot(intTimes,countStds/np.max(countSqrts),'k',
            label=r'%d-%d $\AA$'%(self.rebinnedWvlEdges[0],self.rebinnedWvlEdges[-1]))
        nBins = np.shape(spectrumStds)[1]
        for iBin in xrange(nBins):
            pop.axes.plot(intTimes,
                spectrumStds[:,iBin]/np.max(spectrumSqrts[:,iBin]),
                c=cm.jet((iBin+1.0)/nBins),
                label=r'%d-%d $\AA$'%(self.rebinnedWvlEdges[iBin],
                    self.rebinnedWvlEdges[iBin+1]))
        pop.axes.legend(loc='upper left')
        pop.axes.set_title('Sky Pixel StdDev vs Int Time (%d,%d)'%(row,col))
        pop.draw()

    def showPixelRawPhaseHist(self,row,col):
        phases = np.array([],dtype=np.double)
        for iOb,ob in enumerate(self.obList):
            x = ob.getTimedPacketList(row,col)
            phases=np.append(phases,(np.array(x['peakHeights'],dtype=np.double)-np.array(x['baselines'],dtype=np.double)))
        pop = PopUp(parent=self,title='showPixelRawPhaseHist')
        nBins=np.max(phases)-np.min(phases)
        histPhases,binEdges = np.histogram(phases,bins=nBins)
        pop.axes.step(binEdges[:-1],histPhases,where='post')
        pop.axes.set_xlabel('peak-baseline')
        pop.axes.set_title('Peaks-Baselines')
        pop.draw()

    def showPixelRawBaselineHist(self,row,col):
        baselines = np.array([],dtype=np.double)
        for iOb,ob in enumerate(self.obList):
            x = ob.getTimedPacketList(row,col)
            baselines=np.append(baselines,np.array(x['baselines'],dtype=np.double))
        pop = PopUp(parent=self,title='showPixelRawBaselineHist')
        nBins=np.max(baselines)-np.min(baselines)
        histBaselines,binEdges = np.histogram(baselines,bins=nBins)
        pop.axes.step(binEdges[:-1],histBaselines,where='post')
        pop.axes.set_xlabel('baseline')
        pop.axes.set_title('Baselines')
        pop.draw()

    def showPixelRawPeakHist(self,row,col):
        peaks = np.array([],dtype=np.double)
        for iOb,ob in enumerate(self.obList):
            x = ob.getTimedPacketList(row,col)
            peaks = np.append(peaks,np.array(x['peakHeights'],dtype=np.double))
        pop = PopUp(parent=self,title='showPixelRawPeakHist')
        nBins=np.max(peaks)-np.min(peaks)
        histPeaks,binEdges = np.histogram(peaks,bins=nBins)
        pop.axes.step(binEdges[:-1],histPeaks,where='post')
        pop.axes.set_xlabel('peak')
        pop.axes.set_title('Packet Peaks (No Baseline Subtracted)')
        pop.draw()

    def showPixelFlatWeights(self,row,col):
        pop = PopUp(parent=self,title='showPixelFlatWeights')
        for iFlat,flatInfo in enumerate(self.flatInfos):
            weights = flatInfo['weights'][row,col]
            flatSpectra = flatInfo['spectra'][row,col]
            flatMedians = flatInfo['median']
            deltaFlatSpectra = np.sqrt(flatSpectra)
            deltaWeights = weights*deltaFlatSpectra/flatSpectra
            color=cm.jet((iFlat+1.)/len(self.flatInfos))
            wvlBinCenters = self.wvlBinEdges[:-1]+np.diff(self.wvlBinEdges)/2.
            pop.axes.plot(self.wvlBinEdges[:-1],weights,linestyle='-',label=self.params['flatInfoFiles'][iFlat],color=color,)
            pop.axes.errorbar(wvlBinCenters,weights,linestyle=',',yerr=deltaWeights,color=color)
        pop.axes.set_xlabel(r'$\lambda$ ($\AA$)')
        pop.axes.set_ylabel(r'Weights')
        pop.axes.set_title('Flat Weights')
        pop.axes.legend(loc='lower right')
        pop.draw()

    def showTwilightPixelSpectrum(self,row,col):
        spectrum = self.twilightSpectra[row,col]
        pop = PopUp(parent=self,title='showTwilightPixelSpectrum')
        pop.axes.step(self.wvlBinEdges[:-1],spectrum,where='post')
        pop.axes.set_xlabel(r'$\lambda$ ($\AA$)')
        pop.axes.set_ylabel(r'total counts')
        pop.axes.set_title('twilight spectrum (%d,%d) '%(row,col))
        pop.draw()

    def showTwilightPixelStdVsFlux(self,row,col):
        spectrumVsFluxVsTime = []
        for iOb,ob in enumerate(self.twilightObList):
            spectrumInTime = []
            for sec in range(0,ob.getFromHeader('exptime'),self.twilightIntTime):
                x = ob.getPixelSpectrum(pixelRow=row,pixelCol=col,firstSec=sec,integrationTime=self.twilightIntTime,weighted=True)
                spectrum = x['spectrum']
                binEdges = x['wvlBinEdges']
                spectrum = np.convolve(spectrum,np.ones(self.rebinSpecBins),'same')[self.firstAfterConvolve::self.rebinSpecBins]
                spectrumInTime.append(spectrum)
            spectrumInTime = np.array(spectrumInTime)
            spectrumVsFluxVsTime.append(spectrumInTime)
        spectrumVsFluxVsTime = np.array(spectrumVsFluxVsTime)


        #resulting array indexed as
        #spectrumVsFluxVsTime[iOb][iTimeChunk][iWvlBin]

        #sum over wavelength for total counts
        countsVsFluxVsTime = [np.sum(spectrumInTime,axis=1) for spectrumInTime in spectrumVsFluxVsTime]
        #countsVsFluxVsTime[iFlux][iTimeChunk]

        countStds = [np.std(countsVsTime) for countsVsTime in countsVsFluxVsTime]
        fluxes = [np.median(countsVsTime) for countsVsTime in countsVsFluxVsTime]
        fluxes = np.array(fluxes)
        countStds = np.array(countStds)
        countSqrts = [np.sqrt(np.median(countsVsTime)) for countsVsTime in countsVsFluxVsTime]
        countSqrts = np.array(countSqrts)
        spectrumStds = [np.std(spectrumVsTime,axis=0) for spectrumVsTime in spectrumVsFluxVsTime]
        spectrumSqrts = [np.sqrt(np.median(spectrumVsTime,axis=0)) for spectrumVsTime in spectrumVsFluxVsTime]
        spectrumStds = np.array(spectrumStds)
        spectrumSqrts = np.array(spectrumSqrts)

        pop = PopUp(parent=self,title='showTwilightPixelStdVsFlux')
        pop.axes.set_xlabel('median counts')
        pop.axes.set_ylabel('normalized $\sigma$')
        pop.axes.plot(fluxes,countSqrts/np.max(countSqrts),'k--',
            label=r'$\sqrt{N}$')
        pop.axes.plot(fluxes,countStds/np.max(countSqrts),'k',
            label=r'%d-%d $\AA$'%(self.rebinnedWvlEdges[0],self.rebinnedWvlEdges[-1]))
        nBins = np.shape(spectrumStds)[1]
        for iBin in xrange(nBins):
            pop.axes.plot(fluxes,
                spectrumStds[:,iBin]/np.max(spectrumSqrts[:,iBin]),
                c=cm.jet((iBin+1.0)/nBins),
                label=r'%d-%d $\AA$'%(self.rebinnedWvlEdges[iBin],
                    self.rebinnedWvlEdges[iBin+1]))
        pop.axes.legend(loc='upper left')
        pop.axes.set_title('Normalized Standard Deviation vs Twilight Flux, (%d,%d)'%(row,col))
        pop.draw()

    def getChisq(self,row,col):
        spectrum = self.twilightSpectra[row,col]

        weights = self.flatInfo['weights'][row,col]
        flatSpectra = self.flatInfo['spectra'][row,col]
        flatMedians = self.flatInfo['median']
        deltaFlatSpectra = np.sqrt(flatSpectra)
        deltaWeights = self.flatInfo['deltaWeights'][row,col]#weights*deltaFlatSpectra/flatSpectra
        poissonDeltaSpectra = np.sqrt(spectrum)

        rawSpectrum = spectrum/weights
        deltaRawSpectrum = np.sqrt(rawSpectrum)
        deltaSpectra = spectrum*np.sqrt((deltaWeights/weights)**2+(deltaRawSpectrum/rawSpectrum)**2)
        diffSpectrum = (spectrum-self.averageTwilightSpectrum)
        percentDiffSpectrum = 100.* diffSpectrum/self.averageTwilightSpectrum
        #deltaDiffSpectrum = np.sqrt(deltaSpectra**2+deltaSpectra2**2)
        deltaDiffSpectrum = np.array(deltaSpectra)
        #deltaDiffSpectrum[np.isnan(deltaDiffSpectrum)] = 0
        deltaPercentDiffSpectrum = 100.*deltaDiffSpectrum/self.averageTwilightSpectrum
        nDeltaFromZero = diffSpectrum/deltaDiffSpectrum
        chisqSumTerms =diffSpectrum**2/deltaDiffSpectrum**2
        chisqSumTerms = chisqSumTerms[~np.isnan(chisqSumTerms)]
        chisq = np.sum(chisqSumTerms)
        degreesOfFreedom=sum(~np.isnan(chisqSumTerms))-1
        reducedChisq = chisq/degreesOfFreedom
        return {'reducedChisq':reducedChisq,
                'percentDiffSpectrum':percentDiffSpectrum,
                'deltaPercentDiffSpectrum':deltaPercentDiffSpectrum,
                'nDeltaFromZero':nDeltaFromZero,
                'degreesOfFreedom':degreesOfFreedom,
                'chisq':chisq}
                

    def showTwilightPixelDeviationFromMedian(self,row,col):
        x = self.getChisq(row,col)
        reducedChisq = x['reducedChisq']
        chisq = x['chisq']
        percentDiffSpectrum = x['percentDiffSpectrum']
        deltaPercentDiffSpectrum = x['deltaPercentDiffSpectrum']
        nDeltaFromZero = x['nDeltaFromZero']
        degreesOfFreedom = x['degreesOfFreedom']
        print 'reduced chisq =',reducedChisq
        print 'P-value =',1-chi2.cdf(chisq,degreesOfFreedom)

        pop = PopUp(parent=self,title='showTwilightPixelDeviationFromMedian')
        pop.axes.errorbar(self.wvlBinEdges[:-1],percentDiffSpectrum,linestyle='-',color='k',yerr=deltaPercentDiffSpectrum)
        pop.axes.set_xlabel(r'$\lambda$ ($\AA$)')
        pop.axes.set_ylabel(r'percent difference')
        pop.axes.plot(self.wvlBinEdges[:-1],len(self.wvlBinEdges[:-1])*[0],'gray')
        axes2 = pop.axes.twinx()
        axes2.plot(self.wvlBinEdges[:-1],nDeltaFromZero,'m',alpha=.7)
        align_yaxis(pop.axes,0,axes2,0)
        axes2.set_ylabel(r'(pixelSpectrum-avgSpectrum)/$\sigma$',color='m')
        pop.axes.set_title('Deviation from Avg Spectrum (%d,%d)'%(row,col))
        pop.draw()

        weights = self.flatInfo['weights'][row,col]
        pop = PopUp(parent=self,title='showTwilightPixelDeviationFromMedian')
        pop.axes.step(self.wvlBinEdges[:-1],self.averageTwilightSpectrum/self.wvlBinWidths,'k',label='avg')
        pop.axes.step(self.wvlBinEdges[:-1],self.twilightSpectra[row,col]/self.wvlBinWidths,'b',label='weighted')
        pop.axes.step(self.wvlBinEdges[:-1],(self.twilightSpectra[row,col]/weights)/self.wvlBinWidths,'r',label='raw')
        pop.axes.set_xlabel(r'$\lambda$ ($\AA$)')
        pop.axes.set_ylabel(r'counts per $\AA$')
        pop.axes.set_title('Twilight Spectrum (%d,%d)'%(row,col))
        pop.axes.legend(loc='lower right')
        pop.draw()


def align_yaxis(ax1, v1, ax2, v2):
    """
    adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1
    Taken from http://stackoverflow.com/questions/10481990/matplotlib-axis-with-two-scales-shared-origin
    """
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)
            
        
def main():
    np.seterr(divide='ignore',invalid='ignore')
    app = QApplication(sys.argv)
    form = AppForm()
    form.loadStacks()
    form.prepareForClickPlots()
    form.arrayPlots()
    form.plotWeightedImage()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
