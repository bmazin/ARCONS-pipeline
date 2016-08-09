#!/bin/python
'''
fluxCal.py

Created by Seth Meeker on 11-21-2012
Modified on 02-16-2015 to perform absolute fluxCal with point sources

Opens ARCONS observation of a spectrophotometric standard star and 
associated wavelength cal file, reads in all photons and converts to energies. 
Bins photons to generate a spectrum, then divides this into the known spectrum 
of the object to create a Sensitivity curve.  This curve is then written out to 
h5 file.

Flags are associated with each pixel - see headers/pipelineFlags
for descriptions. Note some flags are set here, others are set
later on when creating photon lists.
'''

import sys,os
import tables
import numpy as np
from scipy import interpolate
from scipy.optimize.minpack import curve_fit
import matplotlib.pyplot as plt
from photometry import LightCurve
from util.FileName import FileName
from util.ObsFile import ObsFile
from util import MKIDStd
from util.readDict import readDict
from util.utils import rebin
from util.utils import gaussianConvolution
from util.utils import makeMovie
from util.utils import fitBlackbody
import hotpix.hotPixels as hp
from scipy.optimize.minpack import curve_fit
from scipy import interpolate
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from headers import pipelineFlags
import figureHeader


class FluxCal:
    def __init__(self,paramFile,plots=False,verbose=False):
        """
        Opens flux file, prepares standard spectrum, and calculates flux factors for the file.
        Method is provided in param file. If 'relative' is selected, an obs file with standard star defocused over
        the entire array is expected, with accompanying sky file to do sky subtraction.
        If any other method is provided, 'absolute' will be done by default, wherein a point source is assumed
        to be present. The obs file is then broken into spectral frames with photometry (psf or aper) performed 
        on each frame to generate the ARCONS observed spectrum.
        """
        self.verbose=verbose
        self.plots = plots

        self.params = readDict()
        self.params.read_from_file(paramFile)
        
        run = self.params['run']
        sunsetDate = self.params['fluxSunsetLocalDate']
        self.fluxTstamp = self.params['fluxTimestamp']
        skyTstamp = self.params['skyTimestamp']
        wvlSunsetDate = self.params['wvlCalSunsetLocalDate']
        wvlTimestamp = self.params['wvlCalTimestamp']
        flatCalFileName = self.params['flatCalFileName']
        needTimeAdjust = self.params['needTimeAdjust']
        self.deadtime = float(self.params['deadtime']) #from firmware pulse detection
        self.timeSpacingCut = self.params['timeSpacingCut']
        bLoadBeammap = self.params.get('bLoadBeammap',False)
        self.method = self.params['method']
        self.objectName = self.params['object']
        self.r = float(self.params['energyResolution'])
        self.photometry = self.params['photometry']
        self.centroidRow = self.params['centroidRow']
        self.centroidCol = self.params['centroidCol']
        self.aperture = self.params['apertureRad']
        self.annulusInner = self.params['annulusInner']
        self.annulusOuter = self.params['annulusOuter']
        self.collectingArea = self.params['collectingArea']
        self.startTime = self.params['startTime']
        self.intTime = self.params['integrationTime']

        fluxFN = FileName(run=run,date=sunsetDate,tstamp=self.fluxTstamp)
        self.fluxFileName = fluxFN.obs()
        self.fluxFile = ObsFile(self.fluxFileName)

        if self.plots:
            self.plotSavePath = os.environ['MKID_PROC_PATH']+os.sep+'fluxCalSolnFiles'+os.sep+run+os.sep+sunsetDate+os.sep+'plots'+os.sep
            if not os.path.exists(self.plotSavePath): os.mkdir(self.plotSavePath)
            if self.verbose: print "Created directory %s"%self.plotSavePath

        obsFNs = [fluxFN]
        self.obsList = [self.fluxFile]

        if self.startTime in ['',None]: self.startTime=0
        if self.intTime in ['',None]: self.intTime=-1

        if self.method=="relative":
            try:
                print "performing Relative Flux Calibration"
                skyFN = FileName(run=run,date=sunsetDate,tstamp=skyTstamp)
                self.skyFileName = skyFN.obs()
                self.skyFile = ObsFile(self.skyFileName)
                obsFNs.append(skyFN)
                self.obsList.append(self.skyFile)
            except:
                print "For relative flux calibration a sky file must be provided in param file"
                self.__del__()
        else:
            self.method='absolute'
            print "performing Absolute Flux Calibration"

        if self.photometry not in ['aperture','PSF']: self.photometry='PSF' #default to PSF fitting if no valid photometry selected

        timeMaskFileNames = [fn.timeMask() for fn in obsFNs]
        timeAdjustFileName = FileName(run=run).timeAdjustments()

        #make filename for output fluxCalSoln file
        self.fluxCalFileName = FileName(run=run,date=sunsetDate,tstamp=self.fluxTstamp).fluxSoln()
        print "Creating flux cal: %s"%self.fluxCalFileName

        if wvlSunsetDate != '':
            wvlCalFileName = FileName(run=run,date=wvlSunsetDate,tstamp=wvlTimestamp).calSoln()
        if flatCalFileName =='':
            flatCalFileName=FileName(obsFile=self.fluxFile).flatSoln()

        #load cal files for flux file and, if necessary, sky file
        for iObs,obs in enumerate(self.obsList):
            if bLoadBeammap:
                print 'loading beammap',os.environ['MKID_BEAMMAP_PATH']
                obs.loadBeammapFile(os.environ['MKID_BEAMMAP_PATH'])
            if wvlSunsetDate != '':
                obs.loadWvlCalFile(wvlCalFileName)
            else:
                obs.loadBestWvlCalFile()

            obs.loadFlatCalFile(flatCalFileName)
            obs.setWvlCutoffs(-1,-1)

            if needTimeAdjust:
                obs.loadTimeAdjustmentFile(timeAdjustFileName)
            timeMaskFileName = timeMaskFileNames[iObs]
            print timeMaskFileName

            if not os.path.exists(timeMaskFileName):
                print 'Running hotpix for ',obs
                hp.findHotPixels(obsFile=obs,outputFileName=timeMaskFileName,fwhm=np.inf,useLocalStdDev=True)
                print "Flux cal/sky file pixel mask saved to %s"%(timeMaskFileName)
            obs.loadHotPixCalFile(timeMaskFileName)
            if self.verbose: print "Loaded hot pixel file %s"%timeMaskFileName

        #get flat cal binning information since flux cal will need to match it
        self.wvlBinEdges = self.fluxFile.flatCalFile.root.flatcal.wavelengthBins.read()
        self.nWvlBins = self.fluxFile.flatWeights.shape[2]
        self.binWidths = np.empty((self.nWvlBins),dtype=float)
        self.binCenters = np.empty((self.nWvlBins),dtype=float)
        for i in xrange(self.nWvlBins):
            self.binWidths[i] = self.wvlBinEdges[i+1]-self.wvlBinEdges[i]
            self.binCenters[i] = (self.wvlBinEdges[i]+(self.binWidths[i]/2.0))

        if self.method=='relative':
            print "Extracting ARCONS flux and sky spectra"
            self.loadRelativeSpectrum()
            print "Flux Spectrum loaded"
            self.loadSkySpectrum()
            print "Sky Spectrum loaded"
        elif self.method=='absolute':
            print "Extracting ARCONS point source spectrum"
            self.loadAbsoluteSpectrum()

        print "Loading standard spectrum"
        try:
            self.loadStdSpectrum(self.objectName)
        except KeyError:
            print "Invalid spectrum object name"
            self.__del__()
            sys.exit()

        print "Generating sensitivity curve"
        self.calculateFactors()
        print "Sensitivity Curve calculated"
        print "Writing fluxCal to file %s"%self.fluxCalFileName
        self.writeFactors(self.fluxCalFileName)
        
        if self.plots: self.makePlots()

        print "Done"

    def __del__(self):
        try:
            self.fluxFile.close()
            self.calFile.close()
        except AttributeError:#fluxFile was never defined
            pass

    def getDeadTimeCorrection(self, obs): #WRONG RIGHT NOW. NEEDS TO HAVE RAW COUNTS SUMMED, NOT CUBE WHICH EXCLUDES NOISE TAIL
        if self.verbose: print "Making raw cube to get dead time correction"
        cubeDict = obs.getSpectralCube(firstSec=self.startTime, integrationTime=self.intTime, weighted=False, fluxWeighted=False)

        cube= np.array(cubeDict['cube'], dtype=np.double)
        wvlBinEdges= cubeDict['wvlBinEdges']
        effIntTime= cubeDict['effIntTime']
        if self.verbose: print "median effective integration time = ", np.median(effIntTime)

        nWvlBins=len(wvlBinEdges)-1
        if self.verbose: print "cube shape ", np.shape(cube)
        if self.verbose: print "effIntTime shape ", np.shape(effIntTime)

        #add third dimension to effIntTime for  broadcasting
        effIntTime = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
        #put cube into counts/s in each pixel
        cube /= effIntTime

        #CALCULATE DEADTIME CORRECTION
        #NEED TOTAL COUNTS PER SECOND FOR EACH PIXEL TO DO PROPERLY
        #ASSUMES SAME CORRECTION FACTOR APPLIED FOR EACH WAVELENGTH, MEANING NO WL DEPENDANCE ON DEAD TIME EFFECT
        DTCorr = np.zeros((np.shape(cube)[0],np.shape(cube)[1]),dtype=float)
        for f in range(0,np.shape(cube)[2]):
            #if self.verbose: print cube[:,:,f]
            #if self.verbose: print '-----------------------'
            DTCorr += cube[:,:,f]
            #if self.verbose: print DTCorr
            #if self.verbose: print '\n=====================\n'
        #Correct for firmware dead time (100us in 2012 ARCONS firmware)
        DTCorrNew=DTCorr/(1-DTCorr*self.deadtime)
        CorrFactors = DTCorrNew/DTCorr #This is what the frames need to be multiplied by to get their true values
        if self.verbose: print "Dead time correction factors: ", CorrFactors
        #add third dimension to CorrFactors for broadcasting
        CorrFactors = np.reshape(CorrFactors,np.shape(CorrFactors)+(1,))
        return CorrFactors

    def loadAbsoluteSpectrum(self):
        '''
        extract the ARCONS measured spectrum of the spectrophotometric standard by breaking data into spectral cube
        and performing photometry (aper or psf) on each spectral frame
        '''
        if self.verbose:print "Making spectral cube"
        cubeDict = self.fluxFile.getSpectralCube(firstSec=self.startTime, integrationTime=self.intTime, weighted=True, fluxWeighted=False)
        cube= np.array(cubeDict['cube'], dtype=np.double)
        effIntTime= cubeDict['effIntTime']
        if self.verbose: print "median effective integration time in flux file cube = ", np.median(effIntTime)
        if self.verbose: print "cube shape ", np.shape(cube)
        if self.verbose: print "effIntTime shape ", np.shape(effIntTime)

        #add third dimension to effIntTime for broadcasting
        effIntTime = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
        #put cube into counts/s in each pixel
        cube /= effIntTime
        #get dead time correction factors
        DTCorr = self.getDeadTimeCorrection(self.fluxFile)
        cube*=DTCorr #cube now in units of counts/s and corrected for dead time

        
        if self.plots and not 'figureHeader' in sys.modules:
            if self.verbose: print "Saving spectral frames as movie..."
            movieCube = np.zeros((self.nWvlBins,np.shape(cube)[0],np.shape(cube)[1]),dtype=float)
            for i in xrange(self.nWvlBins):
                movieCube[i,:,:] = cube[:,:,i]
            makeMovie(movieCube,frameTitles=self.binCenters,cbar=True,outName=self.plotSavePath+'FluxCal_Cube_%s.gif'%(self.objectName), normMin=0, normMax=50)
            if self.verbose: print "Movie saved in %s"%self.plotSavePath
        

        LCplot=False #light curve pop-ups not compatible with FLuxCal plotting 2/18/15
        #if self.photometry=='PSF': LCplot = False
        LC = LightCurve.LightCurve(verbose=self.verbose, showPlot=LCplot)

        self.fluxSpectrum=np.empty((self.nWvlBins),dtype=float)
        self.skySpectrum=np.zeros((self.nWvlBins),dtype=float)

        for i in xrange(self.nWvlBins):
            frame = cube[:,:,i]
            if self.verbose: print "%s photometry on frame %i of cube, central wvl = %f Angstroms"%(self.photometry,i,self.binCenters[i])
            if self.photometry == 'aperture':
                fDict = LC.performPhotometry(self.photometry,frame,[[self.centroidCol,self.centroidRow]],expTime=None,aper_radius = self.aperture, annulus_inner = self.annulusInner, annulus_outer = self.annulusOuter, interpolation="linear")
                self.fluxSpectrum[i] = fDict['flux']
                self.skySpectrum[i] = fDict['skyFlux']
                print "Sky estimate = ", fDict['skyFlux']
            else:
                fDict = LC.performPhotometry(self.photometry,frame,[[self.centroidCol,self.centroidRow]],expTime=None,aper_radius = self.aperture)
                self.fluxSpectrum[i] = fDict['flux']
        
        self.fluxSpectrum=self.fluxSpectrum/self.binWidths/self.collectingArea #spectrum now in counts/s/Angs/cm^2
        self.skySpectrum=self.skySpectrum/self.binWidths/self.collectingArea

        return self.fluxSpectrum, self.skySpectrum

    def loadRelativeSpectrum(self):
        self.fluxSpectra = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        self.fluxEffTime = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                count = self.fluxFile.getPixelCount(iRow,iCol)
                fluxDict = self.fluxFile.getPixelSpectrum(iRow,iCol,weighted=True,firstSec=0,integrationTime=-1)
                self.fluxSpectra[iRow][iCol],self.fluxEffTime[iRow][iCol] = fluxDict['spectrum'],fluxDict['effIntTime']
        self.fluxSpectra = np.array(self.fluxSpectra)
        self.fluxEffTime = np.array(self.fluxEffTime)
        DTCorr = self.getDeadTimeCorrection(self.fluxFile)
        #print "Bin widths = ",self.binWidths
        self.fluxSpectra = self.fluxSpectra/self.binWidths/self.fluxEffTime*DTCorr
        self.fluxSpectrum = self.calculateMedian(self.fluxSpectra) #find median of subtracted spectra across whole array
        return self.fluxSpectrum

    def loadSkySpectrum(self):
        self.skySpectra = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        self.skyEffTime = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)] 
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                count = self.skyFile.getPixelCount(iRow,iCol)
                skyDict = self.skyFile.getPixelSpectrum(iRow,iCol,weighted=True,firstSec=0,integrationTime=-1)
                self.skySpectra[iRow][iCol],self.skyEffTime[iRow][iCol] = skyDict['spectrum'],skyDict['effIntTime']
        self.skySpectra = np.array(self.skySpectra)
        self.skyEffTime = np.array(self.skyEffTime)
        DTCorr = self.getDeadTimeCorrection(self.skyFile)
        self.skySpectra = self.skySpectra/self.binWidths/self.skyEffTime*DTCorr
        self.skySpectrum = self.calculateMedian(self.skySpectra) #find median of subtracted spectra across whole array
        return self.skySpectrum

    def loadStdSpectrum(self, objectName="G158-100"):
        #import the known spectrum of the calibrator and rebin to the histogram parameters given
        #must be imported into array with dtype float so division later does not have error
        std = MKIDStd.MKIDStd()
        a = std.load(objectName)
        a = std.countsToErgs(a) #convert std spectrum to ergs/s/Angs/cm^2 for BB fitting and cleaning
        self.stdWvls = np.array(a[:,0])
        self.stdFlux = np.array(a[:,1]) #std object spectrum in ergs/s/Angs/cm^2

        if self.plots:
            #create figure for plotting standard spectrum modifications      
            self.stdFig = plt.figure()
            self.stdAx = self.stdFig.add_subplot(111)
            plt.xlim(3500,12000)
            plt.plot(self.stdWvls,self.stdFlux*1E15,linewidth=1,color='grey',alpha=0.75)
            
        convX_rev,convY_rev = self.cleanSpectrum(self.stdWvls,self.stdFlux)
        convX = convX_rev[::-1] #convolved spectrum comes back sorted backwards, from long wvls to low which screws up rebinning
        convY = convY_rev[::-1]
        #rebin cleaned spectrum to flat cal's wvlBinEdges
        newa = rebin(convX,convY,self.wvlBinEdges)
        rebinnedWvl = np.array(newa[:,0])
        rebinnedFlux = np.array(newa[:,1])
        if self.plots:
            #plot final resampled spectrum
            plt.plot(convX,convY*1E15,color='blue')
            plt.step(rebinnedWvl,rebinnedFlux*1E15,color = 'black',where='mid')
            plt.legend(['%s Spectrum'%self.objectName,'Blackbody Fit','Gaussian Convolved Spectrum','Rebinned Spectrum'],'upper right', numpoints=1)
            plt.xlabel(ur"Wavelength (\r{A})")
            plt.ylabel(ur"Flux (10$^{-15}$ ergs s$^{-1}$ cm$^{-2}$ \r{A}$^{-1}$)")
            plt.ylim(0.9*min(rebinnedFlux)*1E15, 1.1*max(rebinnedFlux)*1E15)
            plt.savefig(self.plotSavePath+'FluxCal_StdSpectrum_%s.eps'%self.objectName,format='eps')
        
        #convert standard spectrum back into counts/s/angstrom/cm^2
        newa = std.ergsToCounts(newa)
        self.binnedSpectrum = np.array(newa[:,1])


    def cleanSpectrum(self,x,y):
        ##=============== BB Fit to extend spectrum beyond 11000 Angstroms ==================
        fraction = 1.0/3.0
        nirX = np.arange(int(x[(1.0-fraction)*len(x)]),20000)
        T, nirY = fitBlackbody(x,y,fraction=fraction,newWvls=nirX,tempGuess=5600)
        
        if self.plots: plt.plot(nirX,nirY*1E15,linestyle='--',linewidth=2, color="black",alpha=0.5)

        extendedWvl = np.concatenate((x,nirX[nirX>max(x)]))
        extendedFlux = np.concatenate((y,nirY[nirX>max(x)]))
        ##======= Gaussian convolution to smooth std spectrum to MKIDs median resolution ========
        newX, newY = gaussianConvolution(extendedWvl,extendedFlux,xEnMin=0.005,xEnMax=6.0,xdE=0.001,fluxUnits = "lambda",r=self.r,plots=False)
        return newX, newY


    def calculateFactors(self):
        """
        Calculate the sensitivity spectrum: the weighting factors that correct the flat calibrated spectra to the real spectra
        
        For relative calibration:
        First subtract sky spectrum from ARCONS observed spectrum. Then take median of this spectrum as it should be identical 
        across the array, assuming the flat cal has done its job. Then divide this into the known spectrum of the object.
        
        For absolute calibration:
        self.fluxSpectra already has sky subtraction included. Simply divide this spectrum into the known standard spectrum.
        """
        self.subtractedSpectrum = self.fluxSpectrum - self.skySpectrum
        self.subtractedSpectrum = np.array(self.subtractedSpectrum,dtype=float) #cast as floats so division does not fail later

        if self.method=='relative':      
            normWvl = 5500 #Angstroms. Choose an arbitrary wvl to normalize the relative correction at
            ind = np.where(self.wvlBinEdges >= normWvl)[0][0]-1
            self.subtractedSpectrum = self.subtractedSpectrum/(self.subtractedSpectrum[ind]) #normalize
            self.binnedSpectrum = self.binnedSpectrum/(self.binnedSpectrum[ind]) #normalize treated Std spectrum while we are at it

        #Calculate FluxCal factors
        self.fluxFactors = self.binnedSpectrum/self.subtractedSpectrum

        #self.fluxFlags = np.zeros(np.shape(self.fluxFactors),dtype='int')
        self.fluxFlags = np.empty(np.shape(self.fluxFactors),dtype='int')
        self.fluxFlags.fill(pipelineFlags.fluxCal['good'])   #Initialise flag array filled with 'good' flags. JvE 5/1/2013.
        #set factors that will cause trouble to 1
        #self.fluxFlags[self.fluxFactors == np.inf] = 1
        self.fluxFlags[self.fluxFactors == np.inf] = pipelineFlags.fluxCal['infWeight']   #Modified to use flag dictionary - JvE 5/1/2013
        self.fluxFactors[self.fluxFactors == np.inf]=1.0
        self.fluxFlags[np.isnan(self.fluxFactors)] = pipelineFlags.fluxCal['nanWeight']   #Modified to use flag dictionary - JvE 5/1/2013
        self.fluxFactors[np.isnan(self.fluxFactors)]=1.0        
        self.fluxFlags[self.fluxFactors <= 0]=pipelineFlags.fluxCal['LEzeroWeight']   #Modified to use flag dictionary - JvE 5/1/2013
        self.fluxFactors[self.fluxFactors <= 0]=1.0


    def calculateMedian(self, spectra):
        spectra2d = np.reshape(spectra,[self.nRow*self.nCol,self.nWvlBins])
        wvlMedian = np.empty(self.nWvlBins,dtype=float)
        for iWvl in xrange(self.nWvlBins):
            spectrum = spectra2d[:,iWvl]
            goodSpectrum = spectrum[spectrum != 0]#dead pixels need to be taken out before calculating medians
            wvlMedian[iWvl] = np.median(goodSpectrum)
        return wvlMedian


    def makePlots(self):
        """
        Output all debugging plots of ARCONS sky and object spectra, known calibrator spectrum, and sensitivity curve
        """
        scratchDir = os.getenv('MKID_PROC_PATH')
        fluxDir = self.plotSavePath
        fluxCalBase = 'FluxCal_%s'%self.objectName       
        plotFileName = fluxCalBase+".pdf"
        fullFluxPlotFileName = os.path.join(fluxDir,plotFileName)
        
        #uncomment to make some plots for the paper. Proper formatting Will also require figureheader to be imported and for movie making to be turned off
        self.paperFig = plt.figure()
        self.paperAx = self.paperFig.add_subplot(111)
        plt.xlim(4000,11000)
        plt.plot(self.binCenters,self.fluxFactors,linewidth=3,color='black')
        plt.xlabel(ur"Wavelength (\r{A})")
        plt.ylabel(ur"Spectral Calibration Curve")
        plt.ylim(0,150)
        plt.savefig(self.plotSavePath+'FluxCal_Sensitivity_%s.eps'%self.objectName,format='eps')

        #save throughput as a .npz file that other code uses when making paper plots
        np.savez(self.plotSavePath+'%s_%s_throughput.npz'%(self.objectName.strip(),self.fluxTstamp),throughput=1.0/self.fluxFactors,wvls=self.binCenters)

        pp = PdfPages(fullFluxPlotFileName)
        #plt.rcParams['font.size'] = 2

        wvls = self.binCenters

        plt.figure()
        ax1 = plt.subplot(111)
        ax1.set_title('ARCONS median flat cal\'d flux in counts')
        plt.plot(wvls,self.fluxSpectrum)
        pp.savefig()

        plt.figure()
        ax2 = plt.subplot(111)
        ax2.set_title('ARCONS median flat cal\'d sky in counts')
        plt.plot(wvls,self.skySpectrum)
        pp.savefig()

        plt.figure()
        ax3 = plt.subplot(111)
        ax3.set_title('Flux data minus sky in counts')
        plt.plot(wvls,self.subtractedSpectrum)
        pp.savefig()

        plt.figure()
        ax4 = plt.subplot(111)
        ax4.set_title('Std Spectrum of %s'%(self.objectName))
        plt.plot(self.stdWvls,self.stdFlux)
        pp.savefig()

        plt.figure()
        ax5 = plt.subplot(111)
        ax5.set_title('Binned Std Spectrum')
        plt.plot(wvls,self.binnedSpectrum)
        pp.savefig()

        plt.figure()
        ax6 = plt.subplot(111)
        ax6.set_title('Median Sensitivity Spectrum')
        ax6.set_xlim((3500,12000))
        #ax6.set_ylim((0,5))
        plt.plot(wvls,self.fluxFactors)
        pp.savefig()

        plt.figure()
        ax7 = plt.subplot(111)
        ax7.set_title('1/Sensitivity (Throughput)')
        ax7.set_xlim((3500,12000))
        ax7.set_ylim((0,.04))
        plt.plot(wvls,1.0/self.fluxFactors)
        pp.savefig()

        plt.figure()
        ax8 = plt.subplot(111)
        ax8.set_title('Flux Cal\'d ARCONS Spectrum of Std')
        plt.plot(wvls,self.fluxFactors*self.subtractedSpectrum)
        pp.savefig()
        pp.close()
        print "Saved Flux Cal plots to %s"%(fullFluxPlotFileName)
    

    def writeFactors(self,fluxCalFileName):
        """
        Write flux cal weights to h5 file
        """

        if os.path.isabs(fluxCalFileName) == True:
            fullFluxCalFileName = fluxCalFileName
        else:
            scratchDir = os.getenv('MKID_PROC_PATH')
            fluxDir = os.path.join(scratchDir,'fluxCalSolnFiles')
            fullFluxCalFileName = os.path.join(fluxDir,fluxCalFileName)

        try:
            fluxCalFile = tables.openFile(fullFluxCalFileName,mode='w')
        except:
            print 'Error: Couldn\'t create flux cal file, ',fullFluxCalFileName
            return

        calgroup = fluxCalFile.createGroup(fluxCalFile.root,'fluxcal','Table of flux calibration weights by wavelength')
        caltable = tables.Array(calgroup,'weights',object=self.fluxFactors,title='Flux calibration Weights indexed by wavelengthBin')
        flagtable = tables.Array(calgroup,'flags',object=self.fluxFlags,title='Flux cal flags indexed by wavelengthBin. 0 is Good')
        bintable = tables.Array(calgroup,'wavelengthBins',object=self.wvlBinEdges,title='Wavelength bin edges corresponding to third dimension of weights array')
        fluxCalFile.flush()
        fluxCalFile.close()
        print "Finished Flux Cal, written to %s"%(fullFluxCalFileName)


    def cleanSpectrum_old(self,x,y,objectName):
        ''' 
        function to take high resolution spectrum of standard star, extend IR coverage with 
        an exponential tail, then rebin down to ARCONS resolution. This function has since been
        deprecated with the current cleanSpectrum which uses a BB fit to extend IR coverage,
        and does the rebinning using a gaussian convolution. This is left in for reference.
        '''

        #locations and widths of absorption features in Angstroms
        #features = [3890,3970,4099,4340,4860,6564,6883,7619]
        #widths = [50,50,50,50,50,50,50,50]
        #for i in xrange(len(features)):
        #    #check for absorption feature in std spectrum
        #    ind = np.where((x<(features[i]+15)) & (x>(features[i]-15)))[0]
        #    if len(ind)!=0:
        #        ind = ind[len(ind)/2]
        #    #if feature is found (flux is higher on both sides of the specified wavelength where the feature should be)
        #    if y[ind]<y[ind+1] and y[ind]<y[ind-1]:
        #        #cut out width[i] around feature[i]
        #        inds = np.where((x >= features[i]+widths[i]) | (x <= features[i]-widths[i]))
        #        x = x[inds]
        #        y = y[inds]

        #fit a tail to the end of the spectrum to interpolate out to desired wavelength in angstroms
        fraction = 3.0/4.0
        newx = np.arange(int(x[fraction*len(x)]),20000)

        slopeguess = (np.log(y[-1])-np.log(y[fraction*len(x)]))/(x[-1]-x[fraction*len(x)])
        print "Guess at exponential slope is %f"%(slopeguess)
        guess_a, guess_b, guess_c = float(y[fraction*len(x)]), x[fraction*len(x)], slopeguess
        guess = [guess_a, guess_b, guess_c]

        fitx = x[fraction*len(x):]
        fity = y[fraction*len(x):]

        exp_decay = lambda fx, A, x0, t: A * np.exp((fx-x0) * t)

        params, cov = curve_fit(exp_decay, fitx, fity, p0=guess, maxfev=2000)
        A, x0, t= params
        print "A = %s\nx0 = %s\nt = %s\n"%(A, x0, t)
        best_fit = lambda fx: A * np.exp((fx-x0)*t)

        calcx = np.array(newx,dtype=float)
        newy = best_fit(calcx)

        #func = interpolate.splrep(x[fration*len(x):],y[fraction*len(x):],s=smooth)
        #newx = np.arange(int(x[fraction*len(x)]),self.wvlBinEdges[-1])
        #newy = interpolate.splev(newx,func)

        wl = np.concatenate((x,newx[newx>max(x)]))
        flux = np.concatenate((y,newy[newx>max(x)]))

        #new method, rebin data to grid of wavelengths generated from a grid of evenly spaced energy bins
        #R=7.0 at 4500
        #R=E/dE -> dE = R/E
        dE = 0.3936 #eV
        start = 1000 #Angs
        stop = 20000 #Angs
        enBins = ObsFile.makeWvlBins(dE,start,stop)
        rebinned = rebin(wl,flux,enBins)
        re_wl = rebinned[:,0]
        re_flux = rebinned[:,1]
        #plt.plot(re_wl,re_flux,color='r')
        
        re_wl = re_wl[np.isnan(re_flux)==False]
        re_flux = re_flux[np.isnan(re_flux)==False]

        start1 = self.wvlBinEdges[0]
        stop1 = self.wvlBinEdges[-1]
        #regrid downsampled data 

        new_wl = np.arange(start1,stop1)

        #print re_wl
        #print re_flux
        #print new_wl

        #weight=1.0/(re_flux)**(2/1.00)
        print len(re_flux)
        weight = np.ones(len(re_flux))
        #decrease weights near peak
        ind = np.where(re_flux == max(re_flux))[0]
        weight[ind] = 0.3
        for p in [1,2,3]:
            if p==1:
                wt = 0.3
            elif p==2:
                wt = 0.6
            elif p==3:
                wt = 0.7
            try:
                weight[ind+p] = wt
            except IndexError:
                 pass
            try:
                 if ind-p >= 0:
                     weight[ind-p] = wt
            except IndexError:
                pass
        weight[-4:] = 1.0
        #weight = [0.7,1,0.3,0.3,0.5,0.7,1,1,1]
        #print len(weight)
        #weight = re_flux/min(re_flux)
        #weight = 1.0/weight
        #weight = weight/max(weight)
        #print weight
        f = interpolate.splrep(re_wl,re_flux,w=weight,k=3,s=max(re_flux)**1.71)
        new_flux = interpolate.splev(new_wl,f,der=0)
        return new_wl, new_flux


if __name__ == '__main__':
    try:
        paramFile = sys.argv[1]
    except:
        paramFile = '/home/srmeeker/ARCONS-pipeline/params/fluxCal.dict'
    fc = FluxCal(paramFile, plots=True, verbose=True)





































