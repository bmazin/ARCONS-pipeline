#!/bin/python
'''
fluxCal.py

Created by Seth Meeker on 11-21-2012

Opens ARCONS observation of a spread out spectrophotometric standard star and 
associated wavelength cal file, reads in all photons and converts to energies. 
Bins photons to generate a spectrum, then divides this into the known spectrum 
of the object to create a Sensitivity curve.  This curve is then written out to 
h5 file.

Flags are associated with each pixel - see headers/pipelineFlags
for descriptions. Note some flags are set here, others are set
later one when creating photon lists.
'''

import sys,os
import tables
import numpy as np
from scipy import interpolate
from scipy.optimize.minpack import curve_fit
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from util.rebin import rebin
import hotpix.hotPixels as hp
from scipy.optimize.minpack import curve_fit
from scipy import interpolate
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from headers import pipelineFlags


class FluxCal:
    def __init__(self, fluxFileName=None, skyFileName=None, wvlCalFileName=None,flatCalFileName=None, objectName=None, fluxCalFileName=None, plots=False):
        '''
        opens flux files, applies wl cal, makes spectrum of ARCONS data, applies flat cal, calculates flux factors and writes to file
        '''
        #self.wvlBinWidth = 100 #angstroms
        #self.wvlStart = 3000 #angstroms
        #self.wvlStop = 13000 #angstroms
        #self.nWvlBins = int((self.wvlStop - self.wvlStart)/self.wvlBinWidth)

        self.loadFluxFile(fluxFileName)#load file with spectrophotometric standard star data
        self.fluxFile.loadWvlCalFile(wvlCalFileName)#load the wavelength cal to be applied to the flux file
        self.fluxFile.loadFlatCalFile(flatCalFileName)#load FlatCalFile for flux file

        #get flat cal binning information since flux cal will need to match it
        self.wvlBinEdges = self.fluxFile.flatCalFile.root.flatcal.wavelengthBins.read()
        self.nWvlBins = self.fluxFile.flatWeights.shape[2]

        self.loadSkyFile(skyFileName)#load file with sky data that accompanies standard star
        self.skyFile.loadWvlCalFile(wvlCalFileName)#same wvl cal file is used on the sky obs file
        self.skyFile.loadFlatCalFile(flatCalFileName)#same flat cal file is used on the sky file

        print "Getting hot pixel masks for flux file and sky file"
        self.fluxHotPixFile = self.getTimeMaskFileName(fluxFileName)
        if not os.path.exists(self.fluxHotPixFile):
            hp.findHotPixels(fluxFileName,self.fluxHotPixFile)
            print "Flux file pixel mask saved to %s"%(self.fluxHotPixFile)
        self.skyHotPixFile = self.getTimeMaskFileName(skyFileName)
        if not os.path.exists(self.skyHotPixFile):
            hp.findHotPixels(skyFileName,self.skyHotPixFile)
            print "Sky file pixel mask saved to %s"%(self.skyHotPixFile)
        self.skyFile.loadHotPixCalFile(self.skyHotPixFile)
        self.fluxFile.loadHotPixCalFile(self.fluxHotPixFile)
        print "Hot pixel masks loaded"
        print "Loading flux and sky spectra"
        self.loadFluxSpectra()
        print "Flux Spectra loaded"
        self.loadSkySpectra()
        print "Sky Spectra loaded"
        print "Loading standard spectrum"
        try:
            self.loadSpectrum(objectName)
        except KeyError:
            print "Invalid spectrum object name"
            self.__del__()
            sys.exit()
        print "Generating sensitivity curve"
        self.calculateFactors()
        print "Sensitivity Curve calculated"
        print "Writing fluxCal to file"
        self.writeFactors(fluxCalFileName)
        
        if plots == True:
            print "Making plots"
            self.makePlots(fluxFileName,skyFileName, objectName)
        print "Done"

    def __del__(self):
        try:
            self.fluxFile.close()
            self.calFile.close()
        except AttributeError:#fluxFile was never defined
            pass

    def loadFluxFile(self,fluxFileName):
        #open the hdf5 file
        self.fluxFile = ObsFile(fluxFileName)
        self.nRow = self.fluxFile.nRow
        self.nCol = self.fluxFile.nCol
        self.fluxTime = self.fluxFile.getFromHeader("exptime")

    def loadSkyFile(self,skyFileName):
        #open the sky hdf5 file which will be subtracted from the observation file
        self.skyFile = ObsFile(skyFileName)
        skynRow = self.skyFile.nRow
        skynCol = self.skyFile.nCol
        self.skyTime = self.skyFile.getFromHeader("exptime")
        if (skynRow != self.nRow) or (skynCol != self.nCol):
            print "Flux Calibration failed: Sky file not same array dimensions as Flux file"
            sys.exit(1)

    def getTimeMaskFileName(self, obsFileName):
        scratchDir = os.getenv('MKID_PROC_PATH')
        hotPixDir = os.path.join(scratchDir,'timeMasks')
        fileName = obsFileName.split('/')[-1]
        fileNameBase = fileName.split('_')[-1]
        newName = 'timeMask_'+fileNameBase
        fullHotPixFileName = os.path.join(hotPixDir,newName)
        return fullHotPixFileName


    def loadFluxSpectra(self):
        self.fluxSpectra = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        self.fluxEffTime = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                #print iRow,iCol,
                count = self.fluxFile.getPixelCount(iRow,iCol)
                #print count
                #self.fluxSpectra[iRow][iCol],self.wvlBinEdges = self.fluxFile.getPixelSpectrum(iRow,iCol,weighted=True,firstSec=0,integrationTime=-1)
                fluxDict = self.fluxFile.getPixelSpectrum(iRow,iCol,weighted=True,firstSec=0,integrationTime=-1)
                self.fluxSpectra[iRow][iCol],self.wvlBinEdges,self.fluxEffTime[iRow][iCol] = fluxDict['spectrum'],fluxDict['wvlBinEdges'],fluxDict['effIntTime']
        self.fluxSpectra = np.array(self.fluxSpectra)
        self.fluxEffTime = np.array(self.fluxEffTime)
        self.nWvlBins = len(self.wvlBinEdges)-1
        self.binWidths = np.empty(self.nWvlBins)
        for i in xrange(self.nWvlBins):
            self.binWidths[i] = self.wvlBinEdges[i+1]-self.wvlBinEdges[i]
        #print "Bin widths = ",self.binWidths
        self.fluxSpectra = self.fluxSpectra/self.binWidths
        return self.fluxSpectra

    def loadSkySpectra(self):
        self.skySpectra = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        self.skyEffTime = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)] 
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                #print iRow,iCol,
                count = self.skyFile.getPixelCount(iRow,iCol)
                #print count
                #self.skySpectra[iRow][iCol],self.wvlBinEdges = self.skyFile.getPixelSpectrum(iRow,iCol,weighted=True,firstSec=0,integrationTime=-1)
                skyDict = self.skyFile.getPixelSpectrum(iRow,iCol,weighted=True,firstSec=0,integrationTime=-1)
                self.skySpectra[iRow][iCol],self.wvlBinEdges,self.skyEffTime[iRow][iCol] = skyDict['spectrum'],skyDict['wvlBinEdges'],skyDict['effIntTime']
        self.skySpectra = np.array(self.skySpectra)
        self.skyEffTime = np.array(self.skyEffTime)
        self.skySpectra = self.skySpectra/self.binWidths
        return self.skySpectra

    def loadSpectrum(self, objectName='58Aquilae'):
        #import the known spectrum of the calibrator and rebin to the histogram parameters given
        #must be imported into array with dtype float so division later does not have error
        std = MKIDStd.MKIDStd()
        a = std.load(objectName)
        #a = std.normalizeFlux(a)
        x = a[:,0]
        y = a[:,1]
        self.realSpectra = y
        self.realSpectraWvl = x

        x,y = self.cleanSpectrum(self.realSpectraWvl,self.realSpectra,objectName)

        newa = rebin(x,y,self.wvlBinEdges)
        self.binnedSpectra = newa[:,1]

    def cleanSpectrum(self,x,y,objectName):
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

    def calculateFactors(self):
        """
        Calculate the sensitivity spectrum: the weighting factors that correct the flat calibrated spectra to the real spectra
        First subtract sky spectrum from ARCONS observed spectrum. Then divide this into the known spectrum of the object.
        Take the median of these factors, since the flat cal should have removed any pixel-pixel variation. Ideally we only need one master
        set of flux cal weights for the whole array.
        """
        #scale skySpectra to same exposure time as fluxSpectra
        print self.fluxEffTime
        print self.skyEffTime
        intTimeScaled = self.fluxEffTime/self.skyEffTime
        print intTimeScaled
        self.scaledSkySpectra = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                print self.skySpectra[iRow][iCol]
                print intTimeScaled[iRow][iCol]
                self.scaledSkySpectra[iRow][iCol] = self.skySpectra[iRow][iCol] * intTimeScaled[iRow][iCol]
                print self.scaledSkySpectra[iRow][iCol]
        self.scaledSkySpectra = np.array(self.scaledSkySpectra)
        self.subtractedSpectra = self.fluxSpectra - self.scaledSkySpectra
        self.subtractedSpectra = np.array(self.subtractedSpectra,dtype=float) #cast as floats so division does not fail later
        self.subtractedSpectra = self.calculateMedian(self.subtractedSpectra) #find median of subtracted spectra across whole array
        
        normWvl = 5500 #Angstroms
        ind = np.where(self.wvlBinEdges >= normWvl)[0][0]-1
        self.subtractedSpectra = self.subtractedSpectra/(self.subtractedSpectra[ind]) #normalize
        self.binnedSpectra = self.binnedSpectra/(self.binnedSpectra[ind])

        #Calculate FluxCal factors
        self.fluxFactors = self.binnedSpectra/self.subtractedSpectra

        #self.fluxFlags = np.zeros(np.shape(self.fluxFactors),dtype='int')
        self.fluxFlags = np.empty(np.shape(self.fluxFactors),dtype='int')
        self.fluxFlags.fill(pipelineFlags.fluxCal['good'])   #Initialise flag array filled with 'good' flags. JvE 5/1/2013.
        #set factors that will cause trouble to 1
        #self.fluxFlags[self.fluxFactors == np.inf] = 1
        self.fluxFlags[self.fluxFactors == np.inf] = pipelineFlags.fluxCal['infWeight']   #Modified to use flag dictionary - JvE 5/1/2013
        self.fluxFactors[self.fluxFactors == np.inf]=1.0
        self.fluxFactors[np.isnan(self.fluxFactors)]=1.0
        self.fluxFlags[np.isnan(self.fluxFactors)] = pipelineFlags.fluxCal['nanWeight']   #Modified to use flag dictionary - JvE 5/1/2013
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

    def makePlots(self, obsfile, skyfile, target):
        """
        Output all debugging plots of ARCONS sky and object spectra, known calibrator spectrum, and sensitivity curve
        """
        scratchDir = os.getenv('MKID_PROC_PATH')
        fluxDir = os.path.join(scratchDir,'fluxCalSolnFiles')
        fluxCalBase = obsfile.split("/")[-1].split(".")[0]        
        plotFileName = fluxCalBase+".pdf"
        fullFluxPlotFileName = os.path.join(fluxDir,plotFileName)
        pp = PdfPages(fullFluxPlotFileName)
        matplotlib.rcParams['font.size']=6

        #calculate midpoints of wvl bins for plotting
        wvls = np.empty((self.nWvlBins),dtype=float)
        for n in xrange(self.nWvlBins):
            binsize=self.wvlBinEdges[n+1]-self.wvlBinEdges[n]
            wvls[n] = (self.wvlBinEdges[n]+(binsize/2.0))

        medianFluxData = self.calculateMedian(self.fluxSpectra)
        medianSkyData = self.calculateMedian(self.scaledSkySpectra)
        
        plt.figure()

        ax1 = plt.subplot(331)
        ax1.set_title('ARCONS median flat cal\'d flux in counts')
        plt.plot(wvls,medianFluxData)
        #plt.show()
        ax2 = plt.subplot(332)
        ax2.set_title('ARCONS median flat cal\'d sky in counts')
        plt.plot(wvls,medianSkyData)
        #plt.show()
        ax3 = plt.subplot(333)
        ax3.set_title('Flux data minus sky in counts')
        plt.plot(wvls,self.subtractedSpectra)
        ax4 = plt.subplot(334)
        ax4.set_title('Std Spectrum of %s'%(target))
        plt.plot(self.realSpectraWvl,self.realSpectra)
        ax5 = plt.subplot(335)
        ax5.set_title('Binned Std Spectrum')
        plt.plot(wvls,self.binnedSpectra)
        ax6 = plt.subplot(336)
        ax6.set_title('Median Sensitivity Spectrum')
        ax6.set_xlim((3000,13000))
        ax6.set_ylim((0,5))
        plt.plot(wvls,self.fluxFactors)
        ax7 = plt.subplot(337)
        ax7.set_title('Flux Cal\'d ARCONS Spectrum of Std')
        plt.plot(wvls,self.fluxFactors*self.subtractedSpectra)

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

if __name__ == '__main__':

    #paramFile = sys.argv[1]

    if len(sys.argv < 7):
        print "Bad arguments given \n Syntax is: FluxCal(fluxFileName, skyFileName, wvlCalFileName,flatCalFileName, objectName, fluxCalFileName, plots=True) \n note, flux and sky files require the paths, fluxCalFile automatically outputs to INTERM directory"

    fluxFileName=sys.argv[1]
    skyFileName=sys.argv[2]
    wvlCalFileName=sys.argv[3]
    flatCalFileName=sys.argv[4]
    objectName=sys.argv[5]
    fluxCalFileName=sys.argv[6]
    plots=sys.argv[7]

    fc = FluxCal(fluxFileName, skyFileName, wvlCalFileName,flatCalFileName, objectName, fluxCalFileName, plots=True)





































