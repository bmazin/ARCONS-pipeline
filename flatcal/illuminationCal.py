#!/bin/python
"""
Author: Matt Strader        Date:August 19,2012
Opens a twilight flat h5 and makes the spectrum of each pixel.
Then takes the median of each energy over all pixels
A factor is then calculated for each energy in each pixel of its
twilight count rate / median count rate
The factors are written out in an h5 file
"""

import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from functools import partial
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages


from util.popup import PopUp,plotArray,pop
from util.ObsFile import ObsFile
from util.readDict import readDict
from util.FileName import FileName
from util.utils import nearestNRobustMeanFilter
import hotpix.hotPixels as hp

def onscroll_cbar(fig, event):
    if event.inaxes is fig.cbar.ax:
        increment=0.05
        currentClim = fig.cbar.mappable.get_clim()
        if event.button == 'up':
            newClim = (currentClim[0],(1.+increment)*currentClim[1])
        if event.button == 'down':
            newClim = (currentClim[0],(1.-increment)*currentClim[1])
        fig.cbar.mappable.set_clim(newClim)
        fig.canvas.draw()

def onclick_cbar(fig,event):
    if event.inaxes is fig.cbar.ax:
        if event.button == 1:
            fig.oldClim = fig.cbar.mappable.get_clim()
            fig.cbar.mappable.set_clim(fig.oldClim[0],event.ydata*fig.oldClim[1])
            fig.canvas.draw()
        if event.button == 3:
            fig.oldClim = fig.cbar.mappable.get_clim()
            fig.cbar.mappable.set_clim(fig.oldClim[0],1/event.ydata*fig.oldClim[1])
            fig.canvas.draw()

class FlatCal:
    def __init__(self,paramFile):
        """
        opens flat file,sets wavelength binnning parameters, and calculates flat factors for the file
        """
        self.params = readDict()
        self.params.read_from_file(paramFile)

        run = self.params['run']
        sunsetDate = self.params['sunsetDate']
        flatTstamp = self.params['flatTstamp']
        wvlSunsetDate = self.params['wvlSunsetDate']
        wvlTimestamp = self.params['wvlTimestamp']
        obsSequence = self.params['obsSequence']
        self.deadtime = self.params['deadtime'] #from firmware pulse detection
        self.intTime = self.params['intTime']
        self.timeSpacingCut = self.params['timeSpacingCut']
        self.nSigmaClip = self.params['nSigmaClip']
        self.nNearest = self.params['nNearest']

        obsFNs = [FileName(run=run,date=sunsetDate,tstamp=obsTstamp) for obsTstamp in obsSequence]
        self.obsFileNames = [fn.obs() for fn in obsFNs]
        self.obsList = [ObsFile(obsFileName) for obsFileName in self.obsFileNames]
        timeMaskFileNames = [fn.timeMask() for fn in obsFNs]
        timeAdjustFileName = FileName(run=run).timeAdjustments()

        print len(self.obsFileNames), 'flat files to co-add'
        self.flatCalFileName = FileName(run=run,date=sunsetDate,tstamp=flatTstamp).illumSoln()
        if wvlSunsetDate != '':
            wvlCalFileName = FileName(run=run,date=wvlSunsetDate,tstamp=wvlTimestamp).calSoln()
        for iObs,obs in enumerate(self.obsList):
            if wvlSunsetDate != '':
                obs.loadWvlCalFile(wvlCalFileName)
            else:
                obs.loadBestWvlCalFile()
            obs.loadTimeAdjustmentFile(timeAdjustFileName)
            timeMaskFileName = timeMaskFileNames[iObs]
            print timeMaskFileName
            #Temporary step, remove old hotpix file
            #if os.path.exists(timeMaskFileName):
            #    os.remove(timeMaskFileName)
            if not os.path.exists(timeMaskFileName):
                print 'Running hotpix for ',obs
                hp.findHotPixels(self.obsFileNames[iObs],timeMaskFileName,fwhm=np.inf,useLocalStdDev=True)
                print "Flux file pixel mask saved to %s"%(timeMaskFileName)
            obs.loadHotPixCalFile(timeMaskFileName)
        self.wvlFlags = self.obsList[0].wvlFlagTable

        self.nRow = self.obsList[0].nRow
        self.nCol = self.obsList[0].nCol
        print 'files opened'
        #self.wvlBinWidth = params['wvlBinWidth'] #angstroms
        self.energyBinWidth = self.params['energyBinWidth'] #eV
        self.wvlStart = self.params['wvlStart'] #angstroms
        self.wvlStop = self.params['wvlStop'] #angstroms
        self.wvlBinEdges = ObsFile.makeWvlBins(self.energyBinWidth,self.wvlStart,self.wvlStop)
        self.intTime = self.params['intTime']
        self.countRateCutoff = self.params['countRateCutoff']
        self.fractionOfChunksToTrim = self.params['fractionOfChunksToTrim']
        #wvlBinEdges includes both lower and upper limits, so number of bins is 1 less than number of edges
        self.nWvlBins = len(self.wvlBinEdges)-1

        #print 'wrote to',self.flatCalFileName

    def __del__(self):
        pass

    def loadFlatSpectra(self):
        self.spectralCubes = []#each element will be the spectral cube for a time chunk
        self.cubeEffIntTimes = []
        self.frames = []
        for iObs,obs in enumerate(self.obsList):
            print 'obs',iObs
            for firstSec in range(0,obs.getFromHeader('exptime'),self.intTime):
                print 'sec',firstSec
                cubeDict = obs.getSpectralCube(firstSec=firstSec,integrationTime=self.intTime,weighted=False,wvlBinEdges = self.wvlBinEdges,timeSpacingCut = self.timeSpacingCut)
                cube = np.array(cubeDict['cube'],dtype=np.double)
                effIntTime = cubeDict['effIntTime']
                #add third dimension for broadcasting
                effIntTime3d = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
                cube /= effIntTime3d
                
                
                cube[np.isnan(cube)]=0
                frame = np.sum(cube,axis=2) #in counts per sec
                #correct nonlinearity due to deadtime in firmware
                nonlinearFactors = 1. / (1. - frame*self.deadtime)
                nonlinearFactors[np.isnan(nonlinearFactors)]=0.
                frame = frame * nonlinearFactors
                
                nonlinearFactors = np.reshape(nonlinearFactors,np.shape(nonlinearFactors)+(1,))
                cube = cube * nonlinearFactors

                
                self.frames.append(frame)
                self.spectralCubes.append(cube)
                self.cubeEffIntTimes.append(effIntTime3d)
        self.spectralCubes = np.array(self.spectralCubes)
        self.cubeEffIntTimes = np.array(self.cubeEffIntTimes)
        self.countCubes = self.cubeEffIntTimes * self.spectralCubes
        self.spectralCubes = self.intTime * self.spectralCubes # in counts

    def checkCountRates(self):
        medianCountRates = np.array([np.median(frame[frame!=0]) for frame in self.frames])
        boolIncludeFrames = medianCountRates <= self.countRateCutoff
        #boolIncludeFrames = np.logical_and(boolIncludeFrames,medianCountRates >= 200) 
        #mask out frames, or cubes from integration time chunks with count rates too high
        self.spectralCubes = np.array([cube for cube,boolIncludeFrame in zip(self.spectralCubes,boolIncludeFrames) if boolIncludeFrame==True])
        self.frames = [frame for frame,boolIncludeFrame in zip(self.frames,boolIncludeFrames) if boolIncludeFrame==True]
        print 'few enough counts in the chunk',zip(medianCountRates,boolIncludeFrames)

    def calculateWeights(self):
        """
        finds illum cal factors by making an image per wavelength bin, then smoothing it and dividing the mean of the smoothedImage over pixels by the smoothedImage
        """
        cubeWeightsList = []
        self.averageSpectra = []
        deltaWeightsList = []
        self.totalCube = np.sum(self.spectralCubes,axis=0) #sum all cubes
        self.totalFrame = np.sum(self.totalCube,axis=-1)#sum over wvl
        weights = []
        for iWvl in xrange(self.nWvlBins):
            wvlSlice = self.totalCube[:,:,iWvl]
            wvlSlice[wvlSlice == 0] = np.nan
            smoothedWvlSlice = nearestNRobustMeanFilter(wvlSlice,n=self.nNearest,nSigmaClip=self.nSigmaClip)
            wvlIllumWeights = np.mean(smoothedWvlSlice)/smoothedWvlSlice
            weights.append(wvlIllumWeights)
        self.weights = np.array(weights)
        #move the wvl dimension to the end
        self.weights = np.swapaxes(self.weights,0,1)
        self.weights = np.swapaxes(self.weights,1,2)
        self.deltaWeights = np.zeros_like(self.weights)
        self.flags = np.zeros_like(self.weights)
        
            
        
    def plotWeightsWvlSlices(self,verbose=True):
        flatCalPath,flatCalBasename = os.path.split(self.flatCalFileName)
        pdfBasename = os.path.splitext(flatCalBasename)[0]+'_wvlSlices.pdf'
        pdfFullPath = os.path.join(flatCalPath,pdfBasename)
        pp = PdfPages(pdfFullPath)
        nPlotsPerRow = 2
        nPlotsPerCol = 4 
        nPlotsPerPage = nPlotsPerRow*nPlotsPerCol
        iPlot = 0 
        if verbose:
            print 'plotting weights in wavelength sliced images'

        matplotlib.rcParams['font.size'] = 4 
        wvls = self.wvlBinEdges[0:-1]

        cmap = matplotlib.cm.gnuplot2
        cmap.set_bad('0.15')
        for iWvl,wvl in enumerate(wvls):
            if verbose:
                print 'wvl ',iWvl
            if iPlot % nPlotsPerPage == 0:
                fig = plt.figure(figsize=(10,10),dpi=100)

            ax = fig.add_subplot(nPlotsPerCol,nPlotsPerRow,iPlot%nPlotsPerPage+1)
            ax.set_title(r'Weights %.0f $\AA$'%wvl)

            image = self.weights[:,:,iWvl]

            handleMatshow = ax.matshow(image,cmap=cmap,origin='lower',vmax=2.)
            cbar = fig.colorbar(handleMatshow)
        
            if iPlot%nPlotsPerPage == nPlotsPerPage-1:
                pp.savefig(fig)
            iPlot += 1


            if iPlot % nPlotsPerPage == 0:
                fig = plt.figure(figsize=(10,10),dpi=100)

            ax = fig.add_subplot(nPlotsPerCol,nPlotsPerRow,iPlot%nPlotsPerPage+1)
            ax.set_title(r'Twilight Image %.0f $\AA$'%wvl)

            image = self.totalCube[:,:,iWvl]

            nSdev = 3.
            goodImage = image[np.isfinite(image)]
            vmax = np.mean(goodImage)+nSdev*np.std(goodImage)
            handleMatshow = ax.matshow(image,cmap=cmap,origin='lower',vmax=vmax)
            cbar = fig.colorbar(handleMatshow)
        
            if iPlot%nPlotsPerPage == nPlotsPerPage-1:
                pp.savefig(fig)
            iPlot += 1

        pp.savefig(fig)
        pp.close()


    def plotWeightsByPixel(self,verbose=True):
        flatCalPath,flatCalBasename = os.path.split(self.flatCalFileName)
        pdfBasename = os.path.splitext(flatCalBasename)[0]+'.pdf'
        pdfFullPath = os.path.join(flatCalPath,pdfBasename)
        pp = PdfPages(pdfFullPath)
        nPlotsPerRow = 2
        nPlotsPerCol = 4
        nPlotsPerPage = nPlotsPerRow*nPlotsPerCol
        iPlot = 0 
        if verbose:
            print 'plotting weights by pixel at ',pdfFullPath

        matplotlib.rcParams['font.size'] = 4 
        wvls = self.wvlBinEdges[0:-1]
        nCubes = len(self.spectralCubes)

        for iRow in xrange(self.nRow):
            if verbose:
                print 'row',iRow
            for iCol in xrange(self.nCol):
                weights = self.weights[iRow,iCol,:]
                deltaWeights = self.deltaWeights[iRow,iCol,:]
                if iPlot % nPlotsPerPage == 0:
                    fig = plt.figure(figsize=(10,10),dpi=100)

                ax = fig.add_subplot(nPlotsPerCol,nPlotsPerRow,iPlot%nPlotsPerPage+1)
                ax.set_ylim(.5,2.)

                weights = self.weights[iRow,iCol]
                ax.errorbar(wvls,weights,yerr=deltaWeights,label='weights',color='k')
            
                ax.set_title('p %d,%d'%(iRow,iCol))
                ax.set_ylabel('weight')
                ax.set_xlabel(r'$\lambda$ ($\AA$)')
                #ax.plot(wvls,flatSpectrum,label='pixel',alpha=.5)

                #ax.legend(loc='lower left')
                #ax2.legend(loc='lower right')
                if iPlot%nPlotsPerPage == nPlotsPerPage-1 or (iRow == self.nRow-1 and iCol == self.nCol-1):
                    pp.savefig(fig)
                iPlot += 1

                #Put a plot of twilight spectrums for this pixel
                if iPlot % nPlotsPerPage == 0:
                    fig = plt.figure(figsize=(10,10),dpi=100)

                ax = fig.add_subplot(nPlotsPerCol,nPlotsPerRow,iPlot%nPlotsPerPage+1)
                ax.plot(wvls,self.totalCube[iRow,iCol,:],label='spectrum')
            
                ax.set_title('p %d,%d'%(iRow,iCol))
                ax.set_xlabel(r'$\lambda$ ($\AA$)')
                ax.set_ylabel('twilight cps')
                #ax.plot(wvls,flatSpectrum,label='pixel',alpha=.5)

                #ax.legend(loc='lower left')
                #ax2.legend(loc='lower right')
                if iPlot%nPlotsPerPage == nPlotsPerPage-1 or (iRow == self.nRow-1 and iCol == self.nCol-1):
                    pp.savefig(fig)
                    #plt.show()
                iPlot += 1
        pp.close()

    

    def writeWeights(self):
        """
        Writes an h5 file to put calculated flat cal factors in
        """
        if os.path.isabs(self.flatCalFileName) == True:
            fullFlatCalFileName = self.flatCalFileName
        else:
            scratchDir = os.getenv('MKID_PROC_PATH')
            flatDir = os.path.join(scratchDir,'flatCalSolnFiles')
            fullFlatCalFileName = os.path.join(flatDir,self.flatCalFileName)

        try:
            flatCalFile = tables.openFile(fullFlatCalFileName,mode='w')
        except:
            print 'Error: Couldn\'t create flat cal file, ',fullFlatCalFileName
            return
        print 'wrote to',self.flatCalFileName

        calgroup = flatCalFile.createGroup(flatCalFile.root,'flatcal','Table of flat calibration weights by pixel and wavelength')
        caltable = tables.Array(calgroup,'weights',object=self.weights,title='Illumination calibration Weights indexed by pixelRow,pixelCol,wavelengthBin')
        errtable = tables.Array(calgroup,'errors',object=self.deltaWeights,title='Errors in Weights indexed by pixelRow,pixelCol,wavelengthBin')
        flagtable = tables.Array(calgroup,'flags',object=self.flags,title='Illumination cal flags indexed by pixelRow,pixelCol,wavelengthBin. 0 is Good. By default, all are good.')
        bintable = tables.Array(calgroup,'wavelengthBins',object=self.wvlBinEdges,title='Wavelength bin edges corresponding to third dimension of weights array')
        flatCalFile.flush()
        flatCalFile.close()

        npzFileName = os.path.splitext(fullFlatCalFileName)[0]+'.npz'

        #calculate total spectra and medians for programs that expect old format flat cal
        spectra = np.array(np.sum(self.spectralCubes,axis=0))

        wvlAverages = np.zeros(self.nWvlBins)
        spectra2d = np.reshape(spectra,[self.nRow*self.nCol,self.nWvlBins ])
        np.savez(npzFileName,binEdges=self.wvlBinEdges,spectra=spectra,weights=self.weights,deltaWeights=self.deltaWeights,totalFrame=self.totalFrame,totalCube=self.totalCube,spectralCubes=self.spectralCubes,countCubes=self.countCubes,cubeEffIntTimes=self.cubeEffIntTimes )


if __name__ == '__main__':
    paramFile = sys.argv[1]
    flatcal = FlatCal(paramFile)
    flatcal.loadFlatSpectra()
    flatcal.checkCountRates()
    flatcal.calculateWeights()
    flatcal.writeWeights()
    flatcal.plotWeightsWvlSlices()
    flatcal.plotWeightsByPixel()

