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
import glob
from util.ObsFile import ObsFile
from util.readDict import readDict
from util.FileName import FileName
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

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
        wvlSunsetDate = self.params['wvlSunsetDate']
        wvlTimestamp = self.params['wvlTimestamp']
        needTimeAdjust = self.params['needTimeAdjust']
        needHotPix = self.params['needHotPix']
        obsSequence = self.params['obsSequence']

        #flatFileName = params['flatFileName']
        #wvlCalFileName = params['wvlCalFileName']
        #flatCalFileName = params['flatCalFileName']
        obsFNs = [FileName(run=run,date=sunsetDate,tstamp=obsTstamp) for obsTstamp in obsSequence]
        self.obsFileNames = [fn.obs() for fn in obsFNs]
        self.obsList = [ObsFile(obsFileName) for obsFileName in self.obsFileNames]
        timeMaskFileNames = [fn.timeMask() for fn in obsFNs]
        timeAdjustFileName = FileName(run=run).timeAdjustments()

        print len(self.obsFileNames), 'flat files to co-add'
        self.flatCalFileName = FileName(run=run,date=sunsetDate).flatSoln()
        wvlCalFileName = FileName(run=run,date=wvlSunsetDate,tstamp=wvlTimestamp).calSoln()
        for iObs,obs in enumerate(self.obsList):
           obs.loadWvlCalFile(wvlCalFileName)
           obs.loadTimeAdjustmentFile(timeAdjustFileName)
           obs.loadHotPixCalFile(timeMaskFileNames[iObs])

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
        #wvlBinEdges includes both lower and upper limits, so number of bins is 1 less than number of edges
        self.nWvlBins = len(self.wvlBinEdges)-1

        #print 'wrote to',self.flatCalFileName

    def __del__(self):
        pass

    def loadFlatSpectra(self):
        self.spectralCubes = []
        self.frames = []
        for iObs,obs in enumerate(self.obsList):
            print 'obs',iObs
            for firstSec in range(0,obs.getFromHeader('exptime'),self.intTime):
                print 'sec',firstSec
                cubeDict = obs.getSpectralCube(firstSec=firstSec,integrationTime=self.intTime,weighted=False,wvlBinEdges = self.wvlBinEdges)
                cube = np.array(cubeDict['cube'],dtype=np.double)
                effIntTime = cubeDict['effIntTime']
                #add third dimension for broadcasting
                effIntTime = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
                cube /= effIntTime
                frame = np.sum(cube,axis=2)
                frame[np.isnan(frame)]=0
                self.frames.append(frame)
                #self.plotArray(frame)
                self.spectralCubes.append(cube)
        self.spectralCubes = np.array(self.spectralCubes)
        #self.plotArray(self.frames[0])
        print 'done load'

    def checkCountRates(self):
        print 'check'
        medianCountRates = np.array([np.median(frame[frame!=0]) for frame in self.frames])
        boolIncludeFrames = medianCountRates <= self.countRateCutoff
        #mask out frames, or cubes from integration time chunks with count rates too high
        self.spectralCubes = np.array([cube for cube,boolIncludeFrame in zip(self.spectralCubes,boolIncludeFrames) if boolIncludeFrame==True])
        self.frames = [frame for frame,boolIncludeFrame in zip(self.frames,boolIncludeFrames) if boolIncludeFrame==True]
        print boolIncludeFrames

    def calculateWeights(self):
        """
        finds flat cal factors as medians/pixelSpectra for each pixel
        """
        cubeWeightsList = []
        self.medianSpectra = []
        for iCube,cube in enumerate(self.spectralCubes):
            wvlMedians = np.zeros(self.nWvlBins)
            spectra2d = np.reshape(cube,[self.nRow*self.nCol,self.nWvlBins ])
            for iWvl in xrange(self.nWvlBins):
                spectrum = spectra2d[:,iWvl]
                goodSpectrum = spectrum[spectrum != 0]#dead pixels need to be taken out before calculating medians
                wvlMedians[iWvl] = np.median(goodSpectrum)
            weights = np.divide(wvlMedians,cube)
            weights[weights==0] = np.nan
            weights[weights==np.inf] = np.nan
            cubeWeightsList.append(weights)
            self.medianSpectra.append(wvlMedians)
        cubeWeights = np.array(cubeWeightsList)
        cubeWeightsMask = np.isnan(cubeWeights)
        self.maskedCubeWeights = np.ma.array(cubeWeights,mask=cubeWeightsMask,fill_value=1.)

        #sort maskedCubeWeights and rearange spectral cubes the same way
        sortedIndices = np.ma.argsort(self.maskedCubeWeights,axis=0)
        identityIndices = np.ma.indices(np.shape(self.maskedCubeWeights))

        sortedWeights = self.maskedCubeWeights[sortedIndices,identityIndices[1],identityIndices[2],identityIndices[3]]
        spectralCubesReordered = self.spectralCubes[sortedIndices,identityIndices[1],identityIndices[2],identityIndices[3]]

        #trim the beginning and end off the sorted weights for each wvl for each pixel, to exclude extriems from averages
        nCubes = np.shape(self.maskedCubeWeights)[0]
        fractionToRemove=.15 #off both top and bottom
        trimmedWeights = sortedWeights[fractionToRemove*nCubes:(1-fractionToRemove)*nCubes,:,:,:]
        trimmedSpectralCubesReordered = spectralCubesReordered[fractionToRemove*nCubes:(1-fractionToRemove)*nCubes,:,:,:]

        #self.flatWeights = np.ma.average(self.maskedCubeWeights,axis=0,weights=self.spectralCubes)
        self.flatWeights = np.ma.average(trimmedWeights,axis=0,weights=trimmedSpectralCubesReordered)
        self.flatFlags = self.flatWeights.mask
        flagImage = np.shape(self.flatFlags)[2]-np.sum(self.flatFlags,axis=2)
        #self.plotArray(flagImage)

#        X,Y,Z=np.mgrid[0:self.nRow,0:self.nCol,0:self.nWvlBins]
#        Z=self.wvlBinEdges[Z]
#        fig = plt.figure()
#        ax = Axes3D(fig)
#        handleScatter=ax.scatter(X,Y,Z,c=self.flatWeights,vmax=2,vmin=.5)
#        fig.colorbar(handleScatter)
#        plt.show()
        
    def plotWeights(self):
        print 'plotWeights'
        flatCalPath,flatCalBasename = os.path.split(self.flatCalFileName)
        pdfBasename = os.path.splitext(flatCalBasename)[0]+'.pdf'
        pp = PdfPages(os.path.join(flatCalPath,pdfBasename))
        nPlotsPerRow = 2
        nPlotsPerCol = 4
        nPlotsPerPage = nPlotsPerRow*nPlotsPerCol
        iPlot = 0 

        matplotlib.rcParams['font.size'] = 4 
        wvls = self.wvlBinEdges[0:-1]
        nCubes = len(self.maskedCubeWeights)

        for iRow in xrange(self.nRow):
            print iRow
            for iCol in xrange(self.nCol):
                weights = self.flatWeights[iRow,iCol,:]
                if weights.mask.all() == False:
                    if iPlot % nPlotsPerPage == 0:
                        fig = plt.figure(figsize=(10,10),dpi=100)

                    ax = fig.add_subplot(nPlotsPerCol,nPlotsPerRow,iPlot%nPlotsPerPage+1)
                    ax.set_ylim(.5,1.5)

                    for iCube in range(nCubes):
                        cubeWeights = self.maskedCubeWeights[iCube,iRow,iCol]
                        ax.plot(wvls,cubeWeights.data,label='weights %d'%iCube,alpha=.7,color=matplotlib.cm.Paired((iCube+1.)/nCubes))
                    ax.plot(wvls,weights.data,label='weights',color='k')
                
                    ax.set_title('p %d,%d'%(iRow,iCol))
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
                    for iCube in range(nCubes):
                        spectrum = self.spectralCubes[iCube,iRow,iCol]
                        ax.plot(wvls,spectrum,label='spectrum %d'%iCube,alpha=.7,color=matplotlib.cm.Paired((iCube+1.)/nCubes))
                
                    ax.set_title('p %d,%d'%(iRow,iCol))
                    #ax.plot(wvls,flatSpectrum,label='pixel',alpha=.5)

                    #ax.legend(loc='lower left')
                    #ax2.legend(loc='lower right')
                    if iPlot%nPlotsPerPage == nPlotsPerPage-1 or (iRow == self.nRow-1 and iCol == self.nCol-1):
                        pp.savefig(fig)
                        #plt.show()
                    iPlot += 1
        pp.close()
        print 'done plotWeights'

    

    def writeWeights(self):
        """
        Writes an h5 file to put calculated flat cal factors in
        """
        print 'write'
        if os.path.isabs(self.flatCalFileName) == True:
            fullFlatCalFileName = self.flatCalFileName
        else:
            scratchDir = os.getenv('INTERM_PATH')
            flatDir = os.path.join(scratchDir,'flatCalSolnFiles')
            fullFlatCalFileName = os.path.join(flatDir,self.flatCalFileName)

        try:
            flatCalFile = tables.openFile(fullFlatCalFileName,mode='w')
        except:
            print 'Error: Couldn\'t create flat cal file, ',fullFlatCalFileName
            return
        print 'wrote to',self.flatCalFileName

        calgroup = flatCalFile.createGroup(flatCalFile.root,'flatcal','Table of flat calibration weights by pixel and wavelength')
        caltable = tables.Array(calgroup,'weights',object=self.flatWeights.data,title='Flat calibration Weights indexed by pixelRow,pixelCol,wavelengthBin')
        flagtable = tables.Array(calgroup,'flags',object=self.flatFlags,title='Flat cal flags indexed by pixelRow,pixelCol,wavelengthBin. 0 is Good')
        bintable = tables.Array(calgroup,'wavelengthBins',object=self.wvlBinEdges,title='Wavelength bin edges corresponding to third dimension of weights array')
        flatCalFile.flush()
        flatCalFile.close()

        npzFileName = os.path.splitext(fullFlatCalFileName)[0]+'.npz'

        #calculate total spectra and medians for programs that expect old format flat cal
        spectra = np.array(np.sum(self.spectralCubes,axis=0))

        wvlMedians = np.zeros(self.nWvlBins)
        spectra2d = np.reshape(spectra,[self.nRow*self.nCol,self.nWvlBins ])
        for iWvl in xrange(self.nWvlBins):
            spectrum = spectra2d[:,iWvl]
            goodSpectrum = spectrum[spectrum != 0]#dead pixels need to be taken out before calculating medians
            wvlMedians[iWvl] = np.median(goodSpectrum)
        np.savez(npzFileName,median=wvlMedians,medianSpectra=np.array(self.medianSpectra),binEdges=self.wvlBinEdges,spectra=spectra,weights=np.array(self.flatWeights.data))

    def plotArray(self,image,normNSigma=3,title=''):
        self.fig = plt.figure()
        self.axes = self.fig.add_subplot(111)
        handleMatshow = self.axes.matshow(image,cmap=matplotlib.cm.gnuplot2,origin='lower',vmax=np.mean(image)+normNSigma*np.std(image))
        self.fig.cbar = self.fig.colorbar(handleMatshow)
        self.axes.set_title(title)
        cid = self.fig.canvas.mpl_connect('scroll_event', partial(onscroll_cbar, self.fig))
        cid = self.fig.canvas.mpl_connect('button_press_event', partial(onclick_cbar, self.fig))
        plt.show()

if __name__ == '__main__':
    paramFile = sys.argv[1]
    flatcal = FlatCal(paramFile)
    flatcal.loadFlatSpectra()
    flatcal.checkCountRates()
    flatcal.calculateWeights()
    flatcal.plotWeights()
    flatcal.writeWeights()

