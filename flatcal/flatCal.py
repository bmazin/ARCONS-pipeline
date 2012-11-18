#!/bin/python
#Opens a twilight flat h5 and makes the spectrum of each pixel.
#Then takes the median of each energy over all pixels
#A factor is then calculated for each energy in each pixel of its
#twilight count rate / median count rate
#This is added to the flat factor in the photonList

import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile

class FlatCal:
    def __init__(self,flatFileName=None,wvlCalFileName=None,flatCalFileName=None):
        self.wvlBinWidth = 100 #angstroms
        self.wvlStart = 3000 #angstroms
        self.wvlStop = 13000 #angstroms
        self.nWvlBins = int((self.wvlStop - self.wvlStart)/self.wvlBinWidth)
        

        if flatFileName == None:
            self.loadFactors()
        else:
            self.loadFlatFile(flatFileName)
            self.flatFile.loadWvlCalFile(wvlCalFileName)
            self.loadFlatSpectra()
            self.calculateFactors()
            self.writeFactors(flatCalFileName)

    def __del__(self):
        try:
            self.flatFile.close()
            self.calFile.close()
        except AttributeError:#flatFile was never defined
            pass

    def loadFlatFile(self,flatFileName):
        #open the hdf5 file
        self.flatFile = ObsFile(flatFileName)
        self.nRow = self.flatFile.nRow
        self.nCol = self.flatFile.nCol
        
    def loadFlatSpectra(self):
        scratchDir = os.getenv('INTERM_PATH','.')
        flatDir = os.path.join(scratchDir,'flatCalSolnFiles')
        self.spectra = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                print iRow,iCol,
                count = self.flatFile.getPixelCount(iRow,iCol)
                print count
                #wvlList = self.flatFile.getPixelWvlList(iRow,iCol)
                #self.spectra[iRow][iCol],self.wvlBinEdges = np.histogram(wvlList,bins=self.nWvlBins,range=(self.wvlStart,self.wvlStop))
                self.spectra[iRow][iCol],self.wvlBinEdges = self.flatFile.getPixelSpectrum(iRow,iCol,wvlStart=self.wvlStart,wvlStop=self.wvlStop,wvlBinWidth=self.wvlBinWidth,weighted=False,firstSec=0,integrationTime=-1)
        self.spectra = np.array(self.spectra)
        np.save(os.path.join(flatDir,'flatSpectra.npy'),self.spectra)
        self.calculateMedians()
        return self.spectra


    def calculateMedians(self):
        spectra2d = np.reshape(self.spectra,[self.nRow*self.nCol,self.nWvlBins ])
        self.wvlMedians = np.zeros(self.nWvlBins)
        for iWvl in xrange(self.nWvlBins):
            print iWvl,self.wvlBinEdges[iWvl],
            spectrum = spectra2d[:,iWvl]
            goodSpectrum = spectrum[spectrum != 0]
            self.wvlMedians[iWvl] = np.median(goodSpectrum)
            print self.wvlMedians[iWvl]
            
        return self.wvlMedians
        
    def calculateFactors(self):
        fudgedSpectra = self.spectra
        #fudgedSpectra[fudgedSpectra == 0] =1.0
        self.flatFactors = np.divide(self.wvlMedians,self.spectra)
        self.flatFactors[self.flatFactors == np.inf]=1.0
        self.flatFactors[self.flatFactors == 0]=1.0
        np.save('/ScienceData/intermediate/factors.npy',self.flatFactors)

    def writeFactors(self,flatCalFileName):
        scratchDir = os.getenv('INTERM_PATH')
        flatDir = os.path.join(scratchDir,'flatCalSolnFiles')
        fullFlatCalFileName = os.path.join(flatDir,flatCalFileName)
        try:
            flatCalFile = tables.openFile(fullFlatCalFileName,mode='w')
        except:
            print 'Error: Couldn\'t create flat cal file, ',fullFlatCalFileName
            return

        calgroup = flatCalFile.createGroup(flatCalFile.root,'flatcal','Table of flat calibration weights by pixel and wavelength')
        caltable = tables.Array(calgroup,'weights',object=self.flatFactors,title='Flat calibration Weights indexed by pixelRow,pixelCol,wavelengthBin')
        bintable = tables.Array(calgroup,'wavelengthBins',object=self.wvlBinEdges,title='Wavelength bin edges corresponding to third dimension of weights array')
        flatCalFile.flush()
        flatCalFile.close()
