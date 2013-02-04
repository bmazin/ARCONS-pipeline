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
import glob
from util.ObsFile import ObsFile
from util.readDict import readDict
from util.FileName import FileName

class FlatCal:
    def __init__(self,paramFile):
        """
        opens flat file,sets wavelength binnning parameters, and calculates flat factors for the file
        """
        params = readDict()
        params.read_from_file(paramFile)

        run = params['run']
        flatSunsetLocalDate = params['flatSunsetLocalDate']
        flatTimestamp = params['flatTimestamp']

        wvlCalSunsetLocalDate = params['wvlCalSunsetLocalDate']
        wvlCalTimestamp = params['wvlCalTimestamp']

        #flatFileName = params['flatFileName']
        #wvlCalFileName = params['wvlCalFileName']
        #flatCalFileName = params['flatCalFileName']
        flatFileNames = [FileName(run=run,date=flatSunsetLocalDate,tstamp=flatTimestamp).flat()]
        if flatTimestamp == '':
            flatFilePath = os.path.dirname(flatFileNames[0])
            flatFileNames = sorted(glob.glob(os.path.join(flatFilePath,'flat*.h5')))
        print len(flatFileNames), 'flat files to co-add'
        flatCalFileName = FileName(run=run,date=flatSunsetLocalDate,tstamp=flatTimestamp).flatSoln()
        wvlCalFileName = FileName(run=run,date=wvlCalSunsetLocalDate,tstamp=wvlCalTimestamp).calSoln()


        #self.wvlBinWidth = params['wvlBinWidth'] #angstroms
        self.energyBinWidth = params['energyBinWidth'] #eV
        self.wvlStart = params['wvlStart'] #angstroms
        self.wvlStop = params['wvlStop'] #angstroms
        self.setWvlBins(self.energyBinWidth,self.wvlStart,self.wvlStop)
        #self.nWvlBins = int((self.wvlStop - self.wvlStart)/self.wvlBinWidth)
        
        self.loadFlatFile(flatFileNames,wvlCalFileName)
        self.loadFlatSpectra()
        self.calculateFactors()
        self.writeFactors(flatCalFileName)
        print 'wrote to',flatCalFileName

    def __del__(self):
        try:
            for flat in self.flatFiles:
                flat.close()
            self.calFile.close()
        except AttributeError:#flatFile was never defined
            pass

    def calculateFactors(self):
        """
        finds flat cal factors as medians/pixelSpectra for each pixel
        """
        self.flatFactors = np.divide(self.wvlMedians,self.spectra)
        self.flatFlags = np.zeros(np.shape(self.flatFactors),dtype='int')
        #set factors that will cause trouble to 1
        self.flatFlags[self.flatFactors == np.inf] = 1
        self.flatFactors[self.flatFactors == np.inf]=1.0
        self.flatFlags[self.flatFactors == 0]=2
        self.flatFactors[self.flatFactors == 0]=1.0
        #np.save('/ScienceData/intermediate/factors.npy',self.flatFactors)

    def calculateMedians(self):
        spectra2d = np.reshape(self.spectra,[self.nRow*self.nCol,self.nWvlBins ])
        self.wvlMedians = np.zeros(self.nWvlBins)
        for iWvl in xrange(self.nWvlBins):
            spectrum = spectra2d[:,iWvl]
            goodSpectrum = spectrum[spectrum != 0]#dead pixels need to be taken out before calculating medians
            self.wvlMedians[iWvl] = np.median(goodSpectrum)
            
        return self.wvlMedians
    def loadFlatFile(self,flatFileNames,wvlCalFileName):
        #open the hdf5 file
        self.flatFiles = [ObsFile(flatFileName) for flatFileName in flatFileNames]
        self.nRow = self.flatFiles[0].nRow
        self.nCol = self.flatFiles[0].nCol
        for flat in self.flatFiles:
            flat.loadWvlCalFile(wvlCalFileName)#load the wavelength cal to be applied to the flat
        
    def loadFlatSpectra(self):
        self.spectra = [[np.zeros(self.nWvlBins) for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        for iFlat,flat in enumerate(self.flatFiles):
            print 'flat ',iFlat
            for iRow in xrange(self.nRow):
                for iCol in xrange(self.nCol):
                    spectrum,binEdges = flat.getPixelSpectrum(iRow,iCol,wvlBinEdges=self.wvlBinEdges,weighted=False,firstSec=0,integrationTime=-1)
                    self.spectra[iRow][iCol] += spectrum
#                    if iFlat == len(self.flatFiles)-1:
#                        print iRow,iCol,sum(self.spectra[iRow][iCol])
        self.spectra = np.array(self.spectra)
        self.calculateMedians()
        return self.spectra

    def writeFactors(self,flatCalFileName):
        """
        Writes an h5 file to put calculated flat cal factors in
        """
        if os.path.isabs(flatCalFileName) == True:
            fullFlatCalFileName = flatCalFileName
        else:
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
        flagtable = tables.Array(calgroup,'flags',object=self.flatFlags,title='Flat cal flags indexed by pixelRow,pixelCol,wavelengthBin. 0 is Good')
        bintable = tables.Array(calgroup,'wavelengthBins',object=self.wvlBinEdges,title='Wavelength bin edges corresponding to third dimension of weights array')
        flatCalFile.flush()
        flatCalFile.close()

        npzFileName = os.path.splitext(fullFlatCalFileName)[0]+'.npz'
        np.savez(npzFileName,median=self.wvlMedians,binEdges=self.wvlBinEdges,spectra=self.spectra,weights=self.flatFactors)

    def setWvlBins(self,energyBinWidth,wvlStart,wvlStop):
        h = 4.135668e-15 #eV s
        c = 2.998e8 #m/s
        angstromPerMeter = 1e10
        energyStop = 1.0*h*c*angstromPerMeter/wvlStart
        energyStart = 1.0*h*c*angstromPerMeter/wvlStop
        self.nWvlBins = int((energyStop - energyStart)/self.energyBinWidth)
        energyBins = np.linspace(energyStart,energyStop,self.nWvlBins+1)
        self.wvlBinEdges = np.array(h*c*angstromPerMeter/energyBins)[::-1]

if __name__ == '__main__':
    paramFile = sys.argv[1]
    fc = FlatCal(paramFile)
