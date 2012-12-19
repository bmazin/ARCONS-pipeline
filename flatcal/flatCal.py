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
        flatFileName = FileName(run=run,date=flatSunsetLocalDate,tstamp=flatTimestamp).flat()
        print flatFileName
        flatCalFileName = FileName(run=run,date=flatSunsetLocalDate,tstamp=flatTimestamp).flatSoln()
        print flatCalFileName
        wvlCalFileName = FileName(run=run,date=wvlCalSunsetLocalDate,tstamp=wvlCalTimestamp).calSoln()


        self.wvlBinWidth = params['wvlBinWidth'] #angstroms
        self.wvlStart = params['wvlStart'] #angstroms
        self.wvlStop = params['wvlStop'] #angstroms
        self.nWvlBins = int((self.wvlStop - self.wvlStart)/self.wvlBinWidth)
        
        self.loadFlatFile(flatFileName,wvlCalFileName)
        self.loadFlatSpectra()
        self.calculateFactors()
        self.writeFactors(flatCalFileName)

    def __del__(self):
        try:
            self.flatFile.close()
            self.calFile.close()
        except AttributeError:#flatFile was never defined
            pass

    def calculateFactors(self):
        """
        finds flat cal factors as medians/pixelSpectra for each pixel
        """
        self.flatFactors = np.divide(self.wvlMedians,self.spectra)
        #set factors that will cause trouble to 1
        self.flatFactors[self.flatFactors == np.inf]=1.0
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
    def loadFlatFile(self,flatFileName,wvlCalFileName):
        #open the hdf5 file
        self.flatFile = ObsFile(flatFileName)
        self.nRow = self.flatFile.nRow
        self.nCol = self.flatFile.nCol
        self.flatFile.loadWvlCalFile(wvlCalFileName)#load the wavelength cal to be applied to the flat
        
    def loadFlatSpectra(self):
        self.spectra = [[[] for i in xrange(self.nCol)] for j in xrange(self.nRow)]
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                count = self.flatFile.getPixelCount(iRow,iCol)
                self.spectra[iRow][iCol],self.wvlBinEdges = self.flatFile.getPixelSpectrum(iRow,iCol,wvlStart=self.wvlStart,wvlStop=self.wvlStop,wvlBinWidth=self.wvlBinWidth,weighted=False,firstSec=0,integrationTime=-1)
#                if iRow == 45 and iCol == 12:
#                    print iRow,iCol,count,sum(self.spectra[iRow][iCol])
        self.spectra = np.array(self.spectra)
        self.calculateMedians()
        return self.spectra

    def writeFactors(self,flatCalFileName):
        """
        Writes an h5 file to put calculated flat cal factors in
        """
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

if __name__ == '__main__':
    paramFile = sys.argv[1]
    fc = FlatCal(paramFile)
