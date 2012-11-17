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
                wvlList = self.flatFile.getPixelWvlList(iRow,iCol)
                self.spectra[iRow][iCol],self.wvlBinEdges = np.histogram(wvlList,bins=self.nWvlBins,range=(self.wvlStart,self.wvlStop))
                print len(wvlList)
        self.spectra = np.array(self.spectra)
        np.save(os.path.join(flatDir,'flatSpectra.npy'),self.spectra)
        self.calculateMedians()
        return self.spectra

#    def applyWvlCal(self,row,col,pulseHeights):
#        wavelengths = np.array(self.calTable[row][col][0]*pulseHeights + self.calTable[row][col][1],dtype=np.uint32)
#        return wavelengths

        
    def loadSpectraFile(self):
        self.spectra = np.load('/ScienceData/intermediate/flatSpectra.npy')


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

        
    def loadFactors(self):
        self.flatFactors = np.load('/ScienceData/intermediate/factors.npy')
        #plt.matshow(self.flatFactors[:,:,4000])
        #plt.colorbar()
        #plt.title('weights at pulseh=4000')
        #plt.show()

        
    def plotMedian(self):
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(range(self.wvlBinWidth),self.wvlMedians,label='median spectrum')
        plt.xlim([3750,self.wvlBinWidth])
        plt.legend(loc=2)
        plt.show()

    def plotPixelSpectra(self,pixelR,pixelC):
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(range(self.wvlBinWidth),self.wvlMedians,'r-',label='median spectrum')
        ax.plot(np.arange(self.wvlBinWidth),self.spectra[pixelR,pixelC,:],label='pixel spectrum')
        plt.xlim([3750,self.wvlBinWidth])
        plt.legend(loc=2)
        plt.show()
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(range(self.wvlBinWidth),self.wvlMedians,'r-',label='median spectrum')
        ax.plot(np.arange(self.wvlBinWidth),self.pixelCalSpectra,label='cal pixel spectrum')
        plt.xlim([3750,self.wvlBinWidth])
        plt.legend(loc=2)
        plt.show()

    def plotPixelFactors(self,pixelR,pixelC):
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(self.wvlBinEdges[0:self.nWvlBins],self.wvlMedians,'r-',label='median spectrum')
        ax.plot(self.wvlBinEdges[0:self.nWvlBins],self.spectra[pixelR,pixelC,:],label='pixel spectrum')
        ax.plot(self.wvlBinEdges[0:self.nWvlBins],self.flatFactors[pixelR,pixelC,:],label='pixel weights')
        plt.legend(loc=1)
        plt.show()
        
    def displaySec(self,displaySec=0,pixelR=30,pixelC=30):
        secImg = np.zeros((self.nRow,self.nCol))
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                secImg[iRow,iCol] = self.obsFile.getPixelRawCount(iRow,iCol,displaySec,integrationTime=1)
        plt.matshow(secImg,vmax=np.mean(secImg)+2*np.std(secImg))
        plt.colorbar()
        plt.show()

    def displayCalibratedSec(self,displaySec=0,pixelR=30,pixelC=30):
        pulseMask = int(12*'1',2) #bitmask of 12 ones
        nBitsAfterEnergy = 44
        nBitsAfterBaseline = 20
        secImg = np.zeros((self.nRow,self.nCol),dtype='float')
        self.pixelCalSpectra = np.zeros(self.wvlBinWidth)
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                pixel = self.obsBeamImage[iRow][iCol]
                secs = self.obsFile.getNode('/'+pixel)
                for sec in secs:
                    for packet in sec:
                        packet = int(packet) #64 bit photon packet
                        #extract parabolaFitPeak as energy
                        parabolaFitPeak = (packet >> nBitsAfterEnergy) & pulseMask
                        baseline = (packet >> nBitsAfterBaseline) & pulseMask
                        pulseAmplitude = parabolaFitPeak-baseline
                        if pulseAmplitude > 0:
                            wavelength = self.applyWvlCal(row,col,pulseAmplitude)
                            flatFactor = self.flatFactors[iRow,iCol,wavelength]
                            secImg[iRow,iCol]+=flatFactor
                            if (iRow == pixelR and iCol == pixelC):
                                self.pixelCalSpectra[wavelength]+=flatFactor
                        
        print secImg[pixelR,pixelC]
        plt.matshow(secImg,vmax=np.mean(secImg)+2*np.std(secImg))
        plt.colorbar()
        plt.show()


def main():
    np.set_printoptions(threshold=np.nan)
    cal = FlatCal('obs_20120919-131142.h5','calsol_20120917-072537.h5','flatsol_20120919-131142.h5')
    #cal = FlatCal()
    #cal.loadObsFile('obs_20120919-131346.h5')
    #cal.loadObsFile('obs_20120919-131142.h5')
    #cal.loadSpectraFile()
    ob = ObsFile('obs_20120919-131346.h5')
    ob.loadWvlCalFile('calsol_20120917-072537.h5')
    ob.loadFlatCalFile('flatsol_20120919-131142.h5')
    cal.plotPixelFactors(17,12)
    ob.displaySec(17,12,integrationTime=-1)
    ob.displaySec(17,12,weighted=True,integrationTime=-1)
    ob.plotPixelSpectra(17,12)
    ob.plotPixelSpectra(17,12,weighted=True)
    #cal.displayCalibratedSec(pixelR=12,pixelC=35)
    #cal.plotPixelSpectra(12,35)


if __name__ == '__main__':
    main()



