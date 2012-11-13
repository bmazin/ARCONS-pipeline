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

class FlatCal:
    def __init__(self,flatFileName=None):
        self.nEnergyBins = 1000000
        if flatFileName == None:
            self.loadFactors()
        else:
            self.loadFlatFile(flatFileName)
            self.loadFlatSpectra()
            self.calculateFactors()

    def __del__(self):
        try:
            self.flatFile.close()
            self.calFile.close()
        except AttributeError:#flatFile was never defined
            pass

    def loadFlatFile(self,flatFileName):
        self.flatFileName = flatFileName
        #make the full file name by joining the input name to the MKID_DATA_DIR (or . if the environment variable is not defined)
        dataDir = os.getenv('MKID_DATA_DIR','.')
        self.flatFullFileName = os.path.join(dataDir,self.flatFileName)
        print 'full file name is ',self.flatFullFileName
        if (not os.path.exists(self.flatFullFileName)):
            print 'file does not exist: ',self.flatFullFileName
            sys.exit(1)

        #open the hdf5 file
        self.flatFile = tables.openFile(self.flatFullFileName,mode='r')

        #get the header
        header = self.flatFile.root.header.header
        titles = header.colnames
        info = header[0] #header is a table with one row

        #get the beam image.
        try:
            self.beamImage = self.flatFile.getNode('/beammap/beamimage')
        except:
            print 'Can\'t access beamimage'
            sys.exit(2)

        #make arrays for each pixel's spectrum
        beamShape = self.beamImage.shape
        self.nRow = beamShape[0]
        self.nCol = beamShape[1]

    def loadObsFile(self,obsFileName):
        self.obsFileName = obsFileName
        #make the full file name by joining the input name to the MKID_DATA_DIR (or . if the environment variable is not defined)
        dataDir = os.getenv('MKID_DATA_DIR','.')
        self.obsFullFileName = os.path.join(dataDir,self.obsFileName)
        print 'full file name is ',self.obsFullFileName
        if (not os.path.exists(self.obsFullFileName)):
            print 'file does not exist: ',self.obsFullFileName
            sys.exit(1)

        #open the hdf5 file
        self.obsFile = tables.openFile(self.obsFullFileName,mode='r')

        #get the header
        header = self.obsFile.root.header.header
        titles = header.colnames
        info = header[0] #header is a table with one row

        #get the beam image.
        try:
            self.obsBeamImage = self.obsFile.getNode('/beammap/beamimage')
        except:
            print 'Can\'t access beamimage'
            sys.exit(2)

        #make arrays for each pixel's spectrum
        beamShape = self.obsBeamImage.shape
        self.nObsRow = beamShape[0]
        self.nObsCol = beamShape[1]
        
    def loadFlatSpectra(self):
        self.spectra = np.zeros((self.nRow,self.nCol,self.nEnergyBins),dtype=np.uint32)
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                print iRow,iCol,
                pixel = self.beamImage[iRow][iCol]
                secs = self.flatFile.getNode('/'+pixel)
                hgEnergy = self.loadSingleSpectra(secs,iRow,iCol)
                self.spectra[iRow,iCol,:]=hgEnergy
                print sum(hgEnergy)

        np.save('/ScienceData/intermediate/flatSpectra.npy',self.spectra)
        self.calculateMedians()
        return self.spectra

#    def applyWvlCal(self,row,col,pulseHeights):
#        wavelengths = np.array(self.calTable[row][col][0]*pulseHeights + self.calTable[row][col][1],dtype=np.uint32)
#        return wavelengths

    def loadCalFile(self,calFileName):
        intermediateDir = os.getenv('INTERM_PATH','.')
        self.calFullFileName = os.path.join(intermediateDir,calFileName)
        self.calFile = tables.openFile(self.calFullFileName,'r')
        self.nCalCoeffs = 2
        self.calTable = np.zeros([self.nRow,self.nCol,self.nCalCoeffs])
        calsoln = self.calFile.root.wavecal.calsoln
        for row in calsoln:
            print row
            if row['wave_flag'] == 0:
                self.calTable[row['pixelx']][row['pixely']] = row['polyfit']
        
    def loadSpectraFile(self):
        self.spectra = np.load('/ScienceData/intermediate/flatSpectra.npy')


    def calculateMedians(self):
        spectraByEnergy = np.reshape(self.spectra,[self.nRow*self.nCol,self.nEnergyBins])
        self.energyMedians = np.zeros(self.nEnergyBins)
        for iEnergy in xrange(self.nEnergyBins):
            print iEnergy,
            spectrum = spectraByEnergy[:,iEnergy]
            goodSpectrum = spectrum[spectrum != 0]
            self.energyMedians[iEnergy] = np.median(goodSpectrum)
            print self.energyMedians[iEnergy]
            
        return self.energyMedians
        
    def calculateFactors(self):
        fudgedSpectra = self.spectra
        fudgedSpectra[fudgedSpectra == 0] =1.0
        self.flatFactors = np.divide(self.energyMedians,self.spectra)
        self.flatFactors[self.flatFactors == np.inf]=0
        np.save('/ScienceData/intermediate/factors.npy',self.flatFactors)

    def loadFactors(self):
        self.flatFactors = np.load('/ScienceData/intermediate/factors.npy')
        #plt.matshow(self.flatFactors[:,:,4000])
        #plt.colorbar()
        #plt.title('weights at pulseh=4000')
        #plt.show()

        
    def plotMedian(self):
        print np.shape(self.energyMedians)
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(range(self.nEnergyBins),self.energyMedians,label='median spectrum')
        plt.xlim([3750,self.nEnergyBins])
        plt.legend(loc=2)
        plt.show()

    def plotPixelSpectra(self,pixelR,pixelC):
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(range(self.nEnergyBins),self.energyMedians,'r-',label='median spectrum')
        ax.plot(np.arange(self.nEnergyBins),self.spectra[pixelR,pixelC,:],label='pixel spectrum')
        plt.xlim([3750,self.nEnergyBins])
        plt.legend(loc=2)
        plt.show()
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(range(self.nEnergyBins),self.energyMedians,'r-',label='median spectrum')
        ax.plot(np.arange(self.nEnergyBins),self.pixelCalSpectra,label='cal pixel spectrum')
        plt.xlim([3750,self.nEnergyBins])
        plt.legend(loc=2)
        plt.show()

    def plotPixelFactors(self,pixelR,pixelC):
        fig = plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(range(self.nEnergyBins),self.energyMedians,'r-',label='median spectrum')
        ax.plot(np.arange(self.nEnergyBins),self.flatFactors[pixelR,pixelC,:],label='pixel weights')
        plt.xlim([3750,self.nEnergyBins])
        plt.legend(loc=2)
        plt.show()
        
    def loadSingleSpectra(self,secs,row,col):
        hgEnergy = np.zeros(self.nEnergyBins,dtype=np.uint32)
        pulseMask = int(12*'1',2) #bitmask of 12 ones
        nBitsAfterEnergy = 44
        nBitsAfterBaseline = 20
        
        for sec in secs:
            for packet in sec:
               packet = int(packet) #64 bit photon packet
               #extract parabolaFitPeak as energy
               parabolaFitPeak = (packet >> nBitsAfterEnergy) & pulseMask
               baseline = (packet >> nBitsAfterBaseline) & pulseMask
               pulseAmplitude = parabolaFitPeak-baseline
               if pulseAmplitude > 0:
                   wavelength = self.applyWvlCal(row,col,pulseAmplitude)
                   hgEnergy[wavelength] += 1
        return hgEnergy

    def loadWavelengthCal(self,wvlCalFile=None):
        #table converting pulse height (parabolat_fit-baseline) to wavelength in angstroms
        self.wvlTable = np.linspace(400,1100,num=self.nEnergyBins)

    def calibrateFile(self,rawObsFileName):
        obsFile = tables.openFile(rawObsFileName,mode='r')
        self.createPhotonList()

    def createPhotonList():
        photonListFileName = '/ScienceData/intermediate/pl.h5'
        photonListFile = tables.openFile(photonListFileName,mode='w',title='Photon List from %s'%rawObsFileName)

    def displaySec(self,displaySec=0,pixelR=30,pixelC=30):
        secImg = np.zeros((self.nObsRow,self.nObsCol))
        for iRow in xrange(self.nObsRow):
            for iCol in xrange(self.nObsCol):
                pixel = self.obsBeamImage[iRow][iCol]
                secs = self.obsFile.getNode('/'+pixel)
                for sec in secs:
                    secImg[iRow,iCol]+=len(sec)
        print secImg[pixelR,pixelC]
        plt.matshow(secImg,vmax=np.mean(secImg)+2*np.std(secImg))
        plt.colorbar()
        plt.show()

    def displayCalibratedSec(self,displaySec=0,pixelR=30,pixelC=30):
        pulseMask = int(12*'1',2) #bitmask of 12 ones
        nBitsAfterEnergy = 44
        nBitsAfterBaseline = 20
        secImg = np.zeros((self.nObsRow,self.nObsCol),dtype='float')
        self.pixelCalSpectra = np.zeros(self.nEnergyBins)
        for iRow in xrange(self.nObsRow):
            for iCol in xrange(self.nObsCol):
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
    #cal = FlatCal('obs_20120919-131142.h5')
    cal = FlatCal()
    cal.loadFlatFile('obs_20120919-131142.h5')
    cal.loadCalFile('calsol_20120917-072537.h5')
    cal.loadObsFile('obs_20120919-131346.h5')
    cal.loadFlatSpectra()
    #cal.loadSpectraFile()
    cal.calculateMedians()
    cal.calculateFactors()
    cal.plotPixelFactors(17,12)
    cal.displaySec()
    cal.displayCalibratedSec(pixelR=17,pixelC=12)
    cal.plotPixelSpectra(17,12)
    
    cal.plotPixelFactors(12,35)
    cal.displaySec()
    cal.displayCalibratedSec(pixelR=12,pixelC=35)
    cal.plotPixelSpectra(12,35)


if __name__ == '__main__':
    main()



