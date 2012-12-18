#!/bin/python
'''
Author: Matt Strader        Date: August 19, 2012

The class ObsFile is an interface to observation files.  It provides methods for typical ways of accessing and viewing observation data.  It can also load and apply wavelength and flat calibration.  With calibrations loaded, it can write the obs file out as a photon list

Looks for observation files in $MKID_DATA_DIR and calibration files organized in $INTERM_PATH (intermediate or scratch path)
'''

import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.utils import confirm
from headers import ArconsHeaders
from util import utils

class ObsFile:
    h = 4.135668e-15 #eV s
    c = 2.998e8 #m/s
    angstromPerMeter = 1e10
    tickDuration = 1e-6 #s
    nCalCoeffs = 3
    def __init__(self,fileName):
        """
        load the given file with fileName relative to $MKID_DATA_DIR
        """
        self.loadFile(fileName)
        self.wvlCalFile = None #initialize to None for an easy test of whether a cal file has been loaded
        self.flatCalFile = None
        

    def __del__(self):
        """
        Closes the obs file and any cal files that are open
        """
        self.file.close()
        try:
            self.wvlCalFile.close()
            self.flatCalFile.close()
        except:
            pass

    def getFromHeader(self,name):
        return self.info[self.titles.index(name)]


    def loadFile(self,fileName):
        """
        Opens file and loads obs file attributes and beammap
        """
        self.fileName = fileName
        #make the full file name by joining the input name to the MKID_DATA_DIR (or . if the environment variable is not defined)
        dataDir = os.getenv('MKID_DATA_DIR','/')
        self.fullFileName = os.path.join(dataDir,self.fileName)
        if (not os.path.exists(self.fullFileName)):
            print 'file does not exist: ',self.fullFileName
            sys.exit(1)

        #open the hdf5 file
        self.file = tables.openFile(self.fullFileName,mode='r')

        #get the header
        header = self.file.root.header.header
        titles = header.colnames
        info = header[0] #header is a table with one row
        self.titles = titles
        self.info = info

        #get the beam image.
        try:
            self.beamImage = self.file.getNode('/beammap/beamimage').read()
        except:
            print 'Can\'t access beamimage'
            sys.exit(2)

        beamShape = self.beamImage.shape
        self.nRow = beamShape[0]
        self.nCol = beamShape[1]

    def __iter__(self):
        """
        Allows easy iteration over pixels in obs file
        use with 'for pixel in obsFileObject:'
        yields a single pixel h5 dataset
        """
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                pixelLabel = self.beamImage[iRow][iCol]
                pixelData = self.file.getNode('/'+pixelLabel)
                yield pixelData

    def getPixel(self,iRow,iCol):
        """
        retrieves a pixel using the file's attached beammap
        """
        pixelLabel = self.beamImage[iRow][iCol]
        pixelData = self.file.getNode('/'+pixelLabel).read()
        return pixelData

    def getPixelWvlList(self,iRow,iCol,firstSec=0,integrationTime=-1,weighted = False,getTimes=False):
        """
        returns a numpy array of photon wavelengths for a given pixel, integrated from firstSec to firstSec+integrationTime.
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        if getTimes is True, returns timestamps,wavelengths
        """
        if getTimes == False:
            packetList = self.getPixelPacketList(iRow,iCol,firstSec,integrationTime)
            timestamps,parabolaPeaks,baselines = self.parsePhotonPackets(packetList)
        else:
            timestamps,parabolaPeaks,baselines = self.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
            
        pulseHeights = np.array(parabolaPeaks,dtype='double') - np.array(baselines,dtype='double')
        wavelengths = self.convertToWvl(pulseHeights,iRow,iCol)
        if getTimes == False:
            return wavelengths
        else:
            return timestamps,wavelengths
            

    def getPixelCount(self,iRow,iCol,firstSec=0,integrationTime=-1,weighted=False):
        """
        returns the number of photons received in a given pixel from firstSec to firstSec + integrationTime
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
        packetList = self.getPixelPacketList(iRow,iCol,firstSec,integrationTime)
        if weighted == False:
            return len(packetList)
        else:
            weightedSpectrum,binEdges = self.getPixelSpectrum(iRow,iCol,firstSec,integrationTime,weighted=True)
            return sum(weightedSpectrum)


    def getPixelPacketList(self,iRow,iCol,firstSec=0,integrationTime=-1):
        """
        returns a numpy array of 64-bit photon packets for a given pixel, integrated from firstSec to firstSec+integrationTime.
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
        pixelData = self.getPixel(iRow,iCol)
        lastSec = firstSec+integrationTime
        if integrationTime == -1:
            lastSec = len(pixelData)
        pixelData = pixelData[firstSec:lastSec]
        packetList = np.concatenate(pixelData)
        return packetList

    def getTimedPacketList(self,iRow,iCol,firstSec=0,integrationTime=-1):
        pixelData = self.getPixel(iRow,iCol)
        lastSec = firstSec+integrationTime
        if integrationTime == -1:
            lastSec = len(pixelData)
        pixelData = pixelData[firstSec:lastSec]

        timestamps = []
        baselines = []
        peakHeights = []

        for t in range(len(pixelData)):
            times,peaks,bases = self.parsePhotonPackets(pixelData[t])
            times =1e-6*times+t+firstSec
            timestamps.append(times)
            baselines.append(bases)
            peakHeights.append(peaks)
            
        timestamps = np.concatenate(timestamps)
        baselines = np.concatenate(baselines)
        peakHeights = np.concatenate(peakHeights)
        return timestamps,peakHeights,baselines


        

    def displaySec(self,firstSec=0,integrationTime=1,weighted=False):
        """
        plots a time-flattened image of the counts integrated from firstSec to firstSec+integrationTime
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
        secImg = np.zeros((self.nRow,self.nCol))
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                secImg[iRow,iCol] = self.getPixelCount(iRow,iCol,firstSec,integrationTime=integrationTime,weighted=weighted)
        plt.matshow(secImg,vmax=np.mean(secImg)+2*np.std(secImg))
        plt.colorbar()
        plt.show()

    def getPixelSpectrum(self,pixelRow,pixelCol,firstSec=0,integrationTime=-1,weighted=False,wvlStart=3000,wvlStop=13000,wvlBinWidth=100):
        """
        returns a spectral histogram of a given pixel integrated from firstSec to firstSec+integrationTime, and an array giving the cutoff wavelengths used to bin the wavelength values
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied and flat cal wavelength bins are used
        if weighted is False, wavelength bin parameters given are used
        """
        wvlList = self.getPixelWvlList(pixelRow,pixelCol,firstSec,integrationTime)
        nWvlBins = int((wvlStop - wvlStart)/wvlBinWidth)
        spectrum,wvlBinEdges = np.histogram(wvlList,bins=nWvlBins,range=(wvlStart,wvlStop))
        if weighted == False:
            spectrum,wvlBinEdges = np.histogram(wvlList,bins=nWvlBins,range=(wvlStart,wvlStop))
        else:
            spectrum,wvlBinEdges = np.histogram(wvlList,bins=self.flatCalWvlBins)
            spectrum = spectrum * self.flatWeights[pixelRow,pixelCol]
        return spectrum,wvlBinEdges

    def plotPixelSpectra(self,pixelRow,pixelCol,firstSec=0,integrationTime=-1,weighted=False):
        """
        plots the wavelength calibrated spectrum of a given pixel integrated over a given time
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
        spectrum,binEdges = self.getPixelSpectrum(pixelRow,pixelCol,firstSec,integrationTime,weighted=weighted)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.flatCalWvlBins[0:-1],spectrum,label='spectrum for pixel[%d][%d]'%(pixelRow,pixelCol))
        plt.show()

    def convertToWvl(self,pulseHeights,iRow,iCol,excludeBad = True):
        """
        applies wavelength calibration to a list of photon pulse heights
        if excludeBad is True, wavelengths calculated as np.inf are excised from the array returned
        """

        xOffset = self.wvlCalTable[iRow,iCol,0]
        yOffset = self.wvlCalTable[iRow,iCol,1]
        amplitude = self.wvlCalTable[iRow,iCol,2]
        energies = amplitude*(pulseHeights-xOffset)**2+yOffset
        if excludeBad == True:
            energies = energies[energies != 0]
        wavelengths = ObsFile.h*ObsFile.c*ObsFile.angstromPerMeter/energies
        return wavelengths


    def parsePhotonPackets(self,packets):
        """
        Parses an array of uint64 packets with the obs file format
        8 bits - channel
        12 bits - Parabola Fit Peak Height
        12 bits - Sampled Peak Height
        12 bits - Low pass filter baseline
        20 bits - Microsecond timestamp
        """
        nBitsAfterParabolaPeak = 44
        nBitsAfterBaseline = 20
        nBitsInPulseHeight = 12
        nBitsInTimestamp = 20
        pulseMask = int(nBitsInPulseHeight*'1',2) #bitmask of 12 ones
        timestampMask = int(nBitsInTimestamp*'1',2) #bitmask of 12 ones
        packets = np.array(packets,dtype='uint64') #64 bit photon packet
        timestamps = np.bitwise_and(packets,timestampMask)
        parabolaFitPeaks = np.bitwise_and(np.right_shift(packets,nBitsAfterParabolaPeak),pulseMask)
        baselines = np.bitwise_and(np.right_shift(packets,nBitsAfterBaseline),pulseMask)
        return timestamps,parabolaFitPeaks,baselines


    def createPhotonList(self,photonListFileName):
        """
        writes out the photon list for this obs file at $INTERM_PATH/photonListFileName
        *Should be changed to use same filename as obs file but with pl prefix
        """
        tickDuration = 1e-6 # s 
        plTable = self.createEmptyPhotonListFile(photonListFileName)

        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                flag = self.wvlFlagTable[iRow,iCol]
                if flag == 0:#only write photons in good pixels
                    print iRow,iCol
                    wvlError = self.wvlErrorTable[iRow,iCol]
                    flatWeights = self.flatWeights[iRow,iCol]
                    #go through the list of seconds in a pixel dataset
                    for iSec,secData in enumerate(self.getPixel(iRow,iCol)):
                        timestamps,parabolaPeaks,baselines = self.parsePhotonPackets(secData)
                        pulseHeights = np.array(parabolaPeaks,dtype='double') - np.array(baselines,dtype='double')
                        timestamps = iSec + tickDuration*timestamps
                        wavelengths = self.convertToWvl(pulseHeights,iRow,iCol,excludeBad=False)
                        #calculate what wavelength bins each photon falls into to see which flat cal factor should be applied
                        if len(wavelengths) > 0:
                            binIndices = np.digitize(wavelengths,self.flatCalWvlBins[0:-1])
                        else:
                            binIndices = np.array([])
                        for iPhoton in xrange(len(timestamps)):
                            if wavelengths[iPhoton] > 0 and wavelengths[iPhoton] != np.inf and binIndices[iPhoton] < len(flatWeights):
                                #create a new row for the photon list
                                newRow = plTable.row
                                newRow['Xpix'] = iCol
                                newRow['Ypix'] = iRow
                                newRow['ArrivalTime'] = timestamps[iPhoton]
                                newRow['Wavelength'] = wavelengths[iPhoton]
                                newRow['WaveError'] = wvlError
                                newRow['Flag'] = flag
                                newRow['FlatWeight'] = flatWeights[binIndices[iPhoton]]
                                newRow.append()
        plTable.flush()

    def createEmptyPhotonListFile(self,photonListFileName):
        """
        creates a file at $INTERM_PATH/photonLists/photonListFileName
        using header in headers.ArconsHeaders
        """
        scratchDir = os.getenv('INTERM_PATH','/')
        plDir = os.path.join(scratchDir,'photonLists')
        fullPhotonListFileName = os.path.join(scratchDir,photonListFileName)
        if (os.path.exists(fullPhotonListFileName)):
            if confirm('Photon list file  %s exists. Overwrite?'%fullPhotonListFileName,defaultResponse=False) == False:
                exit(0)
        zlibFilter = tables.Filters(complevel=1,complib='zlib',fletcher32=False)
        plFile = tables.openFile(fullPhotonListFileName,mode='w')
        plGroup = plFile.createGroup('/','photons','Group containing photon list')
        plTable = plFile.createTable(plGroup,'photons',ArconsHeaders.PhotonList,'Photon List Data',filters=zlibFilter)
        return plTable



    def loadWvlCalFile(self,wvlCalFileName):
        """
        loads the wavelength cal coefficients from a given file
        """
        scratchDir = '/ScienceData'
        wvlDir = os.path.join(scratchDir,'waveCalSolnFiles')
        fullWvlCalFileName = os.path.join(wvlDir,wvlCalFileName)
        if (not os.path.exists(fullWvlCalFileName)):
            print 'wavelength cal file does not exist: ',fullWvlCalFileName
            return
        self.wvlCalFile = tables.openFile(fullWvlCalFileName,mode='r')
        wvlCalData = self.wvlCalFile.root.wavecal.calsoln
        self.wvlCalTable = np.zeros([self.nRow,self.nCol,ObsFile.nCalCoeffs])
        self.wvlErrorTable = np.zeros([self.nRow,self.nCol])
        self.wvlFlagTable = np.zeros([self.nRow,self.nCol])
        for calPixel in wvlCalData:
            self.wvlFlagTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['wave_flag']
            self.wvlErrorTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['sigma']
            if calPixel['wave_flag'] == 0:
                self.wvlCalTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['polyfit']

        
    def loadFlatCalFile(self,flatCalFileName):
        """
        loads the flat cal factors from the given file
        """
        scratchDir = os.getenv('INTERM_PATH','/')
        flatCalPath = os.path.join(scratchDir,'flatCalSolnFiles')
        fullFlatCalFileName = os.path.join(flatCalPath,flatCalFileName)
        if (not os.path.exists(fullFlatCalFileName)):
            print 'flat cal file does not exist: ',fullFlatCalFileName
            return
        self.flatCalFile = tables.openFile(fullFlatCalFileName,mode='r')
        self.flatWeights = self.flatCalFile.root.flatcal.weights.read()
        self.flatCalWvlBins = self.flatCalFile.root.flatcal.wavelengthBins.read()
        self.nFlatCalWvlBins = self.flatWeights.shape[2]
