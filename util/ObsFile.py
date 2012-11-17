#!/bin/python

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
        Creates an object to parse given file
        """
        self.loadFile(fileName)
        self.wvlCalFile = None
        self.flatCalFile = None

    def __del__(self):
        """
        Closes the reference file
        """
        self.file.close()
        try:
            self.wvlCalFile.close()
        except:
            pass

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

        #get the beam image.
        try:
            self.beamImage = self.file.getNode('/beammap/beamimage')
        except:
            print 'Can\'t access beamimage'
            sys.exit(2)

        beamShape = self.beamImage.shape
        self.nRow = beamShape[0]
        self.nCol = beamShape[1]

    def __iter__(self):
        """
        Allows easy iteration over pixels in obs file
        use with 'for pixel in ObsFileObject:'
        yields a single pixel dataset
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

    def getPixelWvlList(self,iRow,iCol,firstSec=0,integrationTime=-1,weighted = False):
        packetList = self.getPixelPacketList(iRow,iCol,firstSec,integrationTime)
        timestamps,parabolaPeaks,baselines = self.parsePhotonPackets(packetList)
        pulseHeights = np.array(parabolaPeaks,dtype='double') - np.array(baselines,dtype='double')
        wavelengths = self.convertToWvl(pulseHeights,iRow,iCol)
        if weighted == False:
            return wavelengths
        else:
            if len(wavelengths) > 0:
                binIndices = np.digitize(wavelengths,self.flatCalWvlBins[0:-1])
                #Mark photons with wavelengths below first bin and last bin with index 0, which becomes index -1 (below first bin already at 0)
                binIndices[binIndices >= self.nFlatCalWvlBins] = 0
                binIndices -= 1 #adjust so that indices match the lower edge of the wavelength bins
                pixelFlatWeights = self.flatWeights[iRow,iCol,binIndices]
                pixelFlatWeights[binIndices == -1] = 0.0
            else:
                pixelFlatWeights = np.array([])
            return wavelengths,pixelFlatWeights
            

    def getPixelCount(self,iRow,iCol,firstSec=0,integrationTime=-1,weighted=False):
        packetList = self.getPixelPacketList(iRow,iCol,firstSec,integrationTime)
        if weighted == False:
            return len(packetList)
        else:
            wavelengths,weights = self.getPixelWvlList(iRow,iCol,firstSec,integrationTime,weighted=True)
            return sum(weights)


    def getPixelPacketList(self,iRow,iCol,firstSec=0,integrationTime=-1):
        pixelData = self.getPixel(iRow,iCol)
        lastSec = firstSec+integrationTime
        if integrationTime == -1:
            lastSec = len(pixelData)
        pixelData = pixelData[firstSec:lastSec]
        packetList = np.concatenate(pixelData)
        return packetList

    def displaySec(self,pixelRow,pixelCol,firstSec=0,integrationTime=1,weighted=False):
        secImg = np.zeros((self.nRow,self.nCol))
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                secImg[iRow,iCol] = self.getPixelCount(iRow,iCol,firstSec,integrationTime=integrationTime,weighted=weighted)
        plt.matshow(secImg,vmax=np.mean(secImg)+2*np.std(secImg))
        plt.colorbar()
        plt.show()

    def getPixelSpectrum(self,pixelRow,pixelCol,firstSec=0,integrationTime=-1,weighted=False):
        if weighted == False:
            wvlList = self.getPixelWvlList(pixelRow,pixelCol,firstSec,integrationTime)
            spectrum,binEdges = np.histogram(wvlList,bins=self.nFlatCalWvlBins,range=(self.flatCalWvlBins[0],self.flatCalWvlBins[-1]))
        else:
            wvlList,weights = self.getPixelWvlList(pixelRow,pixelCol,firstSec,integrationTime,weighted=True)
            spectrum,binEdges = np.histogram(wvlList,bins=self.nFlatCalWvlBins,range=(self.flatCalWvlBins[0],self.flatCalWvlBins[-1]),weights=weights)
        return spectrum

    def plotPixelSpectra(self,pixelRow,pixelCol,firstSec=0,integrationTime=-1,weighted=False):
        spectrum = self.getPixelSpectrum(pixelRow,pixelCol,firstSec,integrationTime,weighted=weighted)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.flatCalWvlBins[0:-1],spectrum,label='spectrum for pixel[%d][%d]'%(pixelRow,pixelCol))
        plt.show()

    def convertToWvl(self,pulseHeights,iRow,iCol):
        xOffset = self.wvlCalTable[iRow,iCol,0]
        yOffset = self.wvlCalTable[iRow,iCol,1]
        amplitude = self.wvlCalTable[iRow,iCol,2]
        energies = amplitude*(pulseHeights-xOffset)**2+yOffset
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


#    def createPhotonList(self,photonListFileName,wvlCalFileName):
#        h = 4.135668e-15 #eV s
#        c = 2.998e8 #m/s
#        angstromPerMeter = 1e10
#        tickDuration = 1e-6 # s 
#
#        wvlCalTable = self.loadWvlCalFile(wvlCalFileName)
#        plTable = self.createPhotonListFile(photonListFileName)
#
#        for iRow in xrange(self.nRow):
#            for iCol in xrange(self.nCol):
#                for iSec,sec in enumerate(self.getPixel(iRow,iCol)):
#                    for photonPacket in sec:
#                        timestamp,parabolaPeak,baseline = self.parsePhotonPacket(photonPacket)
#                        pulseHeight = parabolaPeak - baseline
#                        xOffset = wvlCalTable[iRow,iCol,0]
#                        amplitude = wvlCalTable[iRow,iCol,2]
#                        yOffset = wvlCalTable[iRow,iCol,1]
#                        energy = amplitude*(pulseHeight-xOffset)**2+yOffset
#                        wavelength = h*c*angstromPerMeter/energy
#                        print wavelength 
#                        if wavelength > 0 and wavelength != np.inf:
#                            newRow = plTable.row
#                            newRow['Xpix'] = iCol
#                            newRow['Ypix'] = iRow
#                            newRow['ArrivalTime'] = iSec+timestamp*tickDuration
#                            newRow['Wavelength'] = wavelength
#                            newRow.append()
#        plTable.flush()

    def createPhotonListFile(self,photonListFileName):
        scratchDir = os.getenv('INTERM_PATH','/')
        fullPhotonListFileName = os.path.join(scratchDir,photonListFileName)
        if (os.path.exists(fullPhotonListFileName)):
            if confirm('Photon list file  %s exists. Overwrite?'%fullPhotonListFileName,resp=False) == False:
                exit(0)
        zlibFilter = tables.Filters(complevel=1,complib='zlib',fletcher32=False)
        plFile = tables.openFile(fullPhotonListFileName,mode='w')
        plGroup = plFile.createGroup('/','photons','Group containing photon list')
        plTable = plFile.createTable(plGroup,'photons',ArconsHeaders.PhotonList,'Photon List Data',filters=zlibFilter)
        return plTable

    def loadWvlCalFile(self,wvlCalFileName):
        scratchDir = '/ScienceData'
        wvlDir = os.path.join(scratchDir,'waveCalSolnFiles')
        fullWvlCalFileName = os.path.join(wvlDir,wvlCalFileName)
        if (not os.path.exists(fullWvlCalFileName)):
            print 'wavelength cal file does not exist: ',fullWvlCalFileName
            return
        self.wvlCalFile = tables.openFile(fullWvlCalFileName,mode='r')
        wvlCalData = self.wvlCalFile.root.wavecal.calsoln
        self.wvlCalTable = np.zeros([self.nRow,self.nCol,ObsFile.nCalCoeffs])
        for calPixel in wvlCalData:
            if calPixel['wave_flag'] == 0:
                self.wvlCalTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['polyfit']
        print 'Loaded Wave Cal file'
        return self.wvlCalTable

        
    def loadFlatCalFile(self,flatCalFileName):
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
