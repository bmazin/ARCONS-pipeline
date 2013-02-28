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
from headers import ArconsHeaders
from util import utils
from interval import interval, inf, imath
from util.FileName import FileName

class ObsFile:
    h = 4.135668e-15 #eV s
    c = 2.998e8 #m/s
    angstromPerMeter = 1e10
    nCalCoeffs = 3
    def __init__(self,fileName):
        """
        load the given file with fileName relative to $MKID_DATA_DIR
        """
        self.loadFile(fileName)
        self.wvlCalFile = None #initialize to None for an easy test of whether a cal file has been loaded
        self.flatCalFile = None
        self.fluxCalFile = None
        self.wvlLowerLimit = None
        self.wvlUpperLimit = None
        

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
        if (os.path.isabs(fileName)):
            self.fileName = os.path.basename(fileName)
            self.fullFileName = fileName
        else:
            self.fileName = fileName
            # make the full file name by joining the input name 
            # to the MKID_DATA_DIR (or . if the environment variable 
            # is not defined)
            dataDir = os.getenv('MKID_DATA_DIR','/')
            self.fullFileName = os.path.join(dataDir,self.fileName)

        if (not os.path.exists(self.fullFileName)):
            print 'file does not exist: ',self.fullFileName
            sys.exit(1)

        #open the hdf5 file
        self.file = tables.openFile(self.fullFileName,mode='r')

        #get the header
        self.header = self.file.root.header.header
        self.titles = self.header.colnames
        self.info = self.header[0] #header is a table with one row

        # Useful information about data format set here.
        # For now, set all of these as constants.
        # If we get data taken with different parameters, straighten
        # that all out here.

        ## These parameters are for LICK2012 and PAL2012 data
        self.tickDuration = 1e-6 #s
        self.ticksPerSec = int(1.0/self.tickDuration)
        self.intervalAll = interval[0.0, (1.0/self.tickDuration)-1]
        self.nonAllocPixelName = '/r0/p250/'
        #  8 bits - channel
        # 12 bits - Parabola Fit Peak Height
        # 12 bits - Sampled Peak Height
        # 12 bits - Low pass filter baseline
        # 20 bits - Microsecond timestamp

        self.nBitsAfterParabolaPeak = 44
        self.nBitsAfterBaseline = 20
        self.nBitsInPulseHeight = 12
        self.nBitsInTimestamp = 20

        #bitmask of 12 ones
        self.pulseMask = int(self.nBitsInPulseHeight*'1',2) 
        #bitmask of 20 ones
        self.timestampMask = int(self.nBitsInTimestamp*'1',2) 

        #get the beam image.
        try:
            self.beamImage = self.file.getNode('/beammap/beamimage').read()
        except:
            print 'Can\'t access beamimage'
            sys.exit(2)

        beamShape = self.beamImage.shape
        self.nRow = beamShape[0]
        self.nCol = beamShape[1]


    def getFromHeader(self,name):
        """
        returns the value for name in the header
        """
        return self.info[self.titles.index(name)]

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

    #def getPixel(self,iRow,iCol):
    #    """
    #    retrieves a pixel using the file's attached beammap
    #    """
    #    pixelLabel = self.beamImage[iRow][iCol]
    #    pixelData = self.file.getNode('/'+pixelLabel).read()
    #    return pixelData

    def getPixel(self, iRow, iCol, firstSec=0, integrationTime= -1):
        """
        Retrieves a pixel using the file's attached beammap.
        If firstSec/integrationTime are provided, only data from the time 
        interval 'firstSec' to firstSec+integrationTime are returned.
        For now firstSec and integrationTime can only be integers.
        If integrationTime is -1, all data after firstSec are returned.
        """
        pixelLabel = self.beamImage[iRow][iCol]
        pixelNode = self.file.getNode('/' + pixelLabel)
        if integrationTime == -1:
            lastSec = pixelNode.nrows
        else:
            lastSec = firstSec + integrationTime
        pixelData = pixelNode.read(firstSec, lastSec)
        return pixelData


    def getPixelWvlList(self,iRow,iCol,firstSec=0,integrationTime=-1,getTimes=False):
        """
        returns a numpy array of photon wavelengths for a given pixel, integrated from firstSec to firstSec+integrationTime.
        if integrationTime is -1, All time after firstSec is used.  
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
            

    def getPixelCount(self,iRow,iCol,firstSec=0,integrationTime=-1,weighted=False,fluxWeighted=False):
        """
        returns the number of photons received in a given pixel from firstSec to firstSec + integrationTime
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
        packetList = self.getPixelPacketList(iRow,iCol,firstSec,integrationTime)
        if weighted == False and self.wvlLowerLimit == None and self.wvlUpperLimit == None:
            return len(packetList)
        else:
            if fluxWeighted==True:
                weightedSpectrum,binEdges = self.getPixelSpectrum(iRow,iCol,firstSec,integrationTime,weighted=True,fluxWeighted=True)
            else:
                weightedSpectrum,binEdges = self.getPixelSpectrum(iRow,iCol,firstSec,integrationTime,weighted=True,fluxWeighted=False)
            return sum(weightedSpectrum)


#    def getPixelPacketList(self,iRow,iCol,firstSec=0,integrationTime=-1):
#        """
#        returns a numpy array of 64-bit photon packets for a given pixel, integrated from firstSec to firstSec+integrationTime.
#        if integrationTime is -1, All time after firstSec is used.  
#        if weighted is True, flat cal weights are applied
#        """
#        pixelData = self.getPixel(iRow,iCol)
#        lastSec = firstSec+integrationTime
#        if integrationTime == -1:
#            lastSec = len(pixelData)
#        pixelData = pixelData[firstSec:lastSec]
#        packetList = np.concatenate(pixelData)
#        return packetList

    def getPixelPacketList(self, iRow, iCol, firstSec=0, integrationTime= -1):
        """
        returns a numpy array of 64-bit photon packets for a given pixel, integrated from firstSec to firstSec+integrationTime.
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
        pixelData = self.getPixel(iRow, iCol, firstSec, integrationTime)
        lastSec = firstSec + integrationTime
        if integrationTime == -1:
            lastSec = len(pixelData)
        packetList = np.concatenate(pixelData)
        return packetList



    def getTimedPacketList(self,iRow,iCol,firstSec=0,integrationTime=-1):
        """
        Parses an array of uint64 packets with the obs file format,and makes timestamps absolute
        inter is an interval of time values to mask out
        returns a list of timestamps,parabolaFitPeaks,baselines
        parses packets from firstSec to firstSec+integrationTime.
        if integrationTime is -1, All time after firstSec is used.  
        """
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
            times =firstSec+self.tickDuration*times+t
            timestamps.append(times)
            baselines.append(bases)
            peakHeights.append(peaks)
            
        timestamps = np.concatenate(timestamps)
        baselines = np.concatenate(baselines)
        peakHeights = np.concatenate(peakHeights)
        return timestamps,peakHeights,baselines


    def getPixelCountImage(self, firstSec=0, integrationTime=-1, weighted=False,fluxWeighted=False):
        """
        Return a time-flattened image of the counts integrated from firstSec to firstSec+integrationTime.
        If integration time is -1, all time after firstSec is used.
        If weighted is True, flat cal weights are applied. JvE 12/28/12
        If fluxWeighted is True, flux cal weights are applied. SM 2/7/13
        """
        secImg = np.zeros((self.nRow, self.nCol))
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                secImg[iRow, iCol] = self.getPixelCount(iRow, iCol, firstSec, integrationTime, weighted,fluxWeighted)
        return secImg
    

    def displaySec(self,firstSec=0,integrationTime=-1,weighted=False,fluxWeighted=False,plotTitle='',nSdevMax=2):
        """
        plots a time-flattened image of the counts integrated from firstSec to firstSec+integrationTime
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
#       secImg = np.zeros((self.nRow,self.nCol))
#       for iRow in xrange(self.nRow):
#           for iCol in xrange(self.nCol):
#               secImg[iRow,iCol] = self.getPixelCount(iRow,iCol,firstSec,integrationTime=integrationTime,weighted=weighted)
        secImg = self.getPixelCountImage(firstSec, integrationTime, weighted, fluxWeighted)
#        plt.matshow(secImg,vmax=np.mean(secImg)+2*np.std(secImg))
#        plt.colorbar()
#        plt.show()
        utils.plotArray(secImg,cbar=True,normMax=np.mean(secImg)+nSdevMax*np.std(secImg),plotTitle=plotTitle)
        
    def getPixelSpectrum(self,pixelRow,pixelCol,firstSec=0,integrationTime=-1,weighted=False,fluxWeighted=False,wvlStart=3000,wvlStop=13000,wvlBinWidth=None,energyBinWidth=None,wvlBinEdges=None):

        """
        returns a spectral histogram of a given pixel integrated from firstSec to firstSec+integrationTime, and an array giving the cutoff wavelengths used to bin the wavelength values
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        if weighted is False, flat cal weights are not applied
        the wavelength bins used depends on the parameters given.
        If energyBinWidth is specified, the wavelength bins use fixed energy bin widths
        If wvlBinWidth is specified, the wavelength bins use fixed wavelength bin widths
        If neither is specified and/or if weighted is True, the flat cal wvlBinEdges is used
        """
        wvlList = self.getPixelWvlList(pixelRow,pixelCol,firstSec,integrationTime)
        
        if self.flatCalFile != None and ((wvlBinEdges == None and energyBinWidth == None and wvlBinWidth == None) or weighted == True):
        #We've loaded a flat cal already, which has wvlBinEdges defined, and no other bin edges parameters are specified to override it.
            spectrum,wvlBinEdges = np.histogram(wvlList,bins=self.flatCalWvlBins)
            if weighted == True:#Need to apply flat weights by wavelenth
                spectrum = spectrum * self.flatWeights[pixelRow,pixelCol]
                if fluxWeighted == True:
                    spectrum = spectrum*self.fluxWeights
        else:
            if weighted == True:
                raise ValueError('when weighted=True, flatCal wvl bins are used, so wvlBinEdges,wvlBinWidth,energyBinWidth,wvlStart,wvlStop should not be specified')
            if wvlBinEdges == None:#We need to construct wvlBinEdges array
                if energyBinWidth != None:#Fixed energy binwidth specified
                    #Construct array with variable wvl binwidths
                    wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth=energyBinWidth,wvlStart=wvlStart,wvlStop=wvlStop)
                    spectrum,wvlBinEdges = np.histogram(wvlList,bins=wvlBinEdges)
                elif wvlBinWidth != None:#Fixed wvl binwidth specified
                    nWvlBins = int((wvlStop - wvlStart)/wvlBinWidth)
                    spectrum,wvlBinEdges = np.histogram(wvlList,bins=nWvlBins,range=(wvlStart,wvlStop))
                else:
                    raise ValueError('getPixelSpectrum needs either wvlBinWidth,wvlBinEnergy, or wvlBinEdges')
            else:#We are given wvlBinEdges array
                spectrum,wvlBinEdges = np.histogram(wvlList,bins=wvlBinEdges)
        return spectrum,wvlBinEdges
                
    def plotPixelSpectra(self,pixelRow,pixelCol,firstSec=0,integrationTime=-1,weighted=False,fluxWeighted=False):
        """
        plots the wavelength calibrated spectrum of a given pixel integrated over a given time
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
        spectrum,binEdges = self.getPixelSpectrum(pixelRow,pixelCol,firstSec,integrationTime,weighted=weighted, fluxWeighted=fluxWeighted)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.flatCalWvlBins[0:-1],spectrum,label='spectrum for pixel[%d][%d]'%(pixelRow,pixelCol))
        plt.show()

    def getApertureSpectrum(self,pixelRow,pixelCol,radius,weighted=True,fluxWeighted=False,hotPixMask=None):
	'''
	Creates a spectrum a group of pixels.  Aperture is defined by pixelRow and pixelCol of
	center, as well as radius.  Wave and flat cals should be loaded before using this
	function.  If no hot pixel mask is applied, taking the median of the sky rather than
	the average to account for high hot pixel counts.
	Will add more options as other pieces of pipeline become more refined.
	'''
	print 'Creating dead pixel mask...'
	deadMask = self.getDeadPixels()
	print 'Creating aperture mask...'
	apertureMask=utils.aperture(pixelCol,pixelRow,radius=radius)
	print 'Creating sky mask...'
	bigMask = utils.aperture(pixelCol,pixelRow,radius=radius*3)
	skyMask = bigMask-apertureMask
	if hotPixMask == None:
	    y_values,x_values= np.where(np.logical_and(apertureMask==0,deadMask==1))
	    y_sky,x_sky=np.where(np.logical_and(skyMask==0,deadMask==1))
	else:
	    y_values,x_values= np.where(np.logical_and(np.logical_and(apertureMask==0,deadMask==1),hotPixMask==0))
	    y_sky,x_sky=np.where(np.logical_and(np.logical_and(skyMask==0,deadMask==1),hotPixMask==0))
	wvlBinEdges = self.getPixelSpectrum(y_values[0],x_values[0],weighted=weighted)[1]
	print 'Creating average sky spectrum...'
	skyspectrum=[]
	for i in range(len(x_sky)):
	    skyspectrum.append(self.getPixelSpectrum(y_sky[i],x_sky[i],weighted=weighted,fluxWeighted=fluxWeighted)[0])
	sky_array = np.zeros(len(skyspectrum[0]))
	for j in range(len(skyspectrum[0])):
	    ispectrum = np.zeros(len(skyspectrum))
	    for i in range(len(skyspectrum)):    
	        ispectrum[i]= skyspectrum[i][j]
	    if hotPixMask==None:
		sky_array[j] = np.median(ispectrum)
	    else:
	        sky_array[j] = np.average(ispectrum)
	print 'Creating sky subtracted spectrum...'
	spectrum=[]
	for i in range(len(x_values)):
	    spectrum.append(self.getPixelSpectrum(y_values[i],x_values[i],weighted=weighted,fluxWeighted=fluxWeighted)[0]-sky_array)
	summed_array = np.zeros(len(spectrum[0]))
	for j in range(len(spectrum[0])):
	    ispectrum = np.zeros(len(spectrum))
	    for i in range(len(spectrum)):    
	        ispectrum[i]= spectrum[i][j]
	    summed_array[j] = np.sum(ispectrum)
	for i in range(len(summed_array)):
	    summed_array[i] /= (wvlBinEdges[i+1]-wvlBinEdges[i])
	return summed_array,wvlBinEdges

    def plotApertureSpectrum(self,pixelRow,pixelCol,radius,weighted=True,fluxWeighted=False,hotPixMask=None):
	summed_array,bin_edges=self.getApertureSpectrum(pixelCol=pixelCol,pixelRow=pixelRow,radius=radius,weighted=weighted,fluxWeighted=fluxWeighted,hotPixMask=hotPixMask)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(bin_edges[12:-2],summed_array[12:-1])
	plt.xlabel('Wavelength ($\AA$)')
	plt.ylabel('Counts')
	plt.show()
	
	

    def setWvlCutoffs(self,wvlLowerLimit=3000,wvlUpperLimit=8000):
        """
        Sets wavelength cutoffs so that if convertToWvl(excludeBad=True) is called
        wavelengths outside these limits are excluded.  To remove limits
        set wvlLowerLimit and/or wvlUpperLimit to None.  To use the wavecal limits
        for each individual pixel, set wvlLowerLimit and/or wvlUpperLimit to -1 
        """
        self.wvlLowerLimit = wvlLowerLimit
        self.wvlUpperLimit = wvlUpperLimit

    def convertToWvl(self,pulseHeights,iRow,iCol,excludeBad = True):
        """
        applies wavelength calibration to a list of photon pulse heights
        if excludeBad is True, wavelengths calculated as np.inf are excised from the array returned, as are wavelengths outside the fit limits of the wavecal
        """

        xOffset = self.wvlCalTable[iRow,iCol,0]
        yOffset = self.wvlCalTable[iRow,iCol,1]
        amplitude = self.wvlCalTable[iRow,iCol,2]
        wvlCalLowerLimit = self.wvlRangeTable[iRow,iCol,0]
        wvlCalUpperLimit = self.wvlRangeTable[iRow,iCol,1]
        energies = amplitude*(pulseHeights-xOffset)**2+yOffset
        if excludeBad == True:
            energies = energies[energies != 0]
        wavelengths = ObsFile.h*ObsFile.c*ObsFile.angstromPerMeter/energies
        if excludeBad == True and self.wvlLowerLimit == -1:
            wavelengths = wavelengths[wvlCalLowerLimit < wavelengths]
        elif excludeBad == True and self.wvlLowerLimit != None:
            wavelengths = wavelengths[self.wvlLowerLimit < wavelengths]
        if excludeBad == True and self.wvlUpperLimit == -1:
            wavelengths = wavelengths[wavelengths < wvlCalUpperLimit]
        elif excludeBad == True and self.wvlUpperLimit != None:
            wavelengths = wavelengths[wavelengths < self.wvlUpperLimit]
#            if len(wavelengths) > 0 and self.flatCalFile != None:
#                #filter out wavelengths without a valid flat weight
#                pixelFlags = self.flatFlags[iRow,iCol]
#                binIndices = np.digitize(wavelengths,self.flatCalWvlBins)-1
#                wavelengths=wavelengths[np.logical_and(binIndices>=0,binIndices<len(pixelFlags))]
#                binIndices=binIndices[np.logical_and(binIndices>=0,binIndices<len(pixelFlags))]
#                flags = pixelFlags[binIndices]
#                wavelengths = wavelengths[flags==1]
        return wavelengths


    def parsePhotonPackets(self,packets, inter=interval(), 
                           doParabolaFitPeaks=True, doBaselines=True):
        """
        Parses an array of uint64 packets with the obs file format
        inter is an interval of time values to mask out
        returns a list of timestamps,parabolaFitPeaks,baselines

        """

        # first special case:  inter masks out everything so return zero-length
        # numpy arrays
        if (inter == self.intervalAll):
            timestamps = np.arange(0)
            parabolaFitPeaks = np.arange(0)
            baselines = np.arange(0)
        else:
            # parse all packets
            packetsAll = np.array(packets,dtype='uint64') #64 bit photon packet
            timestampsAll = np.bitwise_and(packets,self.timestampMask)

            if doParabolaFitPeaks:
                parabolaFitPeaksAll = np.bitwise_and(\
                    np.right_shift(packets,self.nBitsAfterParabolaPeak),\
                        self.pulseMask)
            else:
                parabolaFitPeaksAll = np.arange(0)

            if doBaselines:
                baselinesAll = np.bitwise_and(\
                    np.right_shift(packets,self.nBitsAfterBaseline),\
                        self.pulseMask)
            else:
                baselinesAll = np.arange(0)

            if inter==interval() or len(timestampsAll) == 0:
                # nothing excluded or nothing to exclude
                # so return all unpacked values
                timestamps = timestampsAll
                parabolaFitPeaks = parabolaFitPeaksAll
                baselines = baselinesAll
            else:
                # there is a non-trivial set of times to mask. 
                slices = calculateSlices(inter, timestampsAll)
                timestamps = repackArray(timestampsAll, slices)
                parabolaFitPeaks = repackArray(parabolaFitPeaksAll, slices)
                baselines = repackArray(baselinesAll, slices)
        # return the values filled in above
        return timestamps,parabolaFitPeaks,baselines

    def writePhotonList(self):
        """
        writes out the photon list for this obs file at $INTERM_PATH/photonListFileName
        currently cuts out photons outside the valid wavelength ranges from the wavecal
        """
        plFile = self.createEmptyPhotonListFile()
        plTable = plFile.root.photons.photons
        plFile.copyNode(self.flatCalFile.root.flatcal,newparent=plFile.root,recursive=True)
        plFile.copyNode(self.fluxCalFile.root.fluxcal,newparent=plFile.root,recursive=True)
        plFile.copyNode(self.wvlCalFile.root.wavecal,newparent=plFile.root,recursive=True)
        plFile.copyNode(self.file.root.header,newparent=plFile.root,recursive=True)
        plFile.flush()

        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                flag = self.wvlFlagTable[iRow,iCol]
                if flag == 0:#only write photons in good pixels
                    wvlError = self.wvlErrorTable[iRow,iCol]
                    flatWeights = self.flatWeights[iRow,iCol]
                    flatFlags = self.flatFlags[iRow,iCol]
                    fluxWeights = self.fluxWeights[iRow,iCol]
                    fluxFlags = self.fluxFlags[iRow,iCol]
                    wvlRange = self.wvlRangeTable[iRow,iCol]

                    #go through the list of seconds in a pixel dataset
                    for iSec,secData in enumerate(self.getPixel(iRow,iCol)):
                        timestamps,parabolaPeaks,baselines = self.parsePhotonPackets(secData)
                        pulseHeights = np.array(parabolaPeaks,dtype='double') - np.array(baselines,dtype='double')
                        timestamps = iSec + self.tickDuration*timestamps
                        wavelengths = self.convertToWvl(pulseHeights,iRow,iCol,excludeBad=False)
                        #calculate what wavelength bins each photon falls into to see which flat cal factor should be applied
                        if len(wavelengths) > 0:
                            binIndices = np.digitize(wavelengths,self.flatCalWvlBins)-1

                        else:
                            binIndices = np.array([])
 
                        for iPhoton in xrange(len(timestamps)):
                            if wavelengths[iPhoton] > wvlRange[0] and wavelengths[iPhoton] < wvlRange[1] and binIndices[iPhoton] >= 0 and binIndices[iPhoton] < len(flatWeights):
                                #create a new row for the photon list
                                newRow = plTable.row
                                newRow['Xpix'] = iCol
                                newRow['Ypix'] = iRow
                                newRow['ArrivalTime'] = timestamps[iPhoton]
                                newRow['Wavelength'] = wavelengths[iPhoton]
                                newRow['WaveError'] = wvlError
                                newRow['Flag'] = flatFlags[binIndices[iPhoton]]
                                newRow['FlatWeight'] = flatWeights[binIndices[iPhoton]]
                                newRow['FluxWeight'] = fluxWeights[binIndices[iPhoton]]
                                newRow['FluxFlag'] = fluxFlags[binIndices[iphoton]]
                                newRow.append()
        plTable.flush()

    def createEmptyPhotonListFile(self):
        """
        creates a photonList h5 file 
        using header in headers.ArconsHeaders
        """
            
        fileTimestamp = self.fileName.split('_')[1].split('.')[0]
        fileDate = os.path.basename(os.path.dirname(self.fullFileName))
        run = os.path.basename(os.path.dirname(os.path.dirname(self.fullFileName)))
        fn = FileName(run=run,date=fileDate,tstamp=fileTimestamp)
        fullPhotonListFileName = fn.photonList()
        if (os.path.exists(fullPhotonListFileName)):
            if utils.confirm('Photon list file  %s exists. Overwrite?'%fullPhotonListFileName,defaultResponse=False) == False:
                exit(0)
        zlibFilter = tables.Filters(complevel=1,complib='zlib',fletcher32=False)
        try:
            plFile = tables.openFile(fullPhotonListFileName,mode='w')
            plGroup = plFile.createGroup('/','photons','Group containing photon list')
            plTable = plFile.createTable(plGroup,'photons',ArconsHeaders.PhotonList,'Photon List Data',filters=zlibFilter)
        except:
            plFile.close()
            raise
        return plFile



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
        self.wvlRangeTable = np.zeros([self.nRow,self.nCol,2])
        for calPixel in wvlCalData:
            self.wvlFlagTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['wave_flag']
            self.wvlErrorTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['sigma']
            if calPixel['wave_flag'] == 0:
                self.wvlCalTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['polyfit']
                self.wvlRangeTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['solnrange']

        
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
        self.flatFlags = self.flatCalFile.root.flatcal.flags.read()
        self.flatCalWvlBins = self.flatCalFile.root.flatcal.wavelengthBins.read()
        self.nFlatCalWvlBins = self.flatWeights.shape[2]

    def loadFluxCalFile(self,fluxCalFileName):
        """
        loads the flux cal factors from the given file
        """
        scratchDir = os.getenv('INTERM_PATH','/')
        fluxCalPath = os.path.join(scratchDir,'fluxCalSolnFiles')
        fullFluxCalFileName = os.path.join(fluxCalPath,fluxCalFileName)
        if (not os.path.exists(fullFluxCalFileName)):
            print 'flux cal file does not exist: ',fullFluxCalFileName
            return
        self.fluxCalFile = tables.openFile(fullFluxCalFileName,mode='r')
        self.fluxWeights = self.fluxCalFile.root.fluxcal.weights.read()
        self.fluxFlags = self.fluxCalFile.root.fluxcal.flags.read()
        self.fluxCalWvlBins = self.fluxCalFile.root.fluxcal.wavelengthBins.read()
        self.nFluxCalWvlBins = self.nFlatCalWvlBins

    def getDeadPixels(self,showMe=False,weighted=True):
        """
        returns a mask indicating which pixels had no counts in this observation file
        1's for pixels with counts, 0's for pixels without counts
        if showMe is True, a plot of the mask pops up
        """
        countArray = np.array([[self.getPixelCount(iRow,iCol,weighted=weighted) for iCol in range(self.nCol)] for iRow in range(self.nRow)])
        deadArray = np.ones((self.nRow,self.nCol))
        deadArray[countArray == 0] = 0
        if showMe == True:
            utils.plotArray(deadArray)
        return deadArray 

    def getNonAllocPixels(self,showMe=False):
        """
        returns a mask indicating which pixels had no beammap locations 
        (set to constant /r0/p250/)
        1's for pixels with locations, 0's for pixels without locations
        if showMe is True, a plot of the mask pops up
        """
        nonAllocArray = np.ones((self.nRow,self.nCol))
        nonAllocArray[np.core.defchararray.startswith(self.beamImage,self.nonAllocPixelName)] = 0
        if showMe == True:
            utils.plotArray(nonAllocArray)
        return nonAllocArray

    @staticmethod
    def makeWvlBins(energyBinWidth=.1,wvlStart=3000,wvlStop=13000):
        """
        returns an array of wavlength bin edges, with a fixed energy bin width
        withing the limits given in wvlStart and wvlStop
        Args:
            energyBinWidth: bin width in eV
            wvlStart: Lower wavelength edge in Angstrom
            wvlStop: Upper wavelength edge in Angstrom
        Returns:
            an array of wavelength bin edges that can be used with numpy.histogram(bins=wvlBinEdges)
        """

        #Calculate upper and lower energy limits from wavelengths
        #Note that start and stop switch when going to energy
        energyStop = ObsFile.h * ObsFile.c * ObsFile.angstromPerMeter/wvlStart
        energyStart = ObsFile.h * ObsFile.c * ObsFile.angstromPerMeter/wvlStop
        nWvlBins = int((energyStop - energyStart)/energyBinWidth)
        #Construct energy bin edges
        energyBins = np.linspace(energyStart,energyStop,nWvlBins+1)
        #Convert back to wavelength and reverse the order to get increasing wavelengths
        wvlBinEdges = np.array(ObsFile.h * ObsFile.c * ObsFile.angstromPerMeter/energyBins)
        wvlBinEdges = wvlBinEdges[::-1]
        return wvlBinEdges

def calculateSlices(inter, timestamps):
    """
    return a list of strings, with format i0:i1 for a python array slice
    inter is the interval of values in timestamps to mask out.
    The resulting list of strings indicate elements that are not masked out
    """
    slices = []
    prevInclude = not (timestamps[0] in inter)
    if prevInclude:
        slce = "0:"
    nIncluded = 0
    for i in range(len(timestamps)):
        include = not (timestamps[i] in inter)
        if include:
            nIncluded += 1
        if (prevInclude and not include):
            slce += str(i) # end the current range
            slices.append(slce)
        elif (include and not prevInclude):
            slce = "%d:" % i # begin a new range
        prevInclude = include
    if (prevInclude):
        slce += str(i+1) # end the last range if it is included
        slices.append(slce)
    return slices

def repackArray(array, slices):
    """
    returns a copy of array that includes only the element defined by slices
    """
    nIncluded = 0
    for slce in slices:
        s0 = int(slce.split(":")[0])
        s1 = int(slce.split(":")[1])
        nIncluded += s1 - s0
    retval = np.zeros(nIncluded)
    iPt = 0;
    for slce in slices:
        s0 = int(slce.split(":")[0])
        s1 = int(slce.split(":")[1])
        iPtNew = iPt+s1-s0        
        retval[iPt:iPtNew] = array[s0:s1]
        iPt = iPtNew
    return retval
