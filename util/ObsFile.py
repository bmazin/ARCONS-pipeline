#!/bin/python
'''
Author: Matt Strader        Date: August 19, 2012

The class ObsFile is an interface to observation files.  It provides methods for typical ways of accessing and viewing observation data.  It can also load and apply wavelength and flat calibration.  With calibrations loaded, it can write the obs file out as a photon list

Looks for observation files in $MKID_DATA_DIR and calibration files organized in $INTERM_PATH (intermediate or scratch path)

Class Obsfile:
__init__(self, fileName,verbose=False)
__del__(self)
__iter__(self)
loadFile(self, fileName,verbose=False)
checkIntegrity(self,firstSec=0,integrationTime=-1)
convertToWvl(self, pulseHeights, iRow, iCol, excludeBad=True)
createEmptyPhotonListFile(self)
displaySec(self, firstSec=0, integrationTime= -1, weighted=False,fluxWeighted=False, plotTitle='', nSdevMax=2,scaleByEffInt=False)
getFromHeader(self, name)
getPixel(self, iRow, iCol, firstSec=0, integrationTime= -1)
getPixelWvlList(self,iRow,iCol,firstSec=0,integrationTime=-1,excludeBad=True,dither=True)
getPixelCount(self, iRow, iCol, firstSec=0, integrationTime= -1,weighted=False, fluxWeighted=False, getRawCount=False)
getPixelPacketList(self, iRow, iCol, firstSec=0, integrationTime= -1)
getTimedPacketList_old(self, iRow, iCol, firstSec=0, integrationTime= -1)
getTimedPacketList(self, iRow, iCol, firstSec=0, integrationTime= -1)
getPixelCountImage(self, firstSec=0, integrationTime= -1, weighted=False,fluxWeighted=False, getRawCount=False,scaleByEffInt=False)
getAperturePixelCountImage(self, firstSec=0, integrationTime= -1, y_values=range(46), x_values=range(44), y_sky=[], x_sky=[], apertureMask=np.ones((46,44)), skyMask=np.zeros((46,44)), weighted=False, fluxWeighted=False, getRawCount=False, scaleByEffInt=False)
getSpectralCube(self,firstSec=0,integrationTime=-1,weighted=True,wvlStart=3000,wvlStop=13000,wvlBinWidth=None,energyBinWidth=None,wvlBinEdges=None)
getPixelSpectrum(self, pixelRow, pixelCol, firstSec=0, integrationTime= -1,weighted=False, fluxWeighted=False, wvlStart=3000, wvlStop=13000, wvlBinWidth=None, energyBinWidth=None, wvlBinEdges=None)
getPixelBadTimes(self, pixelRow, pixelCol)
getDeadPixels(self, showMe=False, weighted=True)
getNonAllocPixels(self, showMe=False)
getRoachNum(self,iRow,iCol)
getFrame(self, firstSec=0, integrationTime=-1)
loadCentroidListFile(self, centroidListFileName)
loadFlatCalFile(self, flatCalFileName)
loadFluxCalFile(self, fluxCalFileName)
loadHotPixCalFile(self, hotPixCalFileName, switchOnMask=True)
loadTimeAdjustmentFile(self,timeAdjustFileName,verbose=False)
loadWvlCalFile(self, wvlCalFileName)
makeWvlBins(energyBinWidth=.1, wvlStart=3000, wvlStop=13000)
parsePhotonPackets(self, packets, inter=interval(),doParabolaFitPeaks=True, doBaselines=True)
plotPixelSpectra(self, pixelRow, pixelCol, firstSec=0, integrationTime= -1,weighted=False, fluxWeighted=False)getApertureSpectrum(self, pixelRow, pixelCol, radius1, radius2, weighted=False, fluxWeighted=False, lowCut=3000, highCut=7000,firstSec=0,integrationTime=-1)
plotApertureSpectrum(self, pixelRow, pixelCol, radius1, radius2, weighted=False, fluxWeighted=False, lowCut=3000, highCut=7000, firstSec=0,integrationTime=-1)
setWvlCutoffs(self, wvlLowerLimit=3000, wvlUpperLimit=8000)
switchOffHotPixTimeMask(self)
switchOnHotPixTimeMask(self)
writePhotonList(self)

calculateSlices_old(inter, timestamps)
calculateSlices(inter, timestamps)
repackArray(array, slices)
'''

import sys, os
import warnings
import tables
import numpy as np
from numpy import vectorize
from numpy import ma
import matplotlib.pyplot as plt
from util import utils
from interval import interval, inf, imath
from util.FileName import FileName
from scipy import pi
from tables.nodes import filenode
from headers import TimeMask
import time

class ObsFile:
    h = 4.135668e-15 #eV s
    c = 2.998e8 #m/s
    angstromPerMeter = 1e10
    nCalCoeffs = 3
    def __init__(self, fileName, verbose=False):
        """
        load the given file with fileName relative to $MKID_DATA_DIR
        """
        self.loadFile(fileName,verbose=verbose)
        self.wvlCalFile = None #initialize to None for an easy test of whether a cal file has been loaded
        self.flatCalFile = None
        self.fluxCalFile = None
        self.timeAdjustFile = None
        self.hotPixFile = None
        self.hotPixTimeMask = None
        self.hotPixIsApplied = False
        self.cosmicMaskIsApplied = False
        self.cosmicMask = None # interval of times to mask cosmic ray events
        self.centroidListFile = None
        self.wvlLowerLimit = None
        self.wvlUpperLimit = None
        

    def __del__(self):
        """
        Closes the obs file and any cal files that are open
        """
        try:
            self.wvlCalFile.close()
        except:
            pass
        try:
            self.flatCalFile.close()
        except:
            pass
        try:
            self.fluxCalFile.close()
        except:
            pass
        try:
            self.timeAdjustFile.close()
        except:
            pass
        try:
            self.hotPixFile.close()
        except:
            pass
        try:
            self.centroidListFile.close()
        except:
            pass
        self.file.close()


    def __iter__(self):
        """
        Allows easy iteration over pixels in obs file
        use with 'for pixel in obsFileObject:'
        yields a single pixel h5 dataset

        MJS 3/28
        Warning: if timeAdjustFile is loaded, the data from this
        function will not be corrected for roach delays as in getPixel().
        Use getPixel() instead.
        """
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                pixelLabel = self.beamImage[iRow][iCol]
                pixelData = self.file.getNode('/' + pixelLabel)
                yield pixelData

    def loadFile(self, fileName,verbose=False):
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
            dataDir = os.getenv('MKID_DATA_DIR', '/')
            self.fullFileName = os.path.join(dataDir, self.fileName)

        if (not os.path.exists(self.fullFileName)):
            msg='file does not exist: %s'%self.fullFileName
            if verbose:
                print msg
            raise Exception(msg)
        
        #open the hdf5 file
        self.file = tables.openFile(self.fullFileName, mode='r')

        #get the header
        self.header = self.file.root.header.header
        self.titles = self.header.colnames
        try:
            self.info = self.header[0] #header is a table with one row
        except IndexError as inst:
            if verbose:
                print 'Can\'t read header for ',self.fullFileName
            raise inst

        # Useful information about data format set here.
        # For now, set all of these as constants.
        # If we get data taken with different parameters, straighten
        # that all out here.

        ## These parameters are for LICK2012 and PAL2012 data
        self.tickDuration = 1e-6 #s
        self.ticksPerSec = int(1.0 / self.tickDuration)
        self.intervalAll = interval[0.0, (1.0 / self.tickDuration) - 1]
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
        self.pulseMask = int(self.nBitsInPulseHeight * '1', 2) 
        #bitmask of 20 ones
        self.timestampMask = int(self.nBitsInTimestamp * '1', 2) 

        #get the beam image.
        try:
            self.beamImage = self.file.getNode('/beammap/beamimage').read()
        except Exception as inst:
            if verbose:
                print 'Can\'t access beamimage for ',self.fullFileName
            raise inst

        beamShape = self.beamImage.shape
        self.nRow = beamShape[0]
        self.nCol = beamShape[1]

    def checkIntegrity(self,firstSec=0,integrationTime=-1):
        """
        Checks the obs file for corrupted end-of-seconds
        Corruption is indicated by timestamps greater than 1/tickDuration=1e6
        returns 0 if no corruption found
        """
        corruptedPixels = []
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                packetList = self.getPixelPacketList(iRow,iCol,firstSec,integrationTime)
                timestamps,parabolaPeaks,baselines = self.parsePhotonPackets(packetList)
                if np.any(timestamps > 1./self.tickDuration):
                    print 'Corruption detected in pixel (',iRow,iCol,')'
                    corruptedPixels.append((iRow,iCol))
        corruptionFound = len(corruptedPixels) != 0
        return corruptionFound
#        exptime = self.getFromHeader('exptime')
#        lastSec = firstSec + integrationTime
#        if integrationTime == -1:
#            lastSec = exptime-1
#            
#        corruptedSecs = []
#        for pixelCoord in corruptedPixels:
#            for sec in xrange(firstSec,lastSec):
#                packetList = self.getPixelPacketList(pixelCoord[0],pixelCoord[1],sec,integrationTime=1)
#                timestamps,parabolaPeaks,baselines = self.parsePhotonPackets(packetList)
#                if np.any(timestamps > 1./self.tickDuration):
#                    pixelLabel = self.beamImage[iRow][iCol]
#                    corruptedSecs.append(sec)
#                    print 'Corruption in pixel',pixelLabel, 'at',sec

    def convertToWvl(self, pulseHeights, iRow, iCol, excludeBad=True):
        """
        applies wavelength calibration to a list of photon pulse heights
        if excludeBad is True, wavelengths calculated as np.inf are excised from the array returned, as are wavelengths outside the fit limits of the wavecal
        """

        xOffset = self.wvlCalTable[iRow, iCol, 0]
        yOffset = self.wvlCalTable[iRow, iCol, 1]
        amplitude = self.wvlCalTable[iRow, iCol, 2]
        wvlCalLowerLimit = self.wvlRangeTable[iRow, iCol, 0]
        wvlCalUpperLimit = self.wvlRangeTable[iRow, iCol, 1]
        energies = amplitude * (pulseHeights - xOffset) ** 2 + yOffset
        
        if excludeBad == True:
            energies = energies[energies != 0]
        wavelengths = ObsFile.h * ObsFile.c * ObsFile.angstromPerMeter / energies
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
    
    
    def createEmptyPhotonListFile(self,*nkwargs,**kwargs):
        """
        creates a photonList h5 file using header in headers.ArconsHeaders
        Shifted functionality to photonlist/photlist.py, JvE May 10 2013.
        See that function for input parameters and outputs.
        """
        import photonlist.photlist      #Here instead of at top to avoid circular imports
        photonlist.photlist.createEmptyPhotonListFile(self,*nkwargs,**kwargs)


#    def createEmptyPhotonListFile(self,fileName=None):
#        """
#        creates a photonList h5 file 
#        using header in headers.ArconsHeaders
#        
#        INPUTS:
#            fileName - string, name of file to write to. If not supplied, default is used
#                       based on name of original obs. file and standard directories etc.
#                       (see usil.FileName). Added 4/29/2013, JvE
#        """
#        
#        if fileName is None:    
#            fileTimestamp = self.fileName.split('_')[1].split('.')[0]
#            fileDate = os.path.basename(os.path.dirname(self.fullFileName))
#            run = os.path.basename(os.path.dirname(os.path.dirname(self.fullFileName)))
#            fn = FileName(run=run, date=fileDate, tstamp=fileTimestamp)
#            fullPhotonListFileName = fn.photonList()
#        else:
#            fullPhotonListFileName = fileName
#        if (os.path.exists(fullPhotonListFileName)):
#            if utils.confirm('Photon list file  %s exists. Overwrite?' % fullPhotonListFileName, defaultResponse=False) == False:
#                exit(0)
#        zlibFilter = tables.Filters(complevel=1, complib='zlib', fletcher32=False)
#        try:
#            plFile = tables.openFile(fullPhotonListFileName, mode='w')
#            plGroup = plFile.createGroup('/', 'photons', 'Group containing photon list')
#            plTable = plFile.createTable(plGroup, 'photons', ArconsHeaders.PhotonList, 'Photon List Data', 
#                                         filters=zlibFilter, 
#                                         expectedrows=300000)  #Temporary fudge to see if it helps!
#        except:
#            plFile.close()
#            raise
#        return plFile
        
    def displaySec(self, firstSec=0, integrationTime= -1, weighted=False,
                   fluxWeighted=False, plotTitle='', nSdevMax=2,
                   scaleByEffInt=False, getRawCount=False, fignum=None):
        """
        plots a time-flattened image of the counts integrated from firstSec to firstSec+integrationTime
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        if scaleByEffInt is True, then counts are scaled by effective exposure
        time on a per-pixel basis.
        nSdevMax - max end of stretch scale for display, in # sigmas above the mean.
        getRawCount - if True the raw non-wavelength-calibrated image is
        displayed with no wavelength cutoffs applied (in which case no wavecal
        file need be loaded).
        fignum - as for utils.plotArray (None = new window; False/0 = current window; or 
                 specify target window number).
        """
        secImg = self.getPixelCountImage(firstSec, integrationTime, weighted, fluxWeighted,
                                         getRawCount=getRawCount,scaleByEffInt=scaleByEffInt)['image']
        utils.plotArray(secImg, cbar=True, normMax=np.mean(secImg) + nSdevMax * np.std(secImg),
                        plotTitle=plotTitle, fignum=fignum)

            
    def getFromHeader(self, name):
        """
        Returns a requested entry from the obs file header
        If asked for exptime (exposure time) and some roaches have a timestamp offset
        The returned exposure time will be shortened by the max offset, since ObsFile
        will not retrieve data from seconds in which some roaches do not have data
        """
        entry = self.info[self.titles.index(name)]
        if name=='exptime' and self.timeAdjustFile != None:
            #shorten the effective exptime by the number of seconds that 
            #does not have data from all roaches
            maxDelay = np.max(self.roachDelays)
            entry -= maxDelay
        if name=='unixtime' and self.timeAdjustFile != None:
            #the way getPixel retrieves data accounts for individual roach delay,
            #but shifted everything by np.max(self.roachDelays), relabeling sec maxDelay as sec 0
            #so, add maxDelay to the header start time, so all times will be correct relative to it
            entry += np.max(self.roachDelays)
            entry += self.firmwareDelay
        return entry
        
    def getPixel(self, iRow, iCol, firstSec=0, integrationTime= -1):
        """
        Retrieves a pixel using the file's attached beammap.
        If firstSec/integrationTime are provided, only data from the time 
        interval 'firstSec' to firstSec+integrationTime are returned.
        For now firstSec and integrationTime can only be integers.
        If integrationTime is -1, all data after firstSec are returned.

        MJS 3/28
        Updated so if timeAdjustFile is loaded, data retrieved from roaches
        with a delay will be offset to match other roaches.  Also, if some roaches
        have a delay, seconds in which some roaches don't have data are no longer
        retrieved
        """
        pixelLabel = self.beamImage[iRow][iCol]
        pixelNode = self.file.getNode('/' + pixelLabel)

        if self.timeAdjustFile != None:
            iRoach = self.getRoachNum(iRow,iCol)
            maxDelay = np.max(self.roachDelays)
            #skip over any seconds that don't have data from all roaches
            #and offset by roach delay so all roaches will match
            firstSec += maxDelay-self.roachDelays[iRoach]

            if integrationTime == -1:
                lastSec = pixelNode.nrows-self.roachDelays[iRoach]
            else:
                lastSec = firstSec + integrationTime
        else:
            if integrationTime == -1:
                lastSec = pixelNode.nrows
            else:
                lastSec = firstSec + integrationTime
            
        pixelData = pixelNode.read(firstSec, lastSec)
        #return {'pixelData':pixelData,'firstSec':firstSec,'lastSec':lastSec}
        return pixelData


    def getPixelWvlList(self,iRow,iCol,firstSec=0,integrationTime=-1,excludeBad=True,dither=True,timeSpacingCut=None): #,getTimes=False):
        """
        returns a numpy array of photon wavelengths for a given pixel, integrated from firstSec to firstSec+integrationTime.
        if integrationTime is -1, All time after firstSec is used. 
        Now always accounts for any hot-pixel time masking and returns a 
        dictionary with keys:
            timestamps
            wavelengths
            effIntTime  (effective integration time)
        JvE 3/5/2013
        if excludeBad is True, relevant wavelength cuts are applied to timestamps and wavelengths before returning 
        [if getTimes is True, returns timestamps,wavelengths - OBSELETED - JvE 3/2/2013]

        MJS 3/28/2013
        if dither is True, uniform random values in the range (0,1) will be added to all quantized ADC values read, to remedy the effects of quantization
        """
        
        #if getTimes == False:
        #    packetList = self.getPixelPacketList(iRow,iCol,firstSec,integrationTime)
        #    timestamps,parabolaPeaks,baselines = self.parsePhotonPackets(packetList)
        #    
        #else:
        
        x = self.getTimedPacketList(iRow, iCol, firstSec, integrationTime,timeSpacingCut=timeSpacingCut)
        timestamps, parabolaPeaks, baselines, effIntTime = \
            x['timestamps'], x['peakHeights'], x['baselines'], x['effIntTime']
        parabolaPeaks = np.array(parabolaPeaks,dtype=np.double)
        baselines = np.array(baselines,dtype=np.double)
        if dither==True:
            parabolaPeaks += np.random.random_sample(len(parabolaPeaks))
            baselines += np.random.random_sample(len(baselines))
                    
        pulseHeights = parabolaPeaks - baselines
        xOffset = self.wvlCalTable[iRow,iCol,0]
        yOffset = self.wvlCalTable[iRow,iCol,1]
        amplitude = self.wvlCalTable[iRow,iCol,2]
        wvlCalLowerLimit = self.wvlRangeTable[iRow,iCol,0]
        wvlCalUpperLimit = self.wvlRangeTable[iRow,iCol,1]
        energies = amplitude*(pulseHeights-xOffset)**2+yOffset
        
        wavelengths = ObsFile.h*ObsFile.c*ObsFile.angstromPerMeter/energies
        if excludeBad == True:
            goodMask = ~np.isnan(wavelengths)
            goodMask = np.logical_and(goodMask,wavelengths!=np.inf)
            if self.wvlLowerLimit == -1:
                goodMask = np.logical_and(goodMask,wvlCalLowerLimit < wavelengths)
            elif self.wvlLowerLimit != None:
                goodMask = np.logical_and(goodMask,self.wvlLowerLimit < wavelengths)
            if self.wvlUpperLimit == -1:
                goodMask = np.logical_and(goodMask,wavelengths < wvlCalUpperLimit)
            elif self.wvlUpperLimit != None:
                goodMask = np.logical_and(goodMask,wavelengths < self.wvlUpperLimit)
            wavelengths = wavelengths[goodMask]
            timestamps = timestamps[goodMask]

        return {'timestamps':timestamps, 'wavelengths':wavelengths,
                'effIntTime':effIntTime}
            
    def getPixelCount(self, iRow, iCol, firstSec=0, integrationTime= -1,
                      weighted=False, fluxWeighted=False, getRawCount=False):
        """
        returns the number of photons received in a given pixel from firstSec to firstSec + integrationTime
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        if fluxWeighted is True, flux weights are applied.
        if getRawCount is True, the total raw count for all photon event detections
        is returned irrespective of wavelength calibration, and with no wavelength
        cutoffs (in this case, no wavecal file need have been applied, though 
        bad pixel time-masks *will* still be applied if present and switched 'on'.) 
        Otherwise will now always call getPixelSpectrum (which is also capable 
        of handling hot pixel removal) -- JvE 3/1/2013.
        Updated to return effective exp. times; see below. -- JvE 3/2013. 
        
        OUTPUTS:
        Return value is a dictionary with tags:
            'counts':int, number of photon counts
            'effIntTime':float, effective integration time after time-masking is 
                     accounted for.
        """
        
        if getRawCount is True:
            x = self.getTimedPacketList(iRow, iCol, firstSec=firstSec, integrationTime=integrationTime)
            #x2 = self.getTimedPacketList_old(iRow, iCol, firstSec=firstSec, integrationTime=integrationTime)
            #assert np.array_equal(x['timestamps'],x2['timestamps'])
            #assert np.array_equal(x['effIntTime'],x2['effIntTime'])
            #assert np.array_equal(x['peakHeights'],x2['peakHeights'])
            #assert np.array_equal(x['baselines'],x2['baselines'])
            timestamps, effIntTime = x['timestamps'], x['effIntTime']
            counts = len(timestamps)
            return {'counts':counts, 'effIntTime':effIntTime}

        else:
            pspec = self.getPixelSpectrum(iRow, iCol, firstSec, integrationTime,weighted=weighted, fluxWeighted=fluxWeighted)
            counts = sum(pspec['spectrum'])
            return {'counts':counts, 'effIntTime':pspec['effIntTime']}

    def getPixelPacketList(self, iRow, iCol, firstSec=0, integrationTime= -1):
        """
        returns a numpy array of 64-bit photon packets for a given pixel, integrated from firstSec to firstSec+integrationTime.
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
#        getPixelOutput = self.getPixel(iRow, iCol, firstSec, integrationTime)
#        pixelData = getPixelOutput['pixelData']
        pixelData = self.getPixel(iRow,iCol,firstSec,integrationTime)
        packetList = np.concatenate(pixelData)
        return packetList

    def getTimedPacketList_old(self, iRow, iCol, firstSec=0, integrationTime= -1):
        """
        DEPRECATED VERSION. JvE 3/13/2013
        
        Parses an array of uint64 packets with the obs file format,and makes timestamps absolute
        inter is an interval of time values to mask out [missing?? To be implented? JvE 2/18/13]
        returns a list of timestamps,parabolaFitPeaks,baselines,effectiveIntTime (effective
        integration time after accounting for time-masking.)
        parses packets from firstSec to firstSec+integrationTime.
        if integrationTime is -1, all time after firstSec is used.  
        
        Now updated to take advantage of masking capabilities in parsePhotonPackets
        to allow for correct application of non-integer values in firstSec and
        integrationTime. JvE Feb 27 2013.
        
        CHANGED RETURN VALUES - now returns a dictionary including effective integration
        time (allowing for bad pixel masking), with keys:
        
            'timestamps'
            'peakHeights'
            'baselines'
            'effIntTime'
         
         - JvE 3/5/2013.
        """
        pixelData = self.getPixel(iRow, iCol)
        lastSec = firstSec + integrationTime
        if integrationTime == -1 or lastSec > len(pixelData):
            lastSec = len(pixelData)
        pixelData = pixelData[int(np.floor(firstSec)):int(np.ceil(lastSec))]

        if self.hotPixIsApplied:
            inter = self.getPixelBadTimes(iRow, iCol)
        else:
            inter = interval()
        
        if (type(firstSec) is not int) or (type(integrationTime) is not int):
            #Also exclude times outside firstSec to lastSec. Allows for sub-second
            #(floating point) values in firstSec and integrationTime
            inter = inter | interval([-np.inf, firstSec], [lastSec, np.inf])   #Union the exclusion interval with the excluded time range limits

        #Inter now contains a single 'interval' instance, which contains a list of
        #times to exclude, in seconds, including all times outside the requested
        #integration if necessary.

        #Calculate the total effective time for the integration after removing
        #any 'intervals':
        integrationInterval = interval([firstSec, lastSec])
        maskedIntervals = inter & integrationInterval  #Intersection of the integration and the bad times for this pixel.
        effectiveIntTime = (lastSec - firstSec) - utils.intervalSize(maskedIntervals)

        timestamps = []
        baselines = []
        peakHeights = []

        for t in range(len(pixelData)):
            interTicks = (inter - np.floor(firstSec) - t) * self.ticksPerSec
            times, peaks, bases = self.parsePhotonPackets(pixelData[t], inter=interTicks)
            times = np.floor(firstSec) + self.tickDuration * times + t
            timestamps.append(times)
            baselines.append(bases)
            peakHeights.append(peaks)
            
        timestamps = np.concatenate(timestamps)
        baselines = np.concatenate(baselines)
        peakHeights = np.concatenate(peakHeights)
        return {'timestamps':timestamps, 'peakHeights':peakHeights,
                'baselines':baselines, 'effIntTime':effectiveIntTime}

    def getTimedPacketList(self, iRow, iCol, firstSec=0, integrationTime= -1, timeSpacingCut=None,expTailTimescale=None):
        """
        Parses an array of uint64 packets with the obs file format,and makes timestamps absolute
        (with zero time at beginning of ObsFile).
        Returns a list of:
            timestamps (seconds from start of file),parabolaFitPeaks,baselines,effectiveIntTime (effective
            integration time after accounting for time-masking.)
        parses packets from firstSec to firstSec+integrationTime.
        if integrationTime is -1, all time after firstSec is used.  
        if timeSpacingCut is not None, photons sooner than timeSpacingCut seconds after the last photon are cut.
            Typically we will set timeSpacingCut=1.e-3 (1 ms) to remove effects of photon pile-up
        if expTailTimescale is not None, photons are assumed to exhibit an exponential decay back to baseline with e-fold time
            expTailTimescale, this is used to subtract the exponential tail of one photon from the peakHeight of the next photon
            This also attempts to counter effects of photon pile-up for short (<100 us) dead times.
        
        Now updated to take advantage of masking capabilities in parsePhotonPackets
        to allow for correct application of non-integer values in firstSec and
        integrationTime. JvE Feb 27 2013.
        
        CHANGED RETURN VALUES - now returns a dictionary including effective integration
        time (allowing for bad pixel masking), with keys:
        
            'timestamps'
            'peakHeights'
            'baselines'
            'effIntTime'
         
         - JvE 3/5/2013.
         
         **Modified to increase speed for integrations shorter than the full exposure
         length. JvE 3/13/2013**

         MJS 3/28/2012
         **Modified to add known delays to timestamps from roach delays and firmware delay if timeAdjustFile 
         is loaded**
         
        """
        #pixelData = self.getPixel(iRow, iCol)
        lastSec = firstSec + integrationTime
        #Make sure we include *all* the complete seconds that overlap the requested range
        integerIntTime = int(np.ceil(lastSec)-np.floor(firstSec)) 
        pixelData = self.getPixel(iRow, iCol, firstSec=int(np.floor(firstSec)),
                                          integrationTime=integerIntTime)

        if integrationTime == -1 or integerIntTime > len(pixelData):
            lastSec = int(np.floor(firstSec))+len(pixelData)

        if self.hotPixIsApplied:
            inter = self.getPixelBadTimes(iRow, iCol)
        else:
            inter = interval()

        if self.cosmicMaskIsApplied:
            inter = inter | self.cosmicMask
            
        if (type(firstSec) is not int) or (type(integrationTime) is not int):
            #Also exclude times outside firstSec to lastSec. Allows for sub-second
            #(floating point) values in firstSec and integrationTime in the call to parsePhotonPackets.
            inter = inter | interval([-np.inf, firstSec], [lastSec, np.inf])   #Union the exclusion interval with the excluded time range limits

        #Inter now contains a single 'interval' instance, which contains a list of
        #times to exclude, in seconds, including all times outside the requested
        #integration if necessary.

        #Calculate the total effective time for the integration after removing
        #any 'intervals':
        integrationInterval = interval([firstSec, lastSec])
        maskedIntervals = inter & integrationInterval  #Intersection of the integration and the bad times for this pixel (for calculating eff. int. time)
        effectiveIntTime = (lastSec - firstSec) - utils.intervalSize(maskedIntervals)

        timestamps = []
        baselines = []
        peakHeights = []

        # pixelData is an array of data for this iRow,iCol, at each good time
        for t in range(len(pixelData)):
            interTicks = (inter - (np.floor(firstSec) + t)) * self.ticksPerSec         
            times, peaks, bases = self.parsePhotonPackets(pixelData[t], inter=interTicks)
            times = np.floor(firstSec) + self.tickDuration * times + t
            timestamps.append(times)
            baselines.append(bases)
            peakHeights.append(peaks)
            
        if len(pixelData) > 0:         #Check that concatenate won't barf (check added JvE, 6/17/2013).
            timestamps = np.concatenate(timestamps)
            baselines = np.concatenate(baselines)
            peakHeights = np.concatenate(peakHeights)
        else:
            timestamps = np.array([])
            baselines = np.array([])
            peakHeights = np.array([])

        if expTailTimescale != None and len(timestamps) > 0:
            #find the time between peaks
            timeSpacing = np.diff(timestamps)
            timeSpacing[timeSpacing < 0] = 1.
            timeSpacing = np.append(1.,timeSpacing)#arbitrarily assume the first photon is 1 sec after the one before it
            relPeakHeights = peakHeights-baselines
            
            #assume each peak is riding on the tail of an exponential starting at the peak before it with e-fold time of expTailTimescale
            print 'dt',timeSpacing[0:10]
            expTails = (1.*peakHeights-baselines)*np.exp(-1.*timeSpacing/expTailTimescale)
            print 'expTail',expTails[0:10]
            print 'peak',peakHeights[0:10]
            print 'peak-baseline',1.*peakHeights[0:10]-baselines[0:10]
            print 'expT',np.exp(-1.*timeSpacing[0:10]/expTailTimescale)
            #subtract off this exponential tail
            peakHeights = np.array(peakHeights-expTails,dtype=np.int)
            print 'peak',peakHeights[0:10]
                
            
        if timeSpacingCut != None and len(timestamps) > 0:
            timeSpacing = np.diff(timestamps)
            timeSpacingMask = np.concatenate([[True],timeSpacing >= timeSpacingCut]) #include first photon and photons after who are at least timeSpacingCut after the previous photon
            timestamps = timestamps[timeSpacingMask]
            peakHeights = peakHeights[timeSpacingMask]
            baselines = baselines[timeSpacingMask]


        return {'timestamps':timestamps, 'peakHeights':peakHeights,
                'baselines':baselines, 'effIntTime':effectiveIntTime}

    def getPixelCountImage(self, firstSec=0, integrationTime= -1, weighted=False,
                           fluxWeighted=False, getRawCount=False,
                           scaleByEffInt=False):
        """
        Return a time-flattened image of the counts integrated from firstSec to firstSec+integrationTime.
        If integration time is -1, all time after firstSec is used.
        If weighted is True, flat cal weights are applied. JvE 12/28/12
        If fluxWeighted is True, flux cal weights are applied. SM 2/7/13
        If getRawCount is True then the raw non-wavelength-calibrated image is
        returned with no wavelength cutoffs applied (in which case no wavecal
        file need be loaded). JvE 3/1/13
        If scaleByEffInt is True, any pixels that have 'bad' times masked out
        will have their counts scaled up to match the equivalent integration 
        time requested.
        RETURNS:
            Dictionary with keys:
                'image' - a 2D array representing the image
                'effIntTimes' - a 2D array containing effective integration 
                                times for each pixel.
        """
        secImg = np.zeros((self.nRow, self.nCol))
        effIntTimes = np.zeros((self.nRow, self.nCol), dtype=np.float64)
        effIntTimes.fill(np.nan)   #Just in case an element doesn't get filled for some reason.
        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                pcount = self.getPixelCount(iRow, iCol, firstSec, integrationTime,
                                          weighted, fluxWeighted, getRawCount)
                secImg[iRow, iCol] = pcount['counts']
                effIntTimes[iRow, iCol] = pcount['effIntTime']
        if scaleByEffInt is True:
            if integrationTime == -1:
                totInt = self.getFromHeader('exptime')
            else:
                totInt = integrationTime
            secImg *= (totInt / effIntTimes)                    
        #if getEffInt is True:
        return{'image':secImg, 'effIntTimes':effIntTimes}
        #else:
        #    return secImg

    def getAperturePixelCountImage(self, firstSec=0, integrationTime= -1, y_values=range(46), x_values=range(44), y_sky=[], x_sky=[], apertureMask=np.ones((46,44)), skyMask=np.zeros((46,44)), weighted=False, fluxWeighted=False, getRawCount=False, scaleByEffInt=False):

        """
        Return a time-flattened image of the counts integrated from firstSec to firstSec+integrationTime 
        This aperture version subtracts out the average sky counts/pixel and includes scaling due to circular apertures. GD 5/27/13
        If integration time is -1, all time after firstSec is used.
        If weighted is True, flat cal weights are applied. JvE 12/28/12
        If fluxWeighted is True, flux cal weights are applied. SM 2/7/13
        If getRawCount is True then the raw non-wavelength-calibrated image is
        returned with no wavelength cutoffs applied (in which case no wavecal
        file need be loaded). JvE 3/1/13
        If scaleByEffInt is True, any pixels that have 'bad' times masked out
        will have their counts scaled up to match the equivalent integration 
        time requested.
        RETURNS:
            Dictionary with keys:
                'image' - a 2D array representing the image
                'effIntTimes' - a 2D array containing effective integration 
                                times for each pixel.
        """
        secImg = np.zeros((self.nRow, self.nCol))
        effIntTimes = np.zeros((self.nRow, self.nCol), dtype=np.float64)
        effIntTimes.fill(np.nan)   #Just in case an element doesn't get filled for some reason.
        skyValues=[]
        objValues=[]
        AreaSky=[]
        AreaObj=[]
        for pix in xrange(len(y_sky)):
            pcount = self.getPixelCount(y_sky[pix], x_sky[pix], firstSec, integrationTime,weighted, fluxWeighted, getRawCount)
            skyValue=pcount['counts']*skyMask[y_sky[pix]][x_sky[pix]]
            skyValues.append(skyValue)
            AreaSky.append(skyMask[y_sky[pix]][x_sky[pix]])
        skyCountPerPixel = np.sum(skyValues)/(np.sum(AreaSky))
#        print 'sky count per pixel =',skyCountPerPixel
        for pix in xrange(len(y_values)):
            pcount = self.getPixelCount(y_values[pix], x_values[pix], firstSec, integrationTime,weighted, fluxWeighted, getRawCount)
            secImg[y_values[pix],x_values[pix]] = (pcount['counts']-skyCountPerPixel)*apertureMask[y_values[pix]][x_values[pix]]
            AreaObj.append(apertureMask[y_values[pix]][x_values[pix]])
            effIntTimes[y_values[pix],x_values[pix]] = pcount['effIntTime']
            objValues.append(pcount['counts']*apertureMask[y_values[pix]][x_values[pix]])
        AveObj=np.sum(objValues)/(np.sum(AreaObj))
#        print 'ave obj per pixel (not sub) = ',AveObj
        NumObjPhotons = np.sum(secImg)
#        print 'lightcurve = ',NumObjPhotons
        if scaleByEffInt is True:
            secImg *= (integrationTime / effIntTimes)                    
        #if getEffInt is True:
        return{'image':secImg, 'effIntTimes':effIntTimes, 'SkyCountSubtractedPerPixel':skyCountPerPixel,'lightcurve':NumObjPhotons}
        #else:
        #    return secImg
    
    def getSpectralCube(self,firstSec=0,integrationTime=-1,weighted=True,wvlStart=3000,wvlStop=13000,wvlBinWidth=None,energyBinWidth=None,wvlBinEdges=None,timeSpacingCut=None):
        """
        Return a time-flattened spectral cube of the counts integrated from firstSec to firstSec+integrationTime.
        If integration time is -1, all time after firstSec is used.
        If weighted is True, flat cal weights are applied.
        """
        cube = [[[] for iCol in range(self.nCol)] for iRow in range(self.nRow)]
        effIntTime = np.zeros((self.nRow,self.nCol))

        for iRow in xrange(self.nRow):
            for iCol in xrange(self.nCol):
                x = self.getPixelSpectrum(pixelRow=iRow,pixelCol=iCol,
                                  firstSec=firstSec,integrationTime=integrationTime,
                                  weighted=weighted,wvlStart=wvlStart,wvlStop=wvlStop,
                                  wvlBinWidth=wvlBinWidth,energyBinWidth=energyBinWidth,
                                  wvlBinEdges=wvlBinEdges,timeSpacingCut=timeSpacingCut)
                cube[iRow][iCol] = x['spectrum']
                effIntTime[iRow][iCol] = x['effIntTime']
                wvlBinEdges = x['wvlBinEdges']
        cube = np.array(cube)
        return {'cube':cube,'wvlBinEdges':wvlBinEdges,'effIntTime':effIntTime}

    def getPixelSpectrum(self, pixelRow, pixelCol, firstSec=0, integrationTime= -1,
                         weighted=False, fluxWeighted=False, wvlStart=3000, wvlStop=13000,
                         wvlBinWidth=None, energyBinWidth=None, wvlBinEdges=None,timeSpacingCut=None):
        """
        returns a spectral histogram of a given pixel integrated from firstSec to firstSec+integrationTime,
        and an array giving the cutoff wavelengths used to bin the wavelength values
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        if weighted is False, flat cal weights are not applied
        the wavelength bins used depends on the parameters given.
        If energyBinWidth is specified, the wavelength bins use fixed energy bin widths
        If wvlBinWidth is specified, the wavelength bins use fixed wavelength bin widths
        If neither is specified and/or if weighted is True, the flat cal wvlBinEdges is used
        
        ----
        Updated to return effective integration time for the pixel
        Returns dictionary with keys:
            'spectrum' - spectral histogram of given pixel.
            'wvlBinEdges' - edges of wavelength bins
            'effIntTime' - the effective integration time for the given pixel 
                           after accounting for hot-pixel time-masking.
        JvE 3/5/2013
        ----
        """
        x = self.getPixelWvlList(pixelRow, pixelCol, firstSec, integrationTime,timeSpacingCut=timeSpacingCut)
        wvlList, effIntTime = x['wavelengths'], x['effIntTime']
        
        if self.flatCalFile != None and ((wvlBinEdges == None and energyBinWidth == None and wvlBinWidth == None) or weighted == True):
        #We've loaded a flat cal already, which has wvlBinEdges defined, and no other bin edges parameters are specified to override it.
            spectrum, wvlBinEdges = np.histogram(wvlList, bins=self.flatCalWvlBins)
            if weighted == True:#Need to apply flat weights by wavelenth
                spectrum = spectrum * self.flatWeights[pixelRow, pixelCol]
                if fluxWeighted == True:
                    spectrum = spectrum * self.fluxWeights
        else:
            if weighted == True:
                raise ValueError('when weighted=True, flatCal wvl bins are used, so wvlBinEdges,wvlBinWidth,energyBinWidth,wvlStart,wvlStop should not be specified')
            if wvlBinEdges == None:#We need to construct wvlBinEdges array
                if energyBinWidth != None:#Fixed energy binwidth specified
                    #Construct array with variable wvl binwidths
                    wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth=energyBinWidth, wvlStart=wvlStart, wvlStop=wvlStop)
                    spectrum, wvlBinEdges = np.histogram(wvlList, bins=wvlBinEdges)
                elif wvlBinWidth != None:#Fixed wvl binwidth specified
                    nWvlBins = int((wvlStop - wvlStart) / wvlBinWidth)
                    spectrum, wvlBinEdges = np.histogram(wvlList, bins=nWvlBins, range=(wvlStart, wvlStop))
                else:
                    raise ValueError('getPixelSpectrum needs either wvlBinWidth,wvlBinEnergy, or wvlBinEdges')
            else:#We are given wvlBinEdges array
                spectrum, wvlBinEdges = np.histogram(wvlList, bins=wvlBinEdges)
       
        #if getEffInt is True:
        return {'spectrum':spectrum, 'wvlBinEdges':wvlBinEdges, 'effIntTime':effIntTime}
        #else:
        #    return spectrum,wvlBinEdges
    
    def getApertureSpectrum(self, pixelRow, pixelCol, radius1, radius2, weighted=False,
                            fluxWeighted=False, lowCut=3000, highCut=7000,firstSec=0,integrationTime=-1):
    	'''
    	Creates a spectrum from a group of pixels.  Aperture is defined by pixelRow and pixelCol of
    	center, as well as radius.  Wave and flat cals should be loaded before using this
    	function.  If no hot pixel mask is applied, taking the median of the sky rather than
    	the average to account for high hot pixel counts.
    	Will add more options as other pieces of pipeline become more refined.
    	(Note - not updated to handle loaded hot pixel time-masks - if applied,
    	behaviour may be unpredictable. JvE 3/5/2013).
    	'''
    	print 'Creating dead pixel mask...'
    	deadMask = self.getDeadPixels()
    	print 'Creating wavecal solution mask...'
    	bad_solution_mask = np.zeros((self.nRow, self.nCol))
    	for y in range(self.nRow):
    	    for x in range(self.nCol):
    		if (self.wvlRangeTable[y][x][0] > lowCut or self.wvlRangeTable[y][x][1] < highCut):
    		    bad_solution_mask[y][x] = 1
    	print 'Creating aperture mask...'
    	apertureMask = utils.aperture(pixelCol, pixelRow, radius=radius1)
    	print 'Creating sky mask...'
    	bigMask = utils.aperture(pixelCol, pixelRow, radius=radius2)
    	skyMask = bigMask - apertureMask
    	#if hotPixMask == None:
    	#    y_values, x_values = np.where(np.logical_and(bad_solution_mask == 0, np.logical_and(apertureMask == 0, deadMask == 1)))
    	#    y_sky, x_sky = np.where(np.logical_and(bad_solution_mask == 0, np.logical_and(skyMask == 0, deadMask == 1)))
    	#else:
    	#    y_values, x_values = np.where(np.logical_and(bad_solution_mask == 0, np.logical_and(np.logical_and(apertureMask == 0, deadMask == 1), hotPixMask == 0)))
    	#    y_sky, x_sky = np.where(np.logical_and(bad_solution_mask == 0, np.logical_and(np.logical_and(skyMask == 0, deadMask == 1), hotPixMask == 0)))

        y_values, x_values = np.where(np.logical_and(bad_solution_mask == 0, np.logical_and(apertureMask == 0, deadMask == 1)))
    	y_sky, x_sky = np.where(np.logical_and(bad_solution_mask == 0, np.logical_and(skyMask == 0, deadMask == 1)))

    	#wvlBinEdges = self.getPixelSpectrum(y_values[0], x_values[0], weighted=weighted)['wvlBinEdges']
    	print 'Creating average sky spectrum...'
    	skyspectrum = []
    	for i in range(len(x_sky)):
            specDict = self.getPixelSpectrum(y_sky[i],x_sky[i],weighted=weighted, fluxWeighted=fluxWeighted, firstSec=firstSec, integrationTime=integrationTime)
            self.skySpectrumSingle,wvlBinEdges,self.effIntTime = specDict['spectrum'],specDict['wvlBinEdges'],specDict['effIntTime']
            self.scaledSpectrum = self.skySpectrumSingle/self.effIntTime #scaled spectrum by effective integration time
            #print "Sky spectrum"
            #print self.skySpectrumSingle
            #print "Int time"
            #print self.effIntTime
    	    skyspectrum.append(self.scaledSpectrum)
    	sky_array = np.zeros(len(skyspectrum[0]))
    	for j in range(len(skyspectrum[0])):
    	    ispectrum = np.zeros(len(skyspectrum))
    	    for i in range(len(skyspectrum)):    
    	        ispectrum[i] = skyspectrum[i][j]
            sky_array[j] = np.median(ispectrum)
    	    #if hotPixMask == None:
    	    #    sky_array[j] = np.median(ispectrum)
    	    #else:
    	    #    sky_array[j] = np.average(ispectrum)
    	print 'Creating sky subtracted spectrum...'
    	spectrum = []
    	for i in range(len(x_values)):
            specDict = self.getPixelSpectrum(y_values[i],x_values[i],weighted=weighted, fluxWeighted=fluxWeighted, firstSec=firstSec, integrationTime=integrationTime)
            self.obsSpectrumSingle,wvlBinEdges,self.effIntTime = specDict['spectrum'],specDict['wvlBinEdges'],specDict['effIntTime']   
            self.scaledSpectrum = self.obsSpectrumSingle/self.effIntTime #scaled spectrum by effective integration time
    	    spectrum.append(self.scaledSpectrum - sky_array)

    	    #spectrum.append(self.getPixelSpectrum(y_values[i], x_values[i], weighted=weighted,fluxWeighted=fluxWeighted)['spectrum'] - sky_array)
    	summed_array = np.zeros(len(spectrum[0]))
    	for j in range(len(spectrum[0])):
    	    ispectrum = np.zeros(len(spectrum))
    	    for i in range(len(spectrum)):    
    	        ispectrum[i] = spectrum[i][j]
    	    summed_array[j] = np.sum(ispectrum)
    	for i in range(len(summed_array)):
    	    summed_array[i] /= (wvlBinEdges[i + 1] - wvlBinEdges[i])
    	return summed_array, wvlBinEdges
    
    def getPixelBadTimes(self, pixelRow, pixelCol):
        """
        Get the time interval(s) for which a given pixel is bad (hot/cold,
        whatever, from the hot pixel cal file).
        Returns an 'interval' object (see pyinterval) of bad times (in seconds
        from start of obs file).
        """
        if self.hotPixTimeMask is None:
            raise RuntimeError, 'No hot pixel file loaded'
        inter = interval.union(self.hotPixTimeMask['intervals'][pixelRow, pixelCol])
        return inter

    def getDeadPixels(self, showMe=False, weighted=True):
        """
        returns a mask indicating which pixels had no counts in this observation file
        1's for pixels with counts, 0's for pixels without counts
        if showMe is True, a plot of the mask pops up
        """
        countArray = np.array([[(self.getPixelCount(iRow, iCol, weighted=weighted))['counts'] for iCol in range(self.nCol)] for iRow in range(self.nRow)])
        deadArray = np.ones((self.nRow, self.nCol))
        deadArray[countArray == 0] = 0
        if showMe == True:
            utils.plotArray(deadArray)
        return deadArray 

    def getNonAllocPixels(self, showMe=False):
        """
        returns a mask indicating which pixels had no beammap locations 
        (set to constant /r0/p250/)
        1's for pixels with locations, 0's for pixels without locations
        if showMe is True, a plot of the mask pops up
        """
        nonAllocArray = np.ones((self.nRow, self.nCol))
        nonAllocArray[np.core.defchararray.startswith(self.beamImage, self.nonAllocPixelName)] = 0
        if showMe == True:
            utils.plotArray(nonAllocArray)
        return nonAllocArray

    def getRoachNum(self,iRow,iCol):
        pixelLabel = self.beamImage[iRow][iCol]
        iRoach = int(pixelLabel.split('r')[1][0])
        return iRoach

    def getFrame(self, firstSec=0, integrationTime=-1):
        """
        return a 2d array of numbers with the integrated flux per pixel,
        suitable for use as a frame in util/utils.py function makeMovie

        firstSec=0 is the starting second to include

        integrationTime=-1 is the number of seconds to include, or -1
        to include all to the end of this file


        output: the frame, in photons per pixel, a 2d numpy array of
        np.unint32

        """
        frame = np.zeros((self.nRow,self.nCol),dtype=np.uint32)
        for iRow in range(self.nRow):
            for iCol in range(self.nCol):
                pl = self.getTimedPacketList(iRow,iCol,
                                             firstSec,integrationTime)
                nphoton = pl['timestamps'].size
                frame[iRow][iCol] += nphoton
        return frame
    
    # a different way to get, with the functionality of getTimedPacketList
    def getPackets(self, iRow, iCol, firstSec, integrationTime, fields=()):
        """
        get and parse packets for pixel iRow,iCol starting at firstSec for integrationTime seconds.

        files is a list of strings to indicate what to parse in addition to timestamps
        (All it know how to do right now is peakHeights)

        return a dictionary containing errectiveIntTime (n seconds), the array of timestamps and other fields requested
        """
        parse = {'peakHeights': True, 'baselines': True}
        for key in parse.keys():
            try:
                fields.index(key)
            except ValueError:
                parse[key] = False

        lastSec = firstSec+integrationTime
        # Work out inter, the times to mask
        # start with nothing being masked
        inter = interval()
        # mask the hot pixels if requested
        if self.hotPixIsApplied:
            inter = self.getPixelBadTimes(iRow, iCol)
        # mask cosmics if requested
        if self.cosmicMaskIsApplied:
            inter = inter | self.cosmicMask

        # mask the fraction of the first integer second not requested
        firstSecInt = int(np.floor(firstSec))
        if (firstSec > firstSecInt):
            inter = inter | interval([firstSecInt, firstSec])
        # mask the fraction of the last integer second not requested
        lastSecInt = int(np.ceil(firstSec+integrationTime))
        integrationTimeInt = lastSecInt-firstSecInt
        if (lastSec < lastSecInt):
            inter = inter | interval([lastSec, lastSecInt])
        #Calculate the total effective time for the integration after removing
        #any 'intervals':
        integrationInterval = interval([firstSec, lastSec])
        maskedIntervals = inter & integrationInterval  #Intersection of the integration and the bad times for this pixel.
        effectiveIntTime = integrationTime - utils.intervalSize(maskedIntervals)

        pixelData = self.getPixel(iRow, iCol, firstSec=firstSecInt, 
                                  integrationTime=integrationTimeInt)
        # calculate how long a np array needs to be to hold everything
        nPackets = 0
        for pixels in pixelData:
            nPackets += len(pixels)

        # create empty arrays
        timestamps = np.empty(nPackets, dtype=np.float)
        if parse['peakHeights']: peakHeights=np.empty(nPackets, np.int16)
        if parse['baselines']: baselines=np.empty(nPackets, np.int16)

        # fill in the arrays one second at a time
        ipt = 0
        t = firstSecInt
        for pixels in pixelData:
            iptNext = ipt+len(pixels)
            timestamps[ipt:iptNext] = \
                t + np.bitwise_and(pixels,self.timestampMask)*self.tickDuration
            if parse['peakHeights']:
                peakHeights[ipt:iptNext] = np.bitwise_and(
                    np.right_shift(pixels, self.nBitsAfterParabolaPeak), 
                    self.pulseMask)

            if parse['baselines']:
                baselines[ipt:iptNext] = np.bitwise_and(
                    np.right_shift(pixels, self.nBitsAfterBaseline), 
                    self.pulseMask)

            ipt = iptNext
            t += 1
        # create a mask, "True" mean mask value
        start = time.clock()
        # the call the makeMask dominates the running time
        mask = ObsFile.makeMask(timestamps, inter)
        elapsed = (time.clock() - start)
        #print "In Obsfile:  elapsed due to makeMask=",elapsed
        # compress out all the masked values
        tsMaskedArray = ma.array(timestamps,mask=mask)
        timestamps = ma.compressed(tsMaskedArray)
        #timestamps = ma.compressed(ma.array(timestamps,mask))
        # build up the dictionary of values and return it
        retval =  {"effIntTime": effectiveIntTime,
                "timestamps":timestamps}
        if parse['peakHeights']: 
            retval['peakHeights'] = \
                ma.compressed(ma.array(peakHeights,mask=mask))
        if parse['baselines']: 
            retval['baselines'] = \
                ma.compressed(ma.array(baselines,mask=mask))
        return retval

    @staticmethod
    def makeMask01(timestamps, inter):
        def myfunc(x): return inter.__contains__(x)
        vecfunc = vectorize(myfunc,otypes=[np.bool])
        return vecfunc(timestamps)

    @staticmethod
    def makeMask(timestamps, inter):
        """
        return an array of booleans, the same length as timestamps,
        with that value inter.__contains__(timestamps[i])
        """
        retval = np.empty(len(timestamps),dtype=np.bool)
        ainter = np.array(inter)
        t0s = ainter[:,0]
        t1s = ainter[:,1]

        tMin = t0s[0]
        tMax = t1s[-1]
                       
        for i in range(len(timestamps)):
            ts = timestamps[i]
            if ts < tMin:
                retval[i] = False
            elif ts > tMax:
                retval[i] = False
            else:
                tIndex = np.searchsorted(t0s, ts)
                t0 = t0s[tIndex-1]
                t1 = t1s[tIndex-1]
                if ts < t1:
                    retval[i] = True
                else:
                    retval[i] = False
        return retval

    def loadCentroidListFile(self, centroidListFileName):
        """
        Load an astrometry (centroid list) file into the 
        current obs file instance.
        """
        scratchDir = os.getenv('INTERM_PATH', '/')
        centroidListPath = os.path.join(scratchDir, 'centroidListFiles')
        fullCentroidListFileName = os.path.join(centroidListPath, centroidListFileName)
        if (not os.path.exists(centroidListFileName)):
            print 'Astrometry centroid list file does not exist: ', centroidListFileName
            return
        self.centroidListFile = tables.openFile(centroidListFileName)
        
    def loadFlatCalFile(self, flatCalFileName):
        """
        loads the flat cal factors from the given file
        """
        scratchDir = os.getenv('INTERM_PATH', '/')
        flatCalPath = os.path.join(scratchDir, 'flatCalSolnFiles')
        fullFlatCalFileName = os.path.join(flatCalPath, flatCalFileName)
        if (not os.path.exists(fullFlatCalFileName)):
            print 'flat cal file does not exist: ', fullFlatCalFileName
            return
        self.flatCalFile = tables.openFile(fullFlatCalFileName, mode='r')
        self.flatWeights = self.flatCalFile.root.flatcal.weights.read()
        self.flatFlags = self.flatCalFile.root.flatcal.flags.read()
        self.flatCalWvlBins = self.flatCalFile.root.flatcal.wavelengthBins.read()
        self.nFlatCalWvlBins = self.flatWeights.shape[2]
        
    def loadFluxCalFile(self, fluxCalFileName):
        """
        loads the flux cal factors from the given file
        """
        scratchDir = os.getenv('INTERM_PATH', '/')
        fluxCalPath = os.path.join(scratchDir, 'fluxCalSolnFiles')
        fullFluxCalFileName = os.path.join(fluxCalPath, fluxCalFileName)
        if (not os.path.exists(fullFluxCalFileName)):
            print 'flux cal file does not exist: ', fullFluxCalFileName
            return
        self.fluxCalFile = tables.openFile(fullFluxCalFileName, mode='r')
        self.fluxWeights = self.fluxCalFile.root.fluxcal.weights.read()
        self.fluxFlags = self.fluxCalFile.root.fluxcal.flags.read()
        self.fluxCalWvlBins = self.fluxCalFile.root.fluxcal.wavelengthBins.read()
        self.nFluxCalWvlBins = self.nFlatCalWvlBins

    def loadHotPixCalFile(self, hotPixCalFileName, switchOnMask=True):
        """
        Load a hot pixel time mask from the given file, in a similar way to
        loadWvlCalFile, loadFlatCalFile, etc. Switches on hot pixel
        masking by default.
        Set switchOnMask=False to prevent switching on hot pixel masking.
        """
        import hotpix.hotPixels as hotPixels    #Here instead of at top to prevent circular import problems.

        scratchDir = os.getenv('INTERM_PATH', '/')
        hotPixCalPath = os.path.join(scratchDir, 'hotPixCalFiles')
        fullHotPixCalFileName = os.path.join(hotPixCalPath, hotPixCalFileName)
        if (not os.path.exists(fullHotPixCalFileName)):
            print 'Hot pixel cal file does not exist: ', fullHotPixCalFileName
            return

        self.hotPixFile = tables.openFile(fullHotPixCalFileName)
        self.hotPixTimeMask = hotPixels.readHotPixels(self.hotPixFile)
        
        if (os.path.basename(self.hotPixTimeMask['obsFileName'])
            != os.path.basename(self.fileName)):
            warnings.warn('Mismatch between hot pixel time mask file and obs file. Not loading/applying mask!')
            self.hotPixTimeMask = None
        else:
            if switchOnMask: self.switchOnHotPixTimeMask()

        #print "end of loadHotPixCalFile.  keys=",self.hotPixTimeMask.keys()
        #print "intervals.shape=",self.hotPixTimeMask['intervals'].shape
        #print "one interval"
        #for iRow in range(self.nRow):
            #for iCol in range(self.nCol):
                #print "iRow=",iRow," iCol=",iCol
                #for interval in self.hotPixTimeMask['intervals'][iRow][iCol]:
                    #print "   interval=",interval

    def loadCosmicMask(self, cosmicMaskFileName=None, switchOnCosmicMask=True):
        self.cosmicMask = ObsFile.readCosmicIntervalFromFile(cosmicMaskFileName)
        if switchOnCosmicMask: self.switchOnCosmicTimeMask()

    def setCosmicMask(self, cosmicMask, switchOnCosmicMask=True):
        self.cosmicMask = cosmicMask
        if switchOnCosmicMask: self.switchOnCosmicTimeMask()

    def loadTimeAdjustmentFile(self,timeAdjustFileName,verbose=False):
        """
        loads obsfile specific adjustments to add to all timestamps read
        adjustments are read from timeAdjustFileName
        it is suggested to pass timeAdjustFileName=FileName(run=run).timeAdjustments()
        """

        self.timeAdjustFile = tables.openFile(timeAdjustFileName)
        self.firmwareDelay = self.timeAdjustFile.root.timeAdjust.firmwareDelay.read()[0]['firmwareDelay']
        roachDelayTable = self.timeAdjustFile.root.timeAdjust.roachDelays
        try:
            self.roachDelays = roachDelayTable.readWhere('obsFileName == "%s"'%self.fileName)[0]['roachDelays']
        except IndexError:
            self.timeAdjustFile.close()
            self.timeAdjustFile=None
            del self.firmwareDelay
            if verbose:
                print 'Unable to load time adjustment for '+self.fileName
            raise
                
    def loadWvlCalFile(self, wvlCalFileName):
        """
        loads the wavelength cal coefficients from a given file
        """
        scratchDir = '/Scratch'
        wvlDir = os.path.join(scratchDir, 'waveCalSolnFiles')
        fullWvlCalFileName = os.path.join(wvlDir, wvlCalFileName)
        if (not os.path.exists(fullWvlCalFileName)):
            print 'wavelength cal file does not exist: ', fullWvlCalFileName
            return
        self.wvlCalFile = tables.openFile(fullWvlCalFileName, mode='r')
        wvlCalData = self.wvlCalFile.root.wavecal.calsoln
        self.wvlCalTable = np.zeros([self.nRow, self.nCol, ObsFile.nCalCoeffs])
        self.wvlErrorTable = np.zeros([self.nRow, self.nCol])
        self.wvlFlagTable = np.zeros([self.nRow, self.nCol])
        self.wvlRangeTable = np.zeros([self.nRow, self.nCol, 2])
        for calPixel in wvlCalData:
            self.wvlFlagTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['wave_flag']
            self.wvlErrorTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['sigma']
            if calPixel['wave_flag'] == 0:
                self.wvlCalTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['polyfit']
                self.wvlRangeTable[calPixel['pixelrow']][calPixel['pixelcol']] = calPixel['solnrange']

    @staticmethod
    def makeWvlBins(energyBinWidth=.1, wvlStart=3000, wvlStop=13000):
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
        energyStop = ObsFile.h * ObsFile.c * ObsFile.angstromPerMeter / wvlStart
        energyStart = ObsFile.h * ObsFile.c * ObsFile.angstromPerMeter / wvlStop
        nWvlBins = int((energyStop - energyStart) / energyBinWidth)
        #Construct energy bin edges
        energyBins = np.linspace(energyStart, energyStop, nWvlBins + 1)
        #Convert back to wavelength and reverse the order to get increasing wavelengths
        wvlBinEdges = np.array(ObsFile.h * ObsFile.c * ObsFile.angstromPerMeter / energyBins)
        wvlBinEdges = wvlBinEdges[::-1]
        return wvlBinEdges

    def parsePhotonPackets(self, packets, inter=interval(),
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
            packetsAll = np.array(packets, dtype='uint64') #64 bit photon packet
            timestampsAll = np.bitwise_and(packets, self.timestampMask)

            if doParabolaFitPeaks:
                parabolaFitPeaksAll = np.bitwise_and(\
                    np.right_shift(packets, self.nBitsAfterParabolaPeak), \
                        self.pulseMask)
            else:
                parabolaFitPeaksAll = np.arange(0)

            if doBaselines:
                baselinesAll = np.bitwise_and(\
                    np.right_shift(packets, self.nBitsAfterBaseline), \
                        self.pulseMask)
            else:
                baselinesAll = np.arange(0)

            if inter == interval() or len(timestampsAll) == 0:
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
        return timestamps, parabolaFitPeaks, baselines

    def plotApertureSpectrum(self, pixelRow, pixelCol, radius1, radius2, weighted=False, fluxWeighted=False, lowCut=3000, highCut=7000, firstSec=0,integrationTime=-1):
    	summed_array, bin_edges = self.getApertureSpectrum(pixelCol=pixelCol, pixelRow=pixelRow, radius1=radius1, radius2=radius2, weighted=weighted, fluxWeighted=fluxWeighted, lowCut=lowCut, highCut=highCut, firstSec=firstSec,integrationTime=integrationTime)
    	fig = plt.figure()
    	ax = fig.add_subplot(111)
    	ax.plot(bin_edges[12:-2], summed_array[12:-1])
    	plt.xlabel('Wavelength ($\AA$)')
    	plt.ylabel('Counts')
    	plt.show()

    def plotPixelSpectra(self, pixelRow, pixelCol, firstSec=0, integrationTime= -1,
                         weighted=False, fluxWeighted=False):
        """
        plots the wavelength calibrated spectrum of a given pixel integrated over a given time
        if integrationTime is -1, All time after firstSec is used.  
        if weighted is True, flat cal weights are applied
        """
        spectrum = (self.getPixelSpectrum(pixelRow, pixelCol, firstSec, integrationTime,
                    weighted=weighted, fluxWeighted=fluxWeighted))['spectrum']
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.flatCalWvlBins[0:-1], spectrum, label='spectrum for pixel[%d][%d]' % (pixelRow, pixelCol))
        plt.show()

    def setWvlCutoffs(self, wvlLowerLimit=3000, wvlUpperLimit=8000):
        """
        Sets wavelength cutoffs so that if convertToWvl(excludeBad=True) or getPixelWvlList(excludeBad=True) is called
        wavelengths outside these limits are excluded.  To remove limits
        set wvlLowerLimit and/or wvlUpperLimit to None.  To use the wavecal limits
        for each individual pixel, set wvlLowerLimit and/or wvlUpperLimit to -1 
        NB - changed defaults for lower/upper limits to None (from 3000 and 8000). JvE 2/22/13
        """
        self.wvlLowerLimit = wvlLowerLimit
        self.wvlUpperLimit = wvlUpperLimit
        
    def switchOffHotPixTimeMask(self):
        """
        Switch off hot pixel time masking - bad pixel times will no longer be
        removed (although the mask remains 'loaded' in ObsFile instance).
        """
        self.hotPixIsApplied = False

    def switchOnHotPixTimeMask(self):
        """
        Switch on hot pixel time masking. Subsequent calls to getPixelCountImage
        etc. will have bad pixel times removed.
        """
        if self.hotPixTimeMask is None:
            raise RuntimeError, 'No hot pixel file loaded'
        self.hotPixIsApplied = True

    def switchOffCosmicTimeMask(self):
        """
        Switch off hot pixel time masking - bad pixel times will no longer be
        removed (although the mask remains 'loaded' in ObsFile instance).
        """
        self.cosmicMaskIsApplied = False

    def switchOnCosmicTimeMask(self):
        """
        Switch on cosmic time masking. Subsequent calls to getPixelCountImage
        etc. will have cosmic times removed.
        """
        if self.cosmicMask is None:
            raise RuntimeError, 'No cosmic mask file loaded'
        self.cosmicMaskIsApplied = True

    @staticmethod
    def writeCosmicIntervalToFile(intervals, ticksPerSec, fileName,
                                  beginTime, endTime, stride,
                                  threshold, nSigma, populationMax):
        h5f = tables.openFile(fileName, 'w')

        headerGroup = h5f.createGroup("/", 'Header', 'Header')
        headerTable = h5f.createTable(headerGroup,'Header',
                                      cosmicHeaderDescription, 'Header')
        header = headerTable.row
        header['ticksPerSec'] = ticksPerSec
        header['beginTime'] = beginTime
        header['endTime'] = endTime
        header['stride'] = stride
        header['threshold'] = threshold
        header['nSigma'] = nSigma
        header['populationMax'] = populationMax 
        header.append()
        headerTable.flush()
        headerTable.close()
        tbl = h5f.createTable("/", "cosmicMaskData", TimeMask.TimeMask, 
                              "Cosmic Mask")
        for interval in intervals:
            row = tbl.row
            row['tBegin'] = max(0,int(np.round(interval[0]*ticksPerSec)))
            row['tEnd'] = int(np.round(interval[1]*ticksPerSec))
            row['reason'] = TimeMask.timeMaskReason["cosmic"]
            row.append()
            tbl.flush()
        tbl.close()
        h5f.close()

    @staticmethod
    def readCosmicIntervalFromFile(fileName):
        fid = tables.openFile(fileName, mode='r')
        headerInfo = fid.getNode("/Header","Header")[0]
        ticksPerSec = headerInfo['ticksPerSec']
        table = fid.getNode("/cosmicMaskData")
        enum = table.getEnum('reason')

        retval = interval()
        for i in range(table.nrows):
            temp = (interval[table[i]['tBegin'],table[i]['tEnd']])/ticksPerSec
            retval = retval | temp

        fid.close()
        return retval


    def writePhotonList(self,*nkwargs,**kwargs): #filename=None, firstSec=0, integrationTime=-1):                       
        """
        Write out the photon list for this obs file.
        See photonlist/photlist.py for input parameters and outputs.
        Shifted over to photonlist/, May 10 2013, JvE. All under construction at the moment.
        """        
        import photonlist.photlist      #Here instead of at top to avoid circular imports
        photonlist.photlist.writePhotonList(self,*nkwargs,**kwargs)
        
        
#        writes out the photon list for this obs file at $INTERM_PATH/photonListFileName
#        currently cuts out photons outside the valid wavelength ranges from the wavecal
#       
#        Currently being updated... JvE 4/26/2013.
#        This version should automatically reject time-masked photons assuming a hot pixel mask is
#        loaded and 'switched on'.
#        
#        INPUTS:
#            filename - string, optionally use to specify non-default output file name
#                       for photon list. If not supplied, default name/path is determined
#                       using original obs. file name and standard directory paths (as per
#                       util.FileName). Added 4/29/2013, JvE.
#            firstSec - Start time within the obs. file from which to begin the
#                       photon list (in seconds, from the beginning of the obs. file).
#            integrationTime - Length of exposure time to extract (in sec, starting from
#                       firstSec). -1 to extract to end of obs. file.
#        
#        """
#        
#        if self.flatCalFile is None: raise RuntimeError, "No flat cal. file loaded"
#        if self.fluxCalFile is None: raise RuntimeError, "No flux cal. file loaded"
#        if self.wvlCalFile is None: raise RuntimeError, "No wavelength cal. file loaded"
#        if self.hotPixFile is None: raise RuntimeError, "No hot pixel file loaded"
#        if self.file is None: raise RuntimeError, "No obs file loaded...?"
#        
#        plFile = self.createEmptyPhotonListFile(filename)
#        #try:
#        plTable = plFile.root.photons.photons
#                
#        try:
#            plFile.copyNode(self.flatCalFile.root.flatcal, newparent=plFile.root, newname='flatcal', recursive=True)
#            plFile.copyNode(self.fluxCalFile.root.fluxcal, newparent=plFile.root, newname='fluxcal', recursive=True)
#            plFile.copyNode(self.wvlCalFile.root.wavecal, newparent=plFile.root, newname='wavecal', recursive=True)
#            plFile.copyNode(self.hotPixFile.root, newparent=plFile.root, newname='timemask', recursive=True)
#            plFile.copyNode(self.file.root.beammap, newparent=plFile.root, newname='beammap', recursive=True)
#            plFile.copyNode(self.file.root.header, newparent=plFile.root, recursive=True)
#        except:
#            plFile.flush()
#            plFile.close()
#            raise
#        
#        plFile.flush()
#
#        fluxWeights = self.fluxWeights      #Flux weights are independent of pixel location.
#        #Extend flux weight/flag arrays as for flat weight/flags.
#        fluxWeights = np.hstack((fluxWeights[0],fluxWeights,fluxWeights[-1]))
#        fluxFlags = np.hstack((pipelineFlags.fluxCal['belowWaveCalRange'], 
#                               self.fluxFlags, 
#                               pipelineFlags.fluxCal['aboveWaveCalRange']))
#
#        for iRow in xrange(self.nRow):
#            for iCol in xrange(self.nCol):
#                flag = self.wvlFlagTable[iRow, iCol]
#                if flag == 0:#only write photons in good pixels  ***NEED TO UPDATE TO USE DICTIONARY***
#                    energyError = self.wvlErrorTable[iRow, iCol] #Note wvlErrorTable is in eV !! Assume constant across all wavelengths. Not the best approximation, but a start....
#                    flatWeights = self.flatWeights[iRow, iCol]
#                    #Extend flat weight and flag arrays at beginning and end to include out-of-wavelength-calibration-range photons.
#                    flatWeights = np.hstack((flatWeights[0],flatWeights,flatWeights[-1]))
#                    flatFlags = np.hstack((pipelineFlags.flatCal['belowWaveCalRange'],
#                                           self.flatFlags[iRow, iCol],
#                                           pipelineFlags.flatCal['aboveWaveCalRange']))
#                    
#                    
#                    #wvlRange = self.wvlRangeTable[iRow, iCol]
#
#                    #---------- Replace with call to getPixelWvlList -----------
#                    #go through the list of seconds in a pixel dataset
#                    #for iSec, secData in enumerate(self.getPixel(iRow, iCol)):
#                        
#                    #timestamps, parabolaPeaks, baselines = self.parsePhotonPackets(secData)
#                    #timestamps = iSec + self.tickDuration * timestamps
#                 
#                    #pulseHeights = np.array(parabolaPeaks, dtype='double') - np.array(baselines, dtype='double')
#                    #wavelengths = self.convertToWvl(pulseHeights, iRow, iCol, excludeBad=False)
#                    #------------------------------------------------------------
#
#                    x = self.getPixelWvlList(iRow,iCol,excludeBad=False,dither=True,firstSec=firstSec,
#                                             integrationTime=integrationTime)
#                    timestamps, wavelengths = x['timestamps'], x['wavelengths']     #Wavelengths in Angstroms
#                    
#                    #Convert errors in eV to errors in Angstroms (see notebook, May 7 2013)
#                    wvlErrors = ((( (energyError*units.eV) * (wavelengths*units.Angstrom)**2 ) /
#                                    (constants.h*constants.c) )
#                                 .to(units.Angstrom).value)
#                        
#                    #Calculate what wavelength bin each photon falls into to see which flat cal factor should be applied
#                    if len(wavelengths) > 0:
#                        flatBinIndices = np.digitize(wavelengths, self.flatCalWvlBins)      #- 1 - 
#                    else:
#                        flatBinIndices = np.array([])
#
#                    #Calculate which wavelength bin each photon falls into for the flux cal weight factors.
#                    if len(wavelengths) > 0:
#                        fluxBinIndices = np.digitize(wavelengths, self.fluxCalWvlBins)
#                    else:
#                        fluxBinIndices = np.array([])
#
#                    for iPhoton in xrange(len(timestamps)):
#                        #if wavelengths[iPhoton] > wvlRange[0] and wavelengths[iPhoton] < wvlRange[1] and binIndices[iPhoton] >= 0 and binIndices[iPhoton] < len(flatWeights):
#                        #create a new row for the photon list
#                        newRow = plTable.row
#                        newRow['Xpix'] = iCol
#                        newRow['Ypix'] = iRow
#                        newRow['ArrivalTime'] = timestamps[iPhoton]
#                        newRow['Wavelength'] = wavelengths[iPhoton]
#                        newRow['WaveError'] = wvlErrors[iPhoton]
#                        newRow['FlatFlag'] = flatFlags[flatBinIndices[iPhoton]]
#                        newRow['FlatWeight'] = flatWeights[flatBinIndices[iPhoton]]
#                        newRow['FluxFlag'] = fluxFlags[fluxBinIndices[iPhoton]]
#                        newRow['FluxWeight'] = fluxWeights[fluxBinIndices[iPhoton]]
#                        newRow.append()
#        #finally:
#        plTable.flush()
#        plFile.close()
#


            

def calculateSlices_old(inter, timestamps):
    """
    return a list of strings, with format i0:i1 for a python array slice
    inter is the interval of values in timestamps to mask out.
    The resulting list of strings indicate elements that are not masked out
    
    Renamed to calculateSlices_old and deprecated - JvE 3/8/2013
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
        slce += str(i + 1) # end the last range if it is included
        slices.append(slce)
    return slices

def calculateSlices(inter, timestamps):
    '''
    Hopefully a quicker version of  the original calculateSlices. JvE 3/8/2013
    
    Returns a list of strings, with format i0:i1 for a python array slice
    inter is the interval of values in timestamps to mask out.
    The resulting list of strings indicate elements that are not masked out
    
    inter must be a single pyinterval 'interval' object (can be multi-component)
    timestamps is a 1D array of timestamps (MUST be an *ordered* array).
    
    If inter is a multi-component interval, the components must be unioned and sorted
    (which is the default behaviour when intervals are defined, and is probably
    always the case, so shouldn't be a problem).
    '''
    timerange = interval([timestamps[0],timestamps[-1]])
    slices = []
    slce = "0:"     #Start at the beginning of the timestamps array....
    imax = 0        #Will prevent error if inter is an empty interval
    for eachComponent in inter.components:
        #Check if eachComponent of the interval overlaps the timerange of the 
        #timestamps - if not, skip to the next component.

        if eachComponent & timerange == interval(): continue
        #[
        #Possibly a bit faster to do this and avoid interval package, but not fully tested:
        #if eachComponent[0][1] < timestamps[0] or eachComponent[0][0] > timestamps[-1]: continue
        #]

        imin = np.searchsorted(timestamps, eachComponent[0][0], side='left') #Find nearest timestamp to lower bound
        imax = np.searchsorted(timestamps, eachComponent[0][1], side='right') #Nearest timestamp to upper bound
        #As long as we're not about to create a wasteful '0:0' slice, go ahead 
        #and finish the new slice and append it to the list
        if imin != 0:
            slce += str(imin)
            slices.append(slce)
        slce = str(imax)+":"
    #Finish the last slice at the end of the timestamps array if we're not already there:
    if imax != len(timestamps):
        slce += str(len(timestamps))
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
        iPtNew = iPt + s1 - s0        
        retval[iPt:iPtNew] = array[s0:s1]
        iPt = iPtNew
    return retval

class cosmicHeaderDescription(tables.IsDescription):
    ticksPerSec = tables.Float64Col() # number of ticks per second
    beginTime = tables.Float64Col()   # begin time used to find cosmics (seconds)
    endTime = tables.Float64Col()     # end time used to find cosmics (seconds)
    stride = tables.Int32Col()
    threshold = tables.Float64Col()
    nSigma = tables.Int32Col()
    populationMax = tables.Int32Col()


#Temporary test
if __name__ == "__main__":
    obs = ObsFile(FileName(run='PAL2012', date='20121210', tstamp='20121211-051650').obs())
    obs.loadWvlCalFile(FileName(run='PAL2012',date='20121210',tstamp='20121211-052230').calSoln())
    obs.loadFlatCalFile(FileName(obsFile=obs).flatSoln())
    beforeImg = obs.getPixelCountImage(weighted=False,fluxWeighted=False,scaleByEffInt=True)
