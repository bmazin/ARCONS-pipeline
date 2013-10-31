'''
Author: Julian van Eyken                Date: May 7 2013
For handling calibrated output photon lists.
'''

import numpy as np
import os.path
import warnings
import tables
from astropy import constants, units
from headers import pipelineFlags
import ArconsHeaders
from util import utils
from FileName import FileName

xyPackMult = 100  #Multiplier factor for turning row,col coordinate pairs into single integers for storage in HDF file (see xyPack).


class PhotList(object):
    '''
    Class to hold a fully calibrated photon list. Beginnings of a full
    library of operations similar to the ObsFile class.
    '''
    
    def __init__(self, fileName):
        '''
        Initialise by loading a photon list file.
        '''
        self.file = None
        self.fileName = None
        self.fullFileName = None
        self.nRow = None
        self.nCol = None
        self.startTime = None   #To be implemented!
        self.photTable = None    #To hold the photon list table node (just a shortcut)
        self.loadFile(fileName)
        
    
    def __del__(self):
        '''
        Clean up when object instance is deleted or goes out of scope.
        '''
        self.file.close()
        
    
    #def __iter__(self):
    #    '''
    #    Iterate over individual photons in the list.
    #    To be implemented....
    #    '''
    #
    
    def loadFile(self,fileName):
        '''
        Open up the .h5 photon list file for reading.
        '''
        
        #if (os.path.isabs(fileName)):
        self.fileName = os.path.basename(fileName)
        self.fullFileName = os.path.abspath(fileName)
        #else:
        #    # make the full file name by joining the input name 
        #    # to the MKID_DATA_DIR (or . if the environment variable 
        #    # is not defined)
        #    dataDir = os.getenv('MKID_DATA_DIR', '/')
        #    self.fullFileName = os.path.join(dataDir, self.fileName)

        if (not os.path.exists(self.fullFileName)):
            msg='file does not exist: %s'%self.fullFileName
            #if verbose:
            #    print msg
            raise Exception(msg)

        #open the hdf5 file
        self.file = tables.openFile(self.fullFileName, mode='a')

        #Figure out the number of rows and columns in the detector array.
        self.nRow, self.nCol = self.file.root.beammap.beamimage.shape

        #Set the photon-list node shortcut
        self.photTable = self.file.root.photons.photons

        #get the header
        self.header = self.file.root.header.header
        self.titles = self.header.colnames
        try:
            self.info = self.header[0] #header is a table with one row
        except IndexError as inst:
            if verbose:
                print 'Can\'t read header for ',self.fullFileName
            raise inst
        
        #Can Parse the hot pixel timemasks here if necessary, but may be better
        #left so that it's only parsed when/if needed.
        
    
    def loadData(self):
        self.data = self.photTable.read()

    def getPhotonsForPixel(self):
        '''
        Space holder for now
        '''
        pass
    
    
    def getImageDet(self,firstSec=-np.Inf,integrationTime=-1,wvlMin=-np.Inf,
                 wvlMax=np.Inf):
        '''
        Return a 2D image consisting of photon counts in each pixel, in the detector pixel
        coordinate frame.
        '''
        if integrationTime==-1:
            lastSec=np.Inf
        else:  
            lastSec = firstSec+integrationTime
        
        #Initialise an empty image
        image = np.empty((self.nRow,self.nCol),dtype=float)
        image.fill(np.NaN)
        
        #Step through pixels and fill em up.
        for iRow in range(self.nRow):    #self.nCol):
            print 'Reading pixel row ', iRow
            for iCol in range(self.nCol):
                rowcol = xyPack(iRow,iCol)
                #image[iRow,iCol] = len(self.photTable.getWhereList('(Xpix == iCol) & (Ypix == iRow) & (ArrivalTime >= firstSec) & (ArrivalTime < lastSec)')) # & (Wavelength >= wvlMin) & (Wavelength < wvlMax)' ))
                #image[iRow,iCol] = len(self.photTable.getWhereList(('XYpix==rowcol'))) # & (self.photTable.cols.Ypix == iRow) & (self.photTable.cols.ArrivalTime >= firstSec)
                                                                    #& (self.photTable.cols.ArrivalTime < lastSec) & (self.photTable.cols.Wavelength >= wvlMin) & (self.photTable.cols.Wavelength < wvlMax) ))

                photons = (self.photTable.readWhere('XYpix==rowcol'))
                #image[iRow,iCol] = len(photons.where('ArrivalTime >= firstSec) & ArrivalTime < lastSec)'))
                
                image[iRow,iCol] = np.sum((photons['ArrivalTime'] >= firstSec) & (photons['ArrivalTime'] < lastSec)
                                          & (photons['Wavelength'] >= wvlMin) & (photons['Wavelength'] < wvlMax))
                #image[iRow,iCol] = len(photons)
                #image[iRow,iCol] = len(self.photTable.getWhereList('(ArrivalTime >= firstSec) & (ArrivalTime < lastSec) & (Ypix == iRow) & (Xpix == iCol)')) # & (Wavelength >= wvlMin) & (Wavelength < wvlMax)' ))                
                #image[iRow,iCol] = len(self.photTable.getWhereList('(Ypix == iRow)' ))
                
        #That should be it....
        return image

    def getImageSky(self,firstSec=-np.Inf,integrationTime=-1,wvlMin=-np.Inf,wvlMax=np.Inf):
        '''
        Get a 2D image in the sky frame - i.e, in RA, dec coordinate space, derotated and stacked
        if necessary.
        '''
        
        
def createEmptyPhotonListFile(obsFile,fileName=None):
     """
     creates a photonList h5 file 
     using header in headers.ArconsHeaders
     
     INPUTS:
         fileName - string, name of file to write to. If not supplied, default is used
                    based on name of original obs. file and standard directories etc.
                    (see usil.FileName). Added 4/29/2013, JvE
     """
     
     if fileName is None:    
         fileTimestamp = obsFile.fileName.split('_')[1].split('.')[0]
         fileDate = os.path.basename(os.path.dirname(obsFile.fullFileName))
         run = os.path.basename(os.path.dirname(os.path.dirname(obsFile.fullFileName)))
         fn = FileName(run=run, date=fileDate, tstamp=fileTimestamp)
         fullPhotonListFileName = fn.photonList()
     else:
         fullPhotonListFileName = fileName
     if (os.path.exists(fullPhotonListFileName)):
         if utils.confirm('Photon list file  %s exists. Overwrite?' % fullPhotonListFileName, defaultResponse=False) == False:
             warnings.warn('No photon list file created')
             return
     zlibFilter = tables.Filters(complevel=1, complib='zlib', fletcher32=False)
     try:
         plFile = tables.openFile(fullPhotonListFileName, mode='w')
         plGroup = plFile.createGroup('/', 'photons', 'Group containing photon list')
         plTable = plFile.createTable(plGroup, 'photons', ArconsHeaders.CrabList, 'Photon List Data', 
                                      filters=zlibFilter, 
                                      expectedrows=300000)  #Temporary fudge to see if it helps!
     except:
         plFile.close()
         raise
     return plFile




def writePhotonList(obsFile, filename=None, firstSec=0, integrationTime=-1,rowList=None,colList=None):
    """
    writes out the photon list for this obs file at $INTERM_PATH/photonListFileName
    currently cuts out photons outside the valid wavelength ranges from the wavecal
   
    Currently being updated... JvE 4/26/2013.
    This version should automatically reject time-masked photons assuming a hot pixel mask is
    loaded and 'switched on'.
    
    INPUTS:
        obsFile - the obsFile instance from which the list is to be created.
        filename - string, optionally use to specify non-default output file name
                   for photon list. If not supplied, default name/path is determined
                   using original obs. file name and standard directory paths (as per
                   util.FileName). Added 4/29/2013, JvE.
        firstSec - Start time within the obs. file from which to begin the
                   photon list (in seconds, from the beginning of the obs. file).
        integrationTime - Length of exposure time to extract (in sec, starting from
                   firstSec). -1 to extract to end of obs. file.
    
    """
    
    if obsFile.flatCalFile is None: raise RuntimeError, "No flat cal. file loaded"
    if obsFile.fluxCalFile is None: raise RuntimeError, "No flux cal. file loaded"
    if obsFile.wvlCalFile is None: raise RuntimeError, "No wavelength cal. file loaded"
    if obsFile.hotPixFile is None: raise RuntimeError, "No hot pixel file loaded"
    if obsFile.file is None: raise RuntimeError, "No obs file loaded...?"
    
    #print 'Initialising empty photon list file '+filename
    plFile = createEmptyPhotonListFile(obsFile,filename)
    #try:
    plTable = plFile.root.photons.photons
            
    try:
        plFile.copyNode(obsFile.flatCalFile.root.flatcal, newparent=plFile.root, newname='flatcal', recursive=True)
        plFile.copyNode(obsFile.fluxCalFile.root.fluxcal, newparent=plFile.root, newname='fluxcal', recursive=True)
        plFile.copyNode(obsFile.wvlCalFile.root.wavecal, newparent=plFile.root, newname='wavecal', recursive=True)
        plFile.copyNode(obsFile.hotPixFile.root, newparent=plFile.root, newname='timemask', recursive=True)
        plFile.copyNode(obsFile.file.root.beammap, newparent=plFile.root, newname='beammap', recursive=True)
        plFile.copyNode(obsFile.file.root.header, newparent=plFile.root, recursive=True)
    except:
        plFile.flush()
        plFile.close()
        raise
    
    plFile.flush()

    fluxWeights = obsFile.fluxWeights      #Flux weights are independent of pixel location.
    #Extend flux weight/flag arrays as for flat weight/flags.
    fluxWeights = np.hstack((fluxWeights[0],fluxWeights,fluxWeights[-1]))
    fluxFlags = np.hstack((pipelineFlags.fluxCal['belowWaveCalRange'], 
                           obsFile.fluxFlags, 
                           pipelineFlags.fluxCal['aboveWaveCalRange']))

    for iRow in xrange(obsFile.nRow):
        print 'Writing pixel row ', iRow, '...'
        for iCol in xrange(obsFile.nCol):
            flag = obsFile.wvlFlagTable[iRow, iCol]
            if flag == 0:#only write photons in good pixels  ***NEED TO UPDATE TO USE DICTIONARY***
                    
                if rowList == None or ((iRow,iCol) in zip(rowList,colList)):
                    energyError = obsFile.wvlErrorTable[iRow, iCol] #Note wvlErrorTable is in eV !! Assume constant across all wavelengths. Not the best approximation, but a start....
                    flatWeights = obsFile.flatWeights[iRow, iCol]
                    #Extend flat weight and flag arrays at beginning and end to include out-of-wavelength-calibration-range photons.
                    flatWeights = np.hstack((flatWeights[0],flatWeights,flatWeights[-1]))
                    flatFlags = np.hstack((pipelineFlags.flatCal['belowWaveCalRange'],
                                           obsFile.flatFlags[iRow, iCol],
                                           pipelineFlags.flatCal['aboveWaveCalRange']))
                    
                    
                    wvlRange = obsFile.wvlRangeTable[iRow, iCol]
                    if wvlRange[1] >= 11000:

                        #---------- Replace with call to getPixelWvlList -----------
                        #go through the list of seconds in a pixel dataset
                        #for iSec, secData in enumerate(obsFile.getPixel(iRow, iCol)):
                            
                        #timestamps, parabolaPeaks, baselines = obsFile.parsePhotonPackets(secData)
                        #timestamps = iSec + obsFile.tickDuration * timestamps
                     
                        #pulseHeights = np.array(parabolaPeaks, dtype='double') - np.array(baselines, dtype='double')
                        #wavelengths = obsFile.convertToWvl(pulseHeights, iRow, iCol, excludeBad=False)
                        #------------------------------------------------------------

                        x = obsFile.getPixelWvlList(iRow,iCol,excludeBad=False,dither=True,firstSec=firstSec,
                                                 integrationTime=integrationTime)
                        timestamps, wavelengths = x['timestamps'], x['wavelengths']     #Wavelengths in Angstroms
                        
                        #Convert errors in eV to errors in Angstroms (see notebook, May 7 2013)
                        wvlErrors = ((( (energyError*units.eV) * (wavelengths*units.Angstrom)**2 ) /
                                        (constants.h*constants.c) )
                                     .to(units.Angstrom).value)
                            
                        #Calculate what wavelength bin each photon falls into to see which flat cal factor should be applied
                        if len(wavelengths) > 0:
                            flatBinIndices = np.digitize(wavelengths, obsFile.flatCalWvlBins)      #- 1 - 
                        else:
                            flatBinIndices = np.array([])

                        #Calculate which wavelength bin each photon falls into for the flux cal weight factors.
                        if len(wavelengths) > 0:
                            fluxBinIndices = np.digitize(wavelengths, obsFile.fluxCalWvlBins)
                        else:
                            fluxBinIndices = np.array([])

                        outOfWaveCalRange = np.logical_or(wavelengths < wvlRange[0],wavelengths > wvlRange[1])
                        for iPhoton in xrange(len(timestamps)):
                            #if wavelengths[iPhoton] > wvlRange[0] and wavelengths[iPhoton] < wvlRange[1]:
                                #create a new row for the photon list
                                newRow = plTable.row
                                newRow['xPix'] = iCol
                                newRow['yPix'] = iRow
                                newRow['xyPix'] = xyPack(iRow,iCol)
                                newRow['arrivalTime'] = timestamps[iPhoton]
                                newRow['wavelength'] = wavelengths[iPhoton]
                                newRow['waveError'] = outOfWaveCalRange[iPhoton]
                                newRow['flatFlag'] = flatFlags[flatBinIndices[iPhoton]]
                                newRow['flatWeight'] = flatWeights[flatBinIndices[iPhoton]]
                                newRow['fluxFlag'] = fluxFlags[fluxBinIndices[iPhoton]]
                                newRow['fluxWeight'] = fluxWeights[fluxBinIndices[iPhoton]]
                                newRow.append()
    #finally:
    plTable.flush()
    plFile.close()
    print 'Done.'

def indexPhotonList(photFile):
    '''
    To run pytables indexing on an existing photon list file. Should speed up
    querying a lot...
    
    INPUTS:
        photFile - either a string containing the filename of the .h5 photon list file
                   to update, or a pytables file instance containing the photon list 
                   to index.
                   
    OUTPUTS:
        If photFile is a filename, the file is opened, updated in place, and closed.
        If photFile is a file instance, the file is indexed, but not closed.
    '''
    
    colsToIndex = ['XYpix', 'Xpix', 'Ypix', 'ArrivalTime', 'Wavelength']
    
    if type(photFile) is str:
        photListFile = tables.openFile(photFile, mode='r+')
    elif type(photFile) is tables.file.File:
        photListFile = photFile
    else:
        raise ValueError, 'String or pytables file instance required'
    
    photList = photListFile.root.photons.photons
    
    try:
    #Loop through the columns to index, and index them...
        for eachColName in colsToIndex:
            if photList.colinstances[eachColName].is_indexed:
                warnings.warn('Column is already indexed: '+eachColName+' - skipping.')
            else:
                print 'Indexing column: '+eachColName
                dummy = photList.colinstances[eachColName].createCSIndex()
        
    finally:
        photList.flush()
        if type(photFile is str):
            photListFile.close()
 
 
def xyPack(row,col):
    '''
    Pack an integer pair of pixel coordinates, x, y, into a single integer,
    stored in the 'XYpix' column of the photon list. Used for the purposes
    of access efficiency.
    '''
    return xyPackMult * col + row
    
def xyUnpack(rowCol):
    '''
    Inverse of xyPack. Returns an integer tuple, (row,col).
    '''
    return rowCol % xyPackMult, rowCol // xyPackMult
    
    
