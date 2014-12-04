'''
Author: Julian van Eyken                Date: May 7 2013
For handling calibrated output photon lists.
'''

import time
import numpy as np
import os.path
import warnings
import tables
from astropy import constants, units
import matplotlib.pylab as mpl
from headers import pipelineFlags, ArconsHeaders
from util import utils
from astrometry import CalculateRaDec as crd
from util.FileName import FileName
from hotpix import hotPixels as hp
from beammap import remapPixels


xyPackMult = 100  #Multiplier factor for turning row,col coordinate pairs into single integers for storage in HDF file (see xyPack).


class PhotList(object):
    '''
    Class to hold a fully calibrated photon list. Beginnings of a full
    library of operations similar to the ObsFile class.
    '''
    
    def __init__(self, fileName,mode='r'):
        '''
        Initialise by loading a photon list file.
        mode determines whether to open for reading and/or writing
        '''
        self.file = None           #To contain the photon list file object
        self.fileName = None       #Name of the photon list file
        self.fullFileName = None   #Full path name of the same.
        self.nRow = None           #Number of detector rows 
        self.nCol = None           #Number of detector columns
        self.startTime = None      #To be implemented!
        self.photTable = None      #To reference the photon list table node within self.file (just a shortcut)
        self.hotPixTimeMask = None      #To store the hot pixel info dictionary after a hot (or bad) pixel file is loaded.
        self.loadFile(fileName,mode=mode)
    
    def __del__(self):
        '''
        Clean up when object instance is deleted or goes out of scope.
        '''
        self.file.close()
        
    
    def loadFile(self,fileName,doParseHotPixTimeMask=False,mode='r'):
        '''
        Open up the .h5 photon list file for reading.
        INPUTS:
            fileName - path name of .h5 file to load up
            parseHotPixTimeMask - if True, parse the hotpixel time mask data
                        right away, and store the resulting dictionary
                        in self.hotPixTimeMask (as output by hotPixels.readHotPixels()).
            mode - passed to openFile.  'w' mode is not allowed.
        '''
        
        self.fileName = os.path.basename(fileName)
        self.fullFileName = os.path.abspath(fileName)

        if (not os.path.exists(self.fullFileName)):
            msg='file does not exist: %s'%self.fullFileName
            #if verbose:
            #    print msg
            raise Exception(msg)

        if mode == 'w':
            raise Exception('mode is \'w\'. Do not load a PhotonList with mode=\'w\'. This will delete the file')
        #open the hdf5 file
        self.file = tables.openFile(self.fullFileName, mode=mode)

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
        
        if doParseHotPixTimeMask:
            self.parseHotPixTimeMask()
    
    
    
    def parseHotPixTimeMask(self):
        '''Parse in the hot pixel time mask HDF data.'''        
        print 'Parsing pixel time-mask data'
        self.hotPixTimeMask = hp.readHotPixels(self.file)
        print 'Done parsing'

    
    def getBadWvlCalFlags(self):
        '''
        Check the wavelength cal and figure out which pixels have bad (or no)
        wavelength solutions -- should catch unassigned pixels, most dead pixels, etc.
        INPUTS:
            None (except the photon list instance).
        OUTPUTS:
            Returns a 2D integer array representing the wavelength cal flags for
            each pixel.
        '''
        
        flagArr = np.zeros((self.nRow,self.nCol)) 
        wvlTable = self.file.root.wavecal.calsoln
        for row in wvlTable.iterrows():
            flagArr[row['pixelrow'],row['pixelcol']] = row['wave_flag']
        return flagArr
        
    
    def getImageDet(self,firstSec=0,integrationTime=-1,wvlMin=-np.Inf,
                 wvlMax=np.Inf,doWeighted=True,
                 #newMethod=True
                 ):
        '''
        Return a 2D image consisting of photon counts in each pixel, in the detector pixel
        coordinate frame.
        
        INPUTS:
            firstSec - start time for integration within photon list object (in sec)
            integrationTime - time to integrate for in seconds. If -1, integrates to end of file.
            wvlMin,wvlMax - minimum and maximum wavelengths to include in the integration (Angstroms) 
            doWeighted - True to apply flatcal and fluxcal weights.
            
            ----DEPRECATED, NO LONGER FUNCTIONAL - always uses new method (Jul 9 2014)------
            newMethod - True to use numpy histogram for creating image. Runs faster for short integration times,
                        but just a little slower for full images.
            --------------------------------------------------------------------------------
            
        OUPUTS:
            A 2D image array matching the shape/size of the detector.
            
        '''
        
        if integrationTime==-1:
            lastSec=np.Inf
        else:  
            lastSec = firstSec+integrationTime

#         if newMethod:
        #Should be more efficient....
        print 'Searching photons...'
        #ind = self.photTable.getWhereList('(arrivalTime >= firstSec) & (arrivalTime < lastSec) & (wavelength >= wvlMin) & (wavelength < wvlMax)')
        tic = time.clock()
        photons = self.photTable.readWhere('(arrivalTime >= firstSec) & (arrivalTime < lastSec) & (wavelength >= wvlMin) & (wavelength < wvlMax)')
        xPix = photons['xPix']
        yPix = photons['yPix']
        if doWeighted is True:
            photWeights = photons['flatWeight'] * photons['fluxWeight']   #********EXPERIMENTING WITH ADDING FLUX WEIGHT - NOT FULLY TESTED, BUT SEEMS OKAY....********
        else:
            photWeights = None
        print 'Time taken so far (s): ',time.clock()-tic
        print 'Building 2D histogram'
        image, y,x = np.histogram2d(yPix,xPix,bins=[self.nRow,self.nCol],range=[[-1,self.nRow],[-1,self.nCol]],
                                    weights=photWeights)
        print 'Done - total time taken (s): ',time.clock()-tic

#         else:
#             #Initialise an empty image
#             image = np.empty((self.nRow,self.nCol),dtype=float)
#             image.fill(np.NaN)
# 
#             #Step through pixels and fill em up.
#             print 'Looping through pixels and counting up...'
#             tic = time.clock()
#             for iRow in range(self.nRow):    #self.nCol):
#                 print 'Reading pixel row ', iRow
#                 for iCol in range(self.nCol):
#                     rowcol = xyPack(iRow,iCol)
#                     photons = (self.photTable.readWhere('xyPix==rowcol'))
#                     if doWeighted is True:
#                         photWeights = photons['flatWeight'] * photons['fluxWeight']   #********EXPERIMENTING WITH ADDING FLUX WEIGHT - NOT FULLY TESTED, BUT SEEMS OKAY....********
#                     else:
#                         photWeights = np.ones(len(photons))     #Set all weights to one if doWeighted is not True.
#                     inRange = ((photons['arrivalTime'] >= firstSec) & (photons['arrivalTime'] < lastSec)
#                                 & (photons['wavelength'] >= wvlMin) & (photons['wavelength'] < wvlMax))
#                     
#                     image[iRow,iCol] = np.sum(photWeights[inRange])
#             print 'Done - total time taken (s): ',time.clock()-tic

        #That should be it....
        return image


    def getPhotonsForPixel(self):
        '''
        Space holder for now
        '''
        raise NotImplementedError


    def displayImageDet(self, firstSec=0,integrationTime=-1,wvlMin=-np.Inf,
                 wvlMax=np.Inf, normMin=None, normMax=None, showHotPix=False):
        '''
        Display an image built from the photon list in detector coordinate space
        '''
        im = self.getImageDet(firstSec=firstSec, integrationTime=integrationTime, wvlMin=wvlMin, wvlMax=wvlMax)
        utils.plotArray(im, normMin=normMin, normMax=normMax, cbar=True, fignum=None)
        if showHotPix is True:
            if self.hotPixTimeMask is None:
                self.parseHotPixTimeMask()
            #badPix = hp.getHotPixels(self.hotPixTimeMask, firstSec=firstSec, integrationTime=integrationTime)
            badPix = self.hotPixTimeMask.getHotPixels(firstSec=firstSec, integrationTime=integrationTime)
            x = np.arange(self.nCol)
            y = np.arange(self.nRow)
            xx, yy = np.meshgrid(x, y)
            if np.sum(badPix) > 0: mpl.scatter(xx[badPix], yy[badPix], c='y')


    def getImageSky(self,firstSec=-np.Inf,integrationTime=-1,wvlMin=-np.Inf,wvlMax=np.Inf):
        '''
        Get a 2D image in the sky frame - i.e, in RA, dec coordinate space, derotated and stacked
        if necessary.
        '''
        pass #Place holder for now! (Should point to the RADecImage.py code).



        
def createEmptyPhotonListFile(obsFile,fileName=None,photListDescription=ArconsHeaders.PhotonList):
    """
    creates a photonList h5 file 
    using header in headers.ArconsHeaders
    
    INPUTS:
        fileName - string, name of file to write to. If not supplied, default is used
                    based on name of original obs. file and standard directories etc.
                    (see usil.FileName). Added 4/29/2013, JvE
        photListDescription - a pytables description class that lays out the column structure of the photon table
                for a normal photon list, leave as ArconsHeaders.PhotonList.  For analyzing pulsars, use
                ArconsHeaders.PulsarPhotonList, which has a couple extra pulsar-specfic columns
    """
     
    if fileName is None:    
    #fileTimestamp = obsFile.fileName.split('_')[1].split('.')[0]
    #fileDate = os.path.basename(os.path.dirname(obsFile.fullFileName))
    #run = os.path.basename(os.path.dirname(os.path.dirname(obsFile.fullFileName)))
    #fn = FileName(run=run, date=fileDate, tstamp=fileTimestamp)
    #fullPhotonListFileName = fn.photonList()
        fullPhotonListFileName = FileName(obsFile=obsFile).photonList()
    else:
         fullPhotonListFileName = fileName
    if (os.path.exists(fullPhotonListFileName)):
        if utils.confirm('Photon list file %s exists. Overwrite?' % fullPhotonListFileName, defaultResponse=False) == False:
            warnings.warn('No photon list file created')
            return
    zlibFilter = tables.Filters(complevel=1, complib='zlib', fletcher32=False)
    bloscFilter = tables.Filters(complevel=9, complib='blosc', fletcher32=False)    #May be more efficient to use - needs some experimenting with compression level etc.
    plFile = tables.openFile(fullPhotonListFileName, mode='w')
    try:
        plGroup = plFile.createGroup('/', 'photons', 'Group containing photon list')
        plTable = plFile.createTable(plGroup, 'photons', photListDescription, 'Photon List Data', 
                                      filters=bloscFilter) 
                                      #expectedrows=300000)
    except:
        plFile.close()
        raise
    return plFile




def writePhotonList(obsFile, filename=None, firstSec=0, integrationTime=-1, 
                    doIndex=True, pixRemapFileName=None,
                    photListDescription = ArconsHeaders.PhotonList):
    """
    writes out the photon list for this obs file at $MKID_PROC_PATH/photonListFileName
    currently cuts out photons outside the valid wavelength ranges from the wavecal
   
    Currently being updated... JvE 4/26/2013.
    This version should automatically reject time-masked photons assuming a hot pixel mask is
    loaded and 'switched on'.
    
    Added saving of associated time-adjustment file if present - NEEDS TESTING. JvE 7/8/2013
    
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
        astrometryFileName - name of centroid list file, as created by Paul's code, 
                   - see astrometry/centroidCalc.py. If not provided, attempts to find 
                   the astrometry table already present in the ObsFile instance. (So this
                   keyword should not normally be used or necessary.)  - DEPRECATED, 7/12/2013
                   IN FAVOUR OF PASSING THROUGH ObsFile instance.
        doIndex - boolean. If True, then index the main columns of the photon list table
                   for efficient access.
        pixRemapFileName - optionally set to a pathname for a pixel remapping .h5 file. If set,
                    will use the remap file to correct the x/y detector locations of the pixels
                    according to the map in the file before writing out to the photon list.
                    This is used to correct for errors in the beam map (specifically, for
                    the PAL2012 run). If set, the remap array will also be included in the saved
                    HDF file. BUT NOTE... hot pixel mask files based on the original mapping will 
                    probably not be great masks - mislocated pixels are likely to get masked out
                    before they can be remapped. In an attempt to fix this, pixels that are involved
                    in remapping have their hot pixel masks ignored, in theory. All is a bit of a
                    fudge just to try and get good Crab images out quickly for Matt's paper. Better
                    to remap the original data much earlier on in the whole process, and re-run the
                    calibrations accordingly etc. So this is not intended for use in the long term!
    
    """
    tic = time.clock()
    
    #Calculate unit conversion constant for use later when converting errors in eV
    #to errors in Angstroms
    toAngstromConst = ((( (1.*units.eV) * (1.*units.Angstrom)**2 ) /
                    (constants.h*constants.c) ).to(units.Angstrom).value)
    
    if pixRemapFileName is not None:
        warnings.warn('Pixel remapping during photon-list writing requested - not intended for general use, user beware!')
    
    if obsFile.flatCalFile is None: raise RuntimeError, "No flat cal. file loaded"
    if obsFile.fluxCalFile is None: raise RuntimeError, "No flux cal. file loaded"
    if obsFile.wvlCalFile is None: raise RuntimeError, "No wavelength cal. file loaded"
    if obsFile.hotPixFile is None: raise RuntimeError, "No hot pixel file loaded"
    if obsFile.centroidListFile is None: raise RuntimeError, "No astrometry centroid list file loaded"
    if obsFile.file is None: raise RuntimeError, "No obs file loaded...?"
    
    if filename is None:
        filename=FileName(obsFile=obsFile).photonList()
    
    print 'Initialising empty photon list file'
    plFile = createEmptyPhotonListFile(obsFile,filename,photListDescription=photListDescription)

    plTable = plFile.root.photons.photons
    
    #Copy all the various cal files used into the photon list file.        
    try:
        plFile.copyNode(obsFile.flatCalFile.root.flatcal, newparent=plFile.root, newname='flatcal', recursive=True)
        plFile.copyNode(obsFile.fluxCalFile.root.fluxcal, newparent=plFile.root, newname='fluxcal', recursive=True)
        plFile.copyNode(obsFile.wvlCalFile.root.wavecal, newparent=plFile.root, newname='wavecal', recursive=True)
        plFile.copyNode(obsFile.hotPixFile.root, newparent=plFile.root, newname='timeMask', recursive=True)
        plFile.copyNode(obsFile.file.root.beammap, newparent=plFile.root, newname='beammap', recursive=True)
        plFile.copyNode(obsFile.file.root.header, newparent=plFile.root, recursive=True)
        if obsFile.centroidListFile is not None:
            plFile.copyNode(obsFile.centroidListFile.root, newparent=plFile.root, newname='centroidList', recursive=True)
        if obsFile.timeAdjustFile is not None:
            #If there's any time adjustment file associated, store that too in the photon list 
            #file, and also correct the exptime and unixtime in the output header info 
            #(which is generally not updated in the original obs file).
            plFile.copyNode(obsFile.timeAdjustFile.root,newparent=plFile.root,newname='timeAdjust',recursive=True)
            plFile.root.header.header.cols.exptime[0] = obsFile.getFromHeader('exptime')
            plFile.root.header.header.cols.unixtime[0] = obsFile.getFromHeader('unixtime')
            plFile.root.header.header.cols.jd[0] = obsFile.getFromHeader('jd')
        if pixRemapFileName is not None:
            pixRemapFile = tables.openFile(pixRemapFileName, 'r')
            try:
                plFile.copyNode(pixRemapFile.root, newparent=plFile.root, newname='pixRemap', recursive=True)
            finally:
                pixRemapFile.close()
    except:
        plFile.flush()
        plFile.close()
        raise
    
    plFile.flush()


    try:
        
        fluxWeights = obsFile.fluxWeights      #Flux weights are independent of pixel location.
        #Extend flux weight/flag arrays as for flat weight/flags.
        fluxWeights = np.hstack((fluxWeights[0],fluxWeights,fluxWeights[-1]))
        fluxFlags = np.hstack((pipelineFlags.fluxCal['belowWaveCalRange'], 
                               obsFile.fluxFlags, 
                               pipelineFlags.fluxCal['aboveWaveCalRange']))
    
        #Initialise CalculateRaDec object with the right centroiding file for calculating RA/dec of photons.
        raDecCalcObject = crd.CalculateRaDec(obsFile.centroidListFile)
    
        #Intialise pixel remapping object if necessary:
        if pixRemapFileName is not None:
            pixMap = remapPixels.PixelMap(pixRemapFileName)
            pixMapSourceList,pixMapDestList = pixMap.getRemappedPix()
        else:
            pixMap = None
    
        #Make a numpy structured array dtype from the photon-list description header.
        #Currently uses somewhat undocumented hack in order to easily get a numpy 'dtype' from the photListDescription  class.
        #See http://pytables.github.io/usersguide/libref/declarative_classes.html?highlight=isdescription#tables.IsDescription
        # and specifically http://pytables.github.io/_modules/tables/description.html#dtype_from_descr
        #Newer versions of pytables (some time after v2.3.1?) should have a dtype_from_descr() function to pull this out directly....
        photListDescriptionObject = photListDescription()
        photListDtype = tables.Description(photListDescription.columns)._v_dtype

        wvlGoodFlag = pipelineFlags.waveCal['good']     #Just to avoid the dictionary lookup within the loop.
        belowWaveCalRangeFlag = pipelineFlags.flatCal['belowWaveCalRange']  #Ditto
        aboveWaveCalRangeFlag = pipelineFlags.flatCal['aboveWaveCalRange']  #Ditto
        for iRowRaw in xrange(obsFile.nRow):
            for iColRaw in xrange(obsFile.nCol):
                
                #Remap pixel row/column if necessary
                if pixMap is None:
                    iRowCorr = iRowRaw
                    iColCorr = iColRaw
                else:
                    iRowCorr,iColCorr = pixMap.remapPix(iRowRaw,iColRaw)
                print 'Pixel row, column (raw/remapped): ', iRowRaw, iColRaw, iRowCorr, iColCorr
                #print 'Pixel row, column (remapped): ', iRowCorr, iColCorr 
                
                #Now be careful to use raw and corrected values in the right places....
                flag = obsFile.wvlFlagTable[iRowRaw, iColRaw]
                if flag == wvlGoodFlag:     #only write photons in good pixels
                    energyError = obsFile.wvlErrorTable[iRowRaw, iColRaw] #Note wvlErrorTable is in eV !! Assume constant across all wavelengths. Not the best approximation, but a start....
                    flatWeights = obsFile.flatWeights[iRowRaw, iColRaw]
                    #Extend flat weight and flag arrays at beginning and end to include out-of-wavelength-calibration-range photons.
                    flatWeights = np.hstack((flatWeights[0],flatWeights,flatWeights[-1]))
                    flatFlags = np.hstack((belowWaveCalRangeFlag,
                                           obsFile.flatFlags[iRowRaw, iColRaw],
                                           aboveWaveCalRangeFlag))
                    
                    
                    if pixMap is not None:
                        #Fudge to make it ignore the hot pixel masking for pixels which are to be 
                        #remapped. This is really not the recommended way to do things....
                        if (iRowRaw,iColRaw) in pixMapSourceList:
                            print 'Switching off hot pix mask for this pixel, since it''s being remapped'
                            obsFile.switchOffHotPixTimeMask()
                        else:
                            #print 'Hot pix mask on'
                            obsFile.switchOnHotPixTimeMask()
                    
                    
                    #print 'Reading from obs. file...'
                    x = obsFile.getPixelWvlList(iRowRaw,iColRaw,excludeBad=False,dither=True,firstSec=firstSec,
                                             integrationTime=integrationTime)
                    #print 'Done reading'
                    
                    timestamps, wavelengths = x['timestamps'], x['wavelengths']     #Wavelengths in Angstroms
                    
                    if len(timestamps)==0:
                        print 'No photons for this pixel (probably hot pixel masked)'
                        continue     #If there's no photons for this pixel, don't waste any more time on it.
                    
                    #Convert errors in eV to errors in Angstroms (see notebook, May 7 2013). Prob should really
                    #avoid the fancy units package just to speed things up a bit here....
                    print 'Converting errors to Angstroms'
                    #wvlErrorsOld = ((( (energyError*units.eV) * (wavelengths*units.Angstrom)**2 ) /
                    #                (constants.h*constants.c) )
                    #             .to(units.Angstrom).value)
                    wvlErrors =  energyError * wavelengths**2  * toAngstromConst
                    #IF THIS ASSERTION WORKS, THEN DELETE THE WVLERRORSOLD CALCULATION! SECOND VERSION SHOULD BE 1000x FASTER.... 
                    #assert np.allclose(wvlErrors,wvlErrorsOld,rtol=1e-10,atol=1e-10)
                    #It does work, it seems...
                    
                    print 'Binning wavelengths...'
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
                    print 'Done binning wavelengths'
                    
                    iCols = np.empty(len(timestamps),dtype=int)
                    iRows = np.empty(len(timestamps),dtype=int)
                    iCols.fill(iColCorr)
                    iRows.fill(iRowCorr)
    
                    print 'Calculating RA/Decs...'
                    if obsFile.centroidListFile is not None:
                        ras,decs,hourAngles = raDecCalcObject.getRaDec(timestamps,np.array(iCols),np.array(iRows))
    
#                   if 1==1:     #(Temporary switch to use new code vs. old code - remove when satisfied it's working!)

                    #New method for writing -- avoid looping over each photon.   
                    #Create rows to append for this pixel in memory, *then* write out.
                    #print 'Appending photon block'
                    
                    #Create an empty table ('sub-list') in memory for this pixel:
                    newRows = np.zeros(len(timestamps), dtype=photListDtype)
                    
                    #And fill all the columns:
                    newRows['xPix'] = iColCorr  #***Check to see if it broadcasts properly - may need to broadcast manually***
                    newRows['yPix'] = iRowCorr
                    newRows['xyPix'] = xyPack(iRowCorr,iColCorr)
                    newRows['arrivalTime'] = timestamps
                    newRows['wavelength'] = wavelengths
                    newRows['waveError'] = wvlErrors
                    newRows['flatFlag'] = flatFlags[flatBinIndices]
                    newRows['flatWeight'] = flatWeights[flatBinIndices]
                    newRows['fluxFlag'] = fluxFlags[fluxBinIndices]
                    newRows['fluxWeight'] = fluxWeights[fluxBinIndices]
                    if obsFile.centroidListFile is not None:
                        newRows['ra']=ras
                        newRows['dec']=decs
                        newRows['ha']=hourAngles
                    else:
                        newRows['dec'] = np.nan     #***Check for proper broadcasting!***
                        newRows['ra'] = np.nan      # ""
                        newRows['ha'] = np.nan      # ""
                    plTable.append(newRows)
                    plTable.flush()
                    
#                     else:
#                         #Old row-by-row method.
#                         print 'Appending row-by-row'
#                         for iPhoton in xrange(len(timestamps)):
#                             #if wavelengths[iPhoton] > wvlRange[0] and wavelengths[iPhoton] < wvlRange[1] and binIndices[iPhoton] >= 0 and binIndices[iPhoton] < len(flatWeights):
#                             #create a new row for the photon list
#                             #print 'Photon #',iPhoton
#                             newRow = plTable.row
#                             newRow['xPix'] = iCol
#                             newRow['yPix'] = iRow
#                             newRow['xyPix'] = xyPack(iRow,iCol)
#                             newRow['arrivalTime'] = timestamps[iPhoton]
#                             newRow['wavelength'] = wavelengths[iPhoton]
#                             newRow['waveError'] = wvlErrors[iPhoton]
#                             newRow['flatFlag'] = flatFlags[flatBinIndices[iPhoton]]
#                             newRow['flatWeight'] = flatWeights[flatBinIndices[iPhoton]]
#                             newRow['fluxFlag'] = fluxFlags[fluxBinIndices[iPhoton]]
#                             newRow['fluxWeight'] = fluxWeights[fluxBinIndices[iPhoton]]
#                             if obsFile.centroidListFile is not None:
#                                 newRow['ra']=ras[iPhoton]
#                                 newRow['dec']=decs[iPhoton]
#                                 newRow['ha']=hourAngles[iPhoton]
#                             else:
#                                 newRow['dec'] = np.nan
#                                 newRow['ra'] = np.nan
#                                 newRow['ha'] = np.nan
#                             newRow.append()
                else:
                    print 'Bad wavelength flag, skipping data for this pixel'
    finally:
        print 'Flushing and closing...'
        plTable.flush()
        plFile.close()

    if doIndex is True:
        #Index the table that's just been written out.
        indexPhotonList(filename)
    
    print 'Done.'
    print 'Time taken (s): ', time.clock()-tic


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
    
    colsToIndex = ['xyPix', 'xPix', 'yPix', 'arrivalTime', 'wavelength', 'ra', 'dec']
    
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
        
        #Sort table by xyPix...
        #print 'Sorting table by xyPix'
        #sortedList = photList.copy(sortby='xyPix',overwrite=True,newname='sorted')
        #print 'Reindexing'
        #dummy = sortedList.colinstances['xyPix'].createCSIndex()
        print 'Done.'
    finally:
        photList.flush()
        if type(photFile is str):
            photListFile.close()
            
 
def xyPack(row,col):
    '''
    Pack an integer pair of pixel coordinates, x, y, into a single integer,
    stored in the 'XYpix' column of the photon list. Used for the purposes
    of access efficiency.
    Note - INPUTS MUST BE INTEGER! Otherwise you're gonna get weird answers.
    '''
    return xyPackMult * row + col

    
def xyUnpack(rowCol):
    '''
    Inverse of xyPack. Returns an integer tuple, (row,col).
    '''
    return rowCol // xyPackMult, rowCol % xyPackMult
    
    
