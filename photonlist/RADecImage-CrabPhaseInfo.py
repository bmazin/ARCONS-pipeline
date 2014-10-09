'''
Author: Julian van Eyken                    Date: May 15 2013

Package/class for handling of images created from photon lists that are derotated and
mapped to sky coordinates, and stacked.

(NB - NEED TO DOCUMENT OBJECT PROPERTIES)

'''
import gc
import time
import numpy as np
import tables
import matplotlib.pyplot as mpl
import pyfits
#import hotpix.hotPixels as hp
from util import utils
from util.popup import PopUp, plotArray
import photonlist.photlist as pl
from photonlist import boxer
from astrometry import CalculateRaDec as crd
from headers import pipelineFlags
import astropy.constants

from PyQt4 import QtGui
from PyQt4 import QtCore
from multiprocessing import Process
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

class RADecImage(object):
    h = astropy.constants.h.to('eV s').value  #4.135668e-15 #eV s
    c = astropy.constants.c.to('m/s').value   #'2.998e8 #m/s
    angstromPerMeter = 1e10 
    '''
    Class to hold derotated, integrated and possibly stacked images, in the sky coordinate
    frame.
    
    NOTE - to be able to import boxer for this routine,  
    boxer.so (exact name may vary by platform) needs to be created
    from fortran boxer.f (in the same directory) using command:

    f2py -c -m boxer boxer.f

    This command should work on any platform, I think, though the
    output file may have a different file extension. In any case,
    should all then be importable in python using 'import boxer', 
    and this module should all then run fine (see also readme.txt
    in this directory).
    '''
    
    def __init__(self,photList=None,nPixRA=None,nPixDec=None,cenRA=None,cenDec=None,
                 vPlateScale=0.1, detPlateScale=None, firstSec=0, integrationTime=-1,
                 doWeighted=True,wvlMin=None,wvlMax=None):
                 #expWeightTimeStep=1.0):
        '''
        Initialise a (possibly empty) RA-dec coordinate frame image.

        INPUTS:
            photList: optionally provide a PhotList object from which to create an
                        image (see photonlist.photlist)
            nPixRA, nPixDec: integers, Number of pixels in RA and Dec directions
                        for the virtual image.
            cenRA, cenDec: floats, location of center of virtual image in RA and
                        dec (both in radians)
            vPlateScale: float, plate scale for virtual image (arcseconds per 
                        virtual image pixel). Note that the attribute, self.vPlateScale
                        is stored in *radians* per pixel, as is self.detPlateScale (plate
                        scale for the detector pixels).
            detPlateScale: override the assumed detector plate scale (arcseconds
                        per detector pixel)
            firstSec: float, time from beginning of photon-list file at which to
                        begin integration 
            integrationTime: float, length of time to integrate for in seconds. If
                        -1, integrate to end of photon list.
            doWeighted,wvlMin,wvlMax - passed on to loadImage() if 'photList' keyword
                        is provided. Otherwise ignored.
            
            #### DEPRECATED ##### 
            expWeightTimeStep: float, time step to use when calculating exposure
                        time weights for the virtual pixels (seconds).
            #####################
        '''
        
        self.nPixRA = nPixRA                    #No. of virtual pixels in RA direction
        self.nPixDec = nPixDec                  #No. of virtual pixels in dec. direction
        self.cenRA = cenRA                      #RA location of center of field (radians)
        self.cenDec = cenDec                    #Dec location of center of field (rad)
        self.vPlateScale = vPlateScale*2*np.pi/1296000          #No. of radians on sky per virtual pixel.
        self.imageIsLoaded = False              #Flag to indicate whether an actual image has been loaded yet.
        if detPlateScale is None:
            self.detPlateScale = crd.CalculateRaDec.platescale*2*np.pi/1296000      #Radians per detector pixel. ******For now - but this really needs reading in from the photon list file.
        else:
            self.detPlateScale = detPlateScale
        if nPixRA is not None and nPixDec is not None:
            self.image = np.empty((self.nPixDec,self.nPixRA),dtype=float)   #To take a (possibly stacked) image in virtual
            self.image.fill(np.nan)
            self.effIntTimes = np.empty_like(self.image)    #Effective integration times for each pixel, in seconds.
            self.effIntTimes.fill(np.nan)
            self.expTimeWeights = np.empty((self.nPixDec,self.nPixRA),dtype=float)  #Weights for each pixel in the virtual image to account for effective integration time on each pixel.
            self.expTimeWeights.fill(np.nan)
            self.gridRA = np.empty((self.nPixRA),dtype=float)       #Virtual pixel boundaries in the RA direction
            self.gridRA.fill(np.nan)
            self.gridDec = np.empty((self.nPixDec),dtype=float)     #Virtual pixel boundaries in the dec. direction.
            self.gridDec.fill(np.nan)
            self.totExpTime = np.nan                        #Total exposure time included in current image
            #self.expWeightTimeStep = expWeightTimeStep
        else:
            self.image = None
            self.effIntTimes = None
            self.expTimeWeights = None
            self.gridRA = None
            self.gridDec = None
            self.totExpTime = None
            #self.expWeightTimeStep = expWeightTimeStep
        if (cenRA is not None and cenDec is not None and vPlateScale is not None
            and nPixRA is not None and nPixDec is not None):
            self.setCoordGrid()            
        if photList is not None:
            self.loadImage(photList,firstSec=firstSec,integrationTime=integrationTime,
                           doWeighted=doWeighted,wvlMin=wvlMin,wvlMax=wvlMax)     
    
    
    def setCoordGrid(self):
        '''
        Establish RA and dec coordinates for pixel boundaries in the virtual pixel grid,
        given the number of pixels in each direction (self.nPixRA and self.nPixDec), the
        location of the centre of the array (self.cenRA, self.cenDec), and the plate scale
        (self.vPlateScale). 
        '''
        #self.gridRA = np.empty((self.nPixDec,self.nPixRA),dtype=float)
        #self.gridRA.fill(np.nan)
        #self.gridDec = np.empty((self.nPixDec,self.nPixRA),dtype=float)
        #self.gridDec.fill(np.nan)
        
        #Note - +1's are because these are pixel *boundaries*, not pixel centers:
        self.gridRA = self.cenRA + (self.vPlateScale*(np.arange(self.nPixRA+1) - ((self.nPixRA+1)//2)))
        self.gridDec = self.cenDec + (self.vPlateScale*(np.arange(self.nPixDec+1) - ((self.nPixDec+1)//2)))
    
   
    def loadImage(self,photList,firstSec=0,integrationTime=-1,wvlMin=None,wvlMax=None,
                  doStack=False,        #expWeightTimeStep=None, 
                  savePreStackImage=None, doWeighted=True, spectrum = False, phase = False):  #savePreStackImage is temporary for test purposes
        #error  may not be neccessary here - NEEEEEIIIIIILLLLL
        '''
        
        Build a de-rotated stacked image from a photon list (PhotList) object.
        If the RADecImage instance already contains an image, the new image is added to it.
        
        INPUTS:
            photList - a PhotList object from which to construct the image.
            firstSec - time from start of exposure to start the 'integration' for the image (seconds)
            integrationTime - duration of integration time to include in the image (in seconds; -1 or NaN => to end of exposure)
            wvlMin, wvlMax - min and max wavelengths of photons to include in the image (Angstroms).
            doStack - boolean; if True, then stack the image to be loaded on top of any image data already present.
            
            #### DEPRECATED - NOW GETS TIME STEPS STRAIGHT FROM CENTROID LIST FILES #####
            expWeightTimeStep - see __init__. If set here, overrides any value already set in the RADecImage object.
                                If the new image is being stacked on top of a current image, a new value can be
                                supplied that is different from the current image's value; but only the last value used
                                (i.e. the one supplied) will be stored in the class attribute.
            ################################
            
            wvlMin, wvlMax - set min and max wavelength cutoffs for photons to be loaded in.
            savePreStackImage - temporary fudge, set to a file-name to save the image out to a file prior to stacking.
            doWeighted - if True, includes flat and flux weighting (i.e. flatfielding and spectral response) factors from photons,
                                and rejects photons from pixels where the flatfield is bad at any wavelength within the requested
                                wavelength range (all if wvlMin/wvl Max not specified).
                                ****NOTE - FLUX WEIGHTING NOT FULLY TESTED -- but looks probably okay.****


        Added by Neil - 9/9/14
        
        spectrum - If True, this will add a spectral dimension to the virual image
        phase - If True, this will add a phase dimension to the virtual image
        '''
        
        #posErr = 0.8    #Approx. position error in arcsec (just a fixed estimate for now, will improve later)
        #posErr *= 2*np.pi/(60.*60.*360.)  #Convert to radians
        
        tic = time.clock()
        
        photTable = photList.file.root.photons.photons   #Shortcut to table
        #if expWeightTimeStep is not None:
        #    self.expWeightTimeStep=expWeightTimeStep
        
        if wvlMin is not None and wvlMax is None: wvlMax = np.inf
        if wvlMin is None and wvlMax is not None: wvlMin = 0.0
        
        #Figure out last second of integration
        obsFileExpTime = photList.header.cols.exptime[0]
        if integrationTime==-1 or firstSec+integrationTime > obsFileExpTime:
            lastSec = obsFileExpTime
        else:
            lastSec = firstSec+integrationTime
       
        #If virtual coordinate grid is not yet defined, figure it out.
        if self.gridRA is None or self.gridDec is None:
            #Find RA/dec range needed, taking advantage of the fact that the ra/dec columns are (or should be) indexed....
            print 'Finding RA/dec ranges' 
            self.raMin = photTable.cols.ra[photTable.colindexes['ra'][0]]
            self.raMax = photTable.cols.ra[photTable.colindexes['ra'][-1]]
            self.decMin = photTable.cols.dec[photTable.colindexes['dec'][0]]
            self.decMax = photTable.cols.dec[photTable.colindexes['dec'][-1]]
            self.cenRA = (self.raMin+self.raMax)/2.0
            self.cenDec = (self.decMin+self.decMax)/2.0
            #Set size of virtual grid to accommodate.
            if self.nPixRA is None:
                #+1 for round up; +1 because coordinates are the boundaries of the virtual pixels, not the centers.
                self.nPixRA = int((self.raMax-self.raMin)//self.vPlateScale + 2)     
            if self.nPixDec is None:
                self.nPixDec = int((self.decMax-self.decMin)//self.vPlateScale + 2)
            self.setCoordGrid()
            
        #Short-hand notations for no. of detector and virtual pixels, just for clarity:
        nDPixRow,nDPixCol = photList.nRow,photList.nCol
        nVPixRA,nVPixDec = self.nPixRA,self.nPixDec
        
        #Calculate ratio of virtual pixel area to detector pixel area
        vdPixAreaRatio = (self.vPlateScale/self.detPlateScale)**2
        
        #Make a boolean mask of dead (non functioning for whatever reason) pixels
        #True (1) = good; False (0) = dead 
        #First on the basis of the wavelength cals:
        wvlCalFlagImage = photList.getBadWvlCalFlags()
        deadPixMask = np.where(wvlCalFlagImage == pipelineFlags.waveCal['good'], 1, 0)   #1.0 where flag is good; 0.0 otherwise. (Straight boolean mask would work, but not guaranteed for Python 4....)

        #Next on the basis of the flat cals (or all ones if weighting not requested)
        if doWeighted:
            flatCalFlagArray = photList.file.root.flatcal.flags.read()       # 3D array - nRow * nCol * nWavelength Bins.
            flatWvlBinEdges = photList.file.root.flatcal.wavelengthBins.read()   # 1D array of wavelength bin edges for the flat cal.
            lowerEdges = flatWvlBinEdges[0:-1]
            upperEdges = flatWvlBinEdges[1:]
            if wvlMin is None and wvlMax is None:
                inRange = np.ones(len(lowerEdges),dtype=bool)   # (all bins in range implies all True)
            else:
                inRange = ((lowerEdges >= wvlMin) & (lowerEdges < wvlMax) |
                           (upperEdges >= wvlMin) & (lowerEdges < wvlMax))
            flatCalMask = np.where(np.all(flatCalFlagArray[:,:,inRange]==False, axis=2), 1, 0) # Should be zero where any pixel has a bad flag at any wavelength within the requested range; one otherwise. Spot checked, seems to work.
        else:
            flatCalMask = np.ones((nDPixRow,nDPixCol))
        
        #Get the photons
        print 'Getting photon coords'
        print 'wvlMin, wvlMax: ',wvlMin,wvlMax
        if wvlMin is None:
            assert wvlMin is None and wvlMax is None
            print '(getting all wavelengths)'
            #tic = time.clock()
            photons = photTable.readWhere('(arrivalTime>=firstSec) & (arrivalTime<=lastSec)')
            #print 'v1 time taken (s): ', time.clock()-tic
            #tic = time.clock()
            #photons = np.array([row.fetch_all_fields() for row in photTable.where('(arrivalTime>=firstSec) & (arrivalTime<=lastSec)')])
            #photIndices = photTable.getWhereList('(arrivalTime>=firstSec) & (arrivalTime<=lastSec)')
            #print 'v2 time taken (s): ', time.clock()-tic
            #print 'Doing by second method'
            #tic = time.clock()
            #photons2 = [x for x in photons.iterrows() if (x['arrivalTime']>=firstSec) and (x['arrivalTime']<=lastSec)]
            #print 'Time taken (s): ',time.clock()-tic
        else:
            assert wvlMin is not None and wvlMax is not None
            print '(trimming wavelength range) '
            photons = photTable.readWhere('(arrivalTime>=firstSec) & (arrivalTime<=lastSec) & (wavelength>=wvlMin) & (wavelength<=wvlMax)')
            

        #Filter out photons to be masked out on the basis of detector pixel######
        print 'Finding bad detector pixels...'
        detPixMask = deadPixMask * flatCalMask      #Combine wave cal pixel mask and flat cal mask (should be the same in an ideal world, but not 
        whereBad = np.where(detPixMask == 0)
        badXY = pl.xyPack(whereBad[0],whereBad[1])   #Array of packed x-y values for bad pixels (CHECK X,Y THE RIGHT WAY ROUND!)
        allPhotXY = photons['xyPix']                 #Array of packed x-y values for all photons           
        #Get a boolean array indicating photons whose packed x-y coordinate value is in the 'bad' list.
        toReject = np.where(np.in1d(allPhotXY,badXY))[0]      #Zero to take index array out of the returned 1-element tuple.
        #Chuck out the bad photons
        print 'Rejecting photons from bad pixels...'
        tempphot = np.delete(photons,toReject)      #Create temp. variable to try and handle memory smartly.
        print 'deleting photons/grbage collection'        
        del photons
        gc.collect()
        photons = tempphot
        del tempphot        
        gc.collect()
        del whereBad
        del badXY
        del allPhotXY
        del toReject
        gc.collect()
        photX = photons['xPix']#[interVal[i-1]:interVal[i]]
        photY = photons['yPix']#[interVal[i-1]:interVal[i]]
        #########################################################################
        #photX and block used to start here
        

        if spectrum is not False:
            wvlBinEdges = RADecImage.makeWvlBins()
            if phase is not False:
                phaseBinEdges = RADecImage.makePhaseBins()
                ImageStack = np.zeros((len(self.gridDec)-1,len(self.gridRA)-1,len(wvlBinEdges)-1,len(phaseBinEdges)-1))
            else:
                ImageStack = np.zeros((len(self.gridDec)-1,len(self.gridRA)-1,len(wvlBinEdges)-1))
        
        elif phase is not False:
                phaseBinEdges = RADecImage.makePhaseBins()
                ImageStack = np.zeros((len(self.gridDec)-1,len(self.gridRA)-1,len(phaseBinEdges)-1))

        else:
            ImageStack = np.zeros((len(self.gridDec)-1,len(self.gridRA)-1))

        
        numPhot = len(photons)
        interVal = [0,int(numPhot/10.),int(numPhot/5.),int((3./10.)*numPhot),int((2./5.)*numPhot),int((1./2.)*numPhot), int((3./5.)*numPhot), int((7./10.)*numPhot), int((4./5.)*numPhot), int((9./10.)*numPhot), int(numPhot-1)]   

        print 'beginning the image construction'
        #self.hist_list = []
        #self.wvl_list = []
        #self.phase_list = []
        for i in range(1, len(interVal)):

            #photX = photons['xPix'][interVal[i-1]:interVal[i]]
            #photY = photons['yPix'][interVal[i-1]:interVal[i]]
            photRAs = photons['ra'][interVal[i-1]:interVal[i]] #[:1000000000] #Read all photon coords into an RA and a dec array.
            photDecs = photons['dec'][interVal[i-1]:interVal[i]]
            photHAs = photons['ha'][interVal[i-1]:interVal[i]]       #Along with hour angles...
            photWeights = photons['flatWeight'][interVal[i-1]:interVal[i]] * photons['fluxWeight'][interVal[i-1]:interVal[i]]   #********EXPERIMENTING WITH ADDING FLUX WEIGHT - NOT FULLY TESTED, BUT SEEMS OKAY....********
            print 'INCLUDING FLUX WEIGHTS!'
            photWavelengths = photons['wavelength'][interVal[i-1]:interVal[i]]
            if phase is not False:
        
                photPhases = photons['phase'][interVal[i-1]:interVal[i]]

        
            if wvlMin is not None or wvlMax is not None:
                assert all(photWavelengths>=wvlMin) and all(photWavelengths<=wvlMax)
            print 'Min, max photon wavelengths found: ', np.min(photWavelengths), np.max(photWavelengths)
            nPhot = len(photRAs)
        
        
            #Add uniform random dither to each photon, distributed over a square 
            #area of the same size and orientation as the originating pixel at 
            #the time of observation (assume RA and dec are defined at center of pixel).
            xRand = np.random.rand(nPhot)*self.detPlateScale-self.detPlateScale/2.0
            yRand = np.random.rand(nPhot)*self.detPlateScale-self.detPlateScale/2.0       #Not the same array!
            ditherRAs = xRand*np.cos(photHAs) - yRand*np.sin(photHAs)
            ditherDecs = yRand*np.cos(photHAs) + xRand*np.sin(photHAs)
        
            photRAs=photRAs+ditherRAs
            photDecs=photDecs+ditherDecs
        
            #Make the image for this integration
            if doWeighted:
                print 'Making weighted image'
                if spectrum is not False:
                     #added by Neil
                    if phase is not False:
                        print 'Gathering phase and spectral information'
                        
                        #assert 1==0
                        
                        thisImage, edges = np.histogramdd([photDecs,photRAs,photWavelengths,photPhases],[self.gridDec,self.gridRA,wvlBinEdges, phaseBinEdges],weights=photWeights)  #spectrum and phase
                        
                         
                        self.wvlEdges = edges[2]
                        assert self.wvlEdges.all() == wvlBinEdges.all()
                        self.phaseEdges = edges[3]
                        assert self.phaseEdges.all() == phaseBinEdges.all()

                        
                    else:
                        print 'getting spectral information'
                        thisImage, edges = np.histogramdd([photDecs,photRAs,photWavelengths],[self.gridDec,self.gridRA,wvlBinEdges],weights=photWeights)  #spectrum only
                        self.wvlEdges = edges[2]
                   
                elif phase is not False:
                    'getting phase information'
                    thisImage, edges = np.histogramdd([photDecs,photRAs,photPhases],[self.gridDec,self.gridRA,phaseBinEdges],weights=photWeights)  #phase only
                    self.phaseEdges = edges[2]
    
                else:
                    thisImage,thisGridDec,thisGridRA = np.histogram2d(photDecs,photRAs,[self.gridDec,self.gridRA], weights=photWeights) #neither spectrum nor phase information
        
            else:
            
                print 'Making unweighted image'
                if spectrum is not False:
                    wvlBinEdges = RADecImage.makeWvlBins()
                    if phase is not False:
                        print 'Gathering phase and spectral information'
                        phaseBinEdges = RADecImage.makePhaseBins()
    
                        thisImage, edges = np.histogramdd([photDecs,photRAs,photWavelengths,photPhases],[self.gridDec,self.gridRA,wvlBinEdges, phaseBinEdges]) #spectrum nd phase
                                                
                        self.wvlEdges = edges[2]
                        self.phaseEdges = edges[3]
    
                    else:
                        print 'getting spectral information'
                        thisImage, edges = np.histogramdd([photDecs,photRAs,photWavelengths],[self.gridDec,self.gridRA,wvlBinEdges])  ###watch for anywhere with thisGridRA and thisGridDec
                        self.wvlEdges = edges[2]  #spectrum only
                elif phase is not False:
                    print 'getting phase information'
                    thisImage, edges = np.histogramdd([photDecs,photRAs,photPhases],[self.gridDec,self.gridRA,phaseBinEdges]) #phase only
                    self.phaseEdges = edges[2]
                      
                else:
                    print 'no additional info'
                    thisImage,thisGridDec,thisGridRA = np.histogram2d(photDecs,photRAs,[self.gridDec,self.gridRA]) #nothing
        
            
        
        #Save the time slice images in detector coordinates if image saving is requested.        
            
            
            print 'adding slice'
            ImageStack+=thisImage
            print 'thisImage ', np.shape(thisImage)
            

        thisImage = ImageStack
        del ImageStack
        gc.collect()
        print 'thisImage ', np.shape(thisImage) 
        #------------
        #Time masking
        #------------
        
        if savePreStackImage is not None:
            saveName = 'det-'+savePreStackImage
            print 'Making detector-frame image slice for diagnostics: '+saveName
            detImSlice = np.histogram2d(photY,photX,bins=[photList.nRow,photList.nCol])[0]  #CHECK NROW/NCOL IS THE RIGHT WAY ROUND!
            mpl.imsave(fname=saveName,arr=detImSlice,origin='lower',
                        cmap=mpl.cm.gray,vmin=np.percentile(detImSlice, 0.5), vmax=np.percentile(detImSlice,99.5))
        del photX
        del photY
        gc.collect()
        
        #If hot pixels time-mask data not already parsed in, then parse it.
        if photList.hotPixTimeMask is None:
            photList.parseHotPixTimeMask()      #Loads time mask dictionary into photList.hotPixTimeMask
         
        #And start figuring out the exposure time weights....            
        print 'Calculating effective exposure times'
        
        #First find start/end times of each timestep ('frame') for calculating effective exp. times
        #Use the same timesteps as used in calculating the astrometry.
        
        #tStartFrames = np.arange(start=firstSec,stop=lastSec,
        #                         step=self.expWeightTimeStep)
        #tEndFrames = (tStartFrames+self.expWeightTimeStep).clip(max=lastSec)    #Clip so that the last value doesn't go beyond the end of the exposure.
        if 'centroidlist' in photList.file.root.centroidList:
            tStartFramesAll = np.array(photList.file.root.centroidList.centroidlist.times.read()) #Convert to array, since it's saved as a list.
        else:
            tStartFramesAll = np.array(photList.file.root.centroidList.times.read())   #For back compatibility
        tEndFramesAll = np.append(tStartFramesAll[1:], np.inf)                   #Last frame goes on forever as far as we know at the moment
        withinIntegration = ((tStartFramesAll < lastSec) & (tEndFramesAll > firstSec))
        tStartFrames = tStartFramesAll[withinIntegration].clip(min=firstSec)   #Now clip so that everything is within the requested integration time.
        tEndFrames = tEndFramesAll[withinIntegration].clip(max=lastSec)
        nFrames = len(tStartFrames)
        assert nFrames > 0      #Otherwise we have a problem....
        assert np.all(tStartFrames <= lastSec) and np.all(tEndFrames >= firstSec)
        
        #Get x,y locations of detector pixel corners (2D array of each x,y value, in detector space)
        #Assume definition where integer values represent location of pixel center.
        dPixXmin = np.indices((nDPixRow,nDPixCol))[1] - 0.5
        dPixXmax = np.indices((nDPixRow,nDPixCol))[1] + 0.5
        dPixYmin = np.indices((nDPixRow,nDPixCol))[0] - 0.5
        dPixYmax = np.indices((nDPixRow,nDPixCol))[0] + 0.5
        dPixXminFlat = dPixXmin.flatten()   #Flattened versions of the same since getRaDec() only works on flat arrays.
        dPixXmaxFlat = dPixXmax.flatten()
        dPixYminFlat = dPixYmin.flatten()
        dPixYmaxFlat = dPixYmax.flatten()
        
        #Create (1D) arrays for normalised center locations of virtual pixel grid (=index numbers, representing location of unit squares)
        vPixRANormCen = np.arange(nVPixRA)   #np.indices(nVPixDec,nVPixRA)[1]
        vPixDecNormCen = np.arange(nVPixDec) #np.indices(nVPixDec,nVPixRA)[0]
        
        #Create 1D arrays marking edges of virtual pixels (in 'normalised' space...)
        vPixRANormMin = np.arange(nVPixRA)-0.5
        vPixRANormMax = np.arange(nVPixRA)+0.5
        vPixDecNormMin = np.arange(nVPixDec)-0.5
        vPixDecNormMax = np.arange(nVPixDec)+0.5
        
        #Find origin of virtual array (center of virtual pixel 0,0) in RA/dec space.
        vPixOriginRA = np.mean(self.gridRA[0:2])     
        vPixOriginDec = np.mean(self.gridDec[0:2])
        vPixSize = self.vPlateScale       #Short hand, Length of side of virtual pixel in radians (assume square pixels)
        
        #Make array to take the total exposure times for each virtual pixel at each time step
        vExpTimesStack = np.zeros((nVPixDec,nVPixRA,nFrames))
        #vExpTimesStack2 = np.zeros((nVPixDec,nVPixRA,nFrames))  #FOR TEST PURPOSES
        
        #And one for the total exposure time at each pixel summed over all time steps
        vExpTimes = np.zeros((nVPixDec,nVPixRA))
        
        #Array to hold list of (equal) timestamps for each pixel at each timestep
        #(just for calculating the RA/dec coordinates of the pixel corners)
        frameTimeFlat = np.zeros((nDPixRow*nDPixCol))   #Also flat array for the purposes of getRaDec()
        frameTimeFlat.fill(np.nan)
        
        #Initialise RA/dec calculations of pixel locations for exposure time weighting
        raDecCalcObject = crd.CalculateRaDec(photList.file.root.centroidList)            
         
        #------------ Loop through the time steps ----------
        for iFrame in range(nFrames):
            
            print 'Time slice: ',iFrame+1, '/', nFrames

            #Calculate detector pixel corner locations in RA/dec space (needs to be clockwise in RA/dec space! (checked, gives +ve answers).
            frameTimeFlat.fill(tStartFrames[iFrame])
            dPixRA1,dPixDec1,dummy = raDecCalcObject.getRaDec(frameTimeFlat,dPixXminFlat,dPixYminFlat)      #dPix* should all be flat
            dPixRA2,dPixDec2,dummy = raDecCalcObject.getRaDec(frameTimeFlat,dPixXminFlat,dPixYmaxFlat)   
            dPixRA3,dPixDec3,dummy = raDecCalcObject.getRaDec(frameTimeFlat,dPixXmaxFlat,dPixYmaxFlat)
            dPixRA4,dPixDec4,dummy = raDecCalcObject.getRaDec(frameTimeFlat,dPixXmaxFlat,dPixYminFlat)
            
            #Reshape the flat-array results into arrays matching the detector shape.
            #Default ordering for reshape should just be the reverse of flatten().
            #(Note all this can probably be avoided by just using flat arrays throughout
            # - this is just a bit more intuitive this way at the moment).
            #dPixRA1,dPixDec1 = dPixRA1Flat.reshape(detShape),dPixDec1Flat.reshape(detShape)
            #dPixRA2,dPixDec2 = dPixRA2Flat.reshape(detShape),dPixDec2Flat.reshape(detShape)
            #dPixRA3,dPixDec3 = dPixRA3Flat.reshape(detShape),dPixDec3Flat.reshape(detShape)
            #dPixRA4,dPixDec4 = dPixRA4Flat.reshape(detShape),dPixDec4Flat.reshape(detShape)

            #Normalise to scale where virtual pixel size=1 and origin is the origin of the virtual pixel grid
            dPixNormRA1 = (dPixRA1 - vPixOriginRA)/vPixSize     #dPixNorm* should all be flat.
            dPixNormRA2 = (dPixRA2 - vPixOriginRA)/vPixSize
            dPixNormRA3 = (dPixRA3 - vPixOriginRA)/vPixSize
            dPixNormRA4 = (dPixRA4 - vPixOriginRA)/vPixSize
            dPixNormDec1 = (dPixDec1 - vPixOriginDec)/vPixSize
            dPixNormDec2 = (dPixDec2 - vPixOriginDec)/vPixSize
            dPixNormDec3 = (dPixDec3 - vPixOriginDec)/vPixSize
            dPixNormDec4 = (dPixDec4 - vPixOriginDec)/vPixSize
                
            #Get min and max RA/decs for each of the detector pixels    
            dPixCornersRA = np.array([dPixNormRA1,dPixNormRA2,dPixNormRA3,dPixNormRA4])      #2D array, 4 by nRow*nCol - should be clockwise, I think!
            dPixCornersDec = np.array([dPixNormDec1,dPixNormDec2,dPixNormDec3,dPixNormDec4])
            #dPixCornersRA = np.array([dPixNormRA4,dPixNormRA3,dPixNormRA2,dPixNormRA1])      #2D array, 4 by nRow*nCol - reversed, but gives -ve results, so prob. anti-clockwise....
            #dPixCornersDec = np.array([dPixNormDec4,dPixNormDec3,dPixNormDec2,dPixNormDec1])
            dPixRANormMin = dPixCornersRA.min(axis=0)     #Flat 1D array, nRow * nCol
            dPixRANormMax = dPixCornersRA.max(axis=0)
            dPixDecNormMin = dPixCornersDec.min(axis=0)
            dPixDecNormMax = dPixCornersDec.max(axis=0)

            #Get array of effective exposure times for each detector pixel based on the hot pixel time mask
            #Multiply by the bad pixel mask and the flatcal mask so that non-functioning pixels have zero exposure time.
            #Flatten the array in the same way as the previous arrays (1D array, nRow*nCol elements).
            #detExpTimes = (hp.getEffIntTimeImage(photList.hotPixTimeMask, integrationTime=tEndFrames[iFrame]-tStartFrames[iFrame],
            #                                     firstSec=tStartFrames[iFrame]) * detPixMask).flatten()
            detExpTimes = (photList.hotPixTimeMask.getEffIntTimeImage(firstSec=tStartFrames[iFrame], integrationTime=tEndFrames[iFrame]-tStartFrames[iFrame]) * detPixMask).flatten()
                     
            #Loop over the detector pixels.... (should be faster than looping over virtual pixels)
            for iDPix in np.arange(nDPixRow * nDPixCol):
                    #Find the pixels which are likely to be overlapping (note - could do this as a sorted search to make things faster)
                    maybeOverlappingRA = np.where((dPixRANormMax[iDPix] > vPixRANormMin) & (dPixRANormMin[iDPix] < vPixRANormMax))[0]
                    maybeOverlappingDec = np.where((dPixDecNormMax[iDPix] > vPixDecNormMin) & (dPixDecNormMin[iDPix] < vPixDecNormMax))[0]
                    
                    for overlapLocRA in maybeOverlappingRA:
                        for overlapLocDec in maybeOverlappingDec:
                            overlapFrac = boxer.boxer(overlapLocDec,overlapLocRA,dPixCornersDec[:,iDPix],dPixCornersRA[:,iDPix])
                            expTimeToAdd = overlapFrac*detExpTimes[iDPix]
                            vExpTimesStack[overlapLocDec,overlapLocRA,iFrame] += expTimeToAdd
           
            
        #------------ End loop through time steps ----------
                
        #Sum up the exposure times from each frame:
        if spectrum is not False:
            if phase is not False:
                print 'vExpTimesStacl ', np.shape(vExpTimesStack)
                vExpTimes = np.sum(vExpTimesStack,axis=2)
                print 'vExpTimes ', np.shape(vExpTimes)
                vExpTimes = vExpTimes[:,:,np.newaxis,np.newaxis]#(len(vExpTimes),len(vExpTimes[0]),1,1)
                print 'vExpTimes ', np.shape(vExpTimes)
                vExpTimes = np.tile(vExpTimes, (1,1, len(thisImage[0][0][:]), len(thisImage[0][0][0][:])))
                print 'vExpTimes ', np.shape(vExpTimes)
            else:
                vExpTimes = np.sum(vExpTimesStack,axis=2)
                vExpTimes = vExpTimes[:,:,np.newaxis]
                vExpTimes = np.tile(vExpTimes, (1,1,len(thisImage[0][0][:])))
        
        elif phase is not False:
            vExpTimes = np.sum(vExpTimesStack,axis=2)
            vExpTimes = vExpTimes[:,:,np.newaxis]#np.atleast_3d(vExpTimes)
            vExpTimes = np.tile(vExpTimes, (1,1,len(thisImage[0][0][:])))
        else:
            vExpTimes = np.sum(vExpTimesStack, axis = 2)

        print 'vExpTimes ', np.shape(vExpTimes)
        #Check that wherever the exposure time is zero, there are no photons that have not been rejected
        #assert np.all(thisImage[vExpTimes==0] == 0)
        #assert 1==0
        
        if savePreStackImage is not None:
            print 'Saving exp.time weighted pre-stacked image to '+savePreStackImage
            print 'cmap: ', mpl.cm.gray
            

            if spectrum is not False:
                if phase is not False:
                    imToSave = thisImage/vExpTimes
                    
                    print 'imToSave ', np.shape(imToSave)
                    imToSave = np.sum(imToSave, axis = 3)
                    print 'imToSave ', np.shape(imToSave)
                    imToSave = np.sum(imToSave, axis = 2)
                    print 'imToSave ', np.shape(imToSave)
                else:
                    imToSave = thisImage/vExpTimes
                    imToSave = np.sum(imToSave, axis = 2) #need to make this more general... NEil is spectrum not False....
            elif phase is not False:
                imToSave = thisImage/vExpTimes
                imToSave = np.sum(imToSave, axis = 2)            
            else:
                imToSave = thisImage/vExpTimes
           

            mpl.imsave(fname=savePreStackImage,arr=imToSave,origin='lower',cmap=mpl.cm.gray,
                       vmin=np.percentile(imToSave, 0.5), vmax=np.percentile(imToSave,99.5))
        
        if self.imageIsLoaded is False or doStack is False:
            self.image = thisImage           #For now, let's keep it this way.... Since weighting does odd things.
           
            self.effIntTimes = vExpTimes
            self.totExpTime = lastSec-firstSec
            self.expTimeWeights = self.totExpTime/self.effIntTimes
            self.vExpTimesStack = vExpTimesStack                   #TEMPORARY FOR DEBUGGING PURPOSES
            self.imageIsLoaded = True
        else:
            assert self.imageIsLoaded == True
            print 'Stacking'
            self.image += thisImage
            
            self.effIntTimes += vExpTimes
            self.totExpTime += lastSec-firstSec
            self.expTimeWeights = self.totExpTime/self.effIntTimes

        print 'Image load done. Time taken (s): ', time.clock()-tic


    def display(self,normMin=None,normMax=None,expWeight=True,pclip=None,colormap=mpl.cm.gnuplot2,
                image=None, logScale=False, fileName=None, DDArr = False, showColorBar = True):
        '''
        Display the current image. Currently just a short-cut to utils.plotArray,
        but needs updating to mark RA and Dec on the axes.
        
        INPUTS:
            normMin, normMax: Minimum and maximum values for the color-scale stretch
            expWeight: if True, scale each virtual pixel by its exposure time
            pclip: Apply a percentile clip to the image at the lower [pclip] and upper [100-pclip]'th
                    percentiles (where pclip is in percent). Overrides normMin and normMax.
            colorMap: as for plotArray, can specify the matplotlib color map to use.
            image: if set to a 2D array, displays that array instead of the default image.
            logScale: if True, display the intensities on a log scale.
            fileName: if a string, save the plot to this filename. If anything else (inc. None),
                      display to screen.
           

            Added by Neil - 9/1/14

            DDArr - If the histogrammed virtual Image has more than two spacial dimensions, set DDArr = True and the image will be displayed as in 2D space.
            showColorBar - set True for an adjustable color bar on the display. 
        '''
        
        if expWeight:
            toDisplay = np.copy(self.image*self.expTimeWeights)
        else:
            toDisplay = np.copy(self.image)
        
        if DDArr is not False:
            if len(np.shape(toDisplay)) == 3:

                toDisplay = np.sum(toDisplay, axis = 2)
            elif len(np.shape(toDisplay)) == 4:
                toDisplay = np.sum(toDisplay, axis = 3)
                toDisplay = np.sum(toDisplay, axis = 2)
                
        else:
            toDisplay = toDisplay        
        
        if logScale is True: toDisplay = np.log10(toDisplay)
        
        if image is not None: toDisplay = np.copy(image)
        

        
        #####################
        #Rotate by 180deg so north is up and east is left - just a fudge for now, need to
        #sort this out properly.
        toDisplay = np.rot90(np.rot90(toDisplay))
        #####################
        
        if pclip:
            normMin = np.percentile(toDisplay[np.isfinite(toDisplay)],q=pclip)
            normMax = np.percentile(toDisplay[np.isfinite(toDisplay)],q=100.0-pclip)

        #Display NaNs as zeros so it looks better
        toDisplay[np.isnan(toDisplay)] = 0
        
        self.toDisplay = toDisplay
        #Find the coordinates of the centers of the virtual pixels in degrees
        #raMin = (self.gridRA[0:-1] + self.gridRA[1:])/2.0 / np.pi * 180.
        #dec = (self.gridDec[0:-1] + self.gridDec[1:])/2.0 / np.pi * 180.
        
        #utils.plotArray(toDisplay,cbar=True,normMin=normMin,normMax=normMax,
        #                colormap=colormap,plotFileName=fileName,
        #                showMe=(type(fileName) is not str))
        #fig = mpl.figure()
        #ax = fig.add_subplot(111)
        #ax.set_xticklabels(ax.get_xticklabels(),rotation=90)
        
        # neil adding stuff chnging mpl to ax
        self.fig = mpl.figure()

        mpl.xticks(rotation=90)
        mpl.xlabel('R.A. (deg)')
        mpl.ylabel('Dec. (deg)')        
        
        self.ax = mpl.gca() #Neil
    
        self.ax.ticklabel_format(style='plain',useOffset=False)  #used to be mpl.tick
        self.ax.tick_params(direction='out')                        #mpl.tick 
        self.im = self.ax.imshow(self.toDisplay,vmin=normMin,vmax=normMax, extent=(180./np.pi)*   #mpl.imshow
                    np.array([self.gridRA[0],self.gridRA[-1],self.gridDec[0],self.gridDec[-1]]),
                    origin='upper', cmap=colormap)
        
        if showColorBar:
            self.fig.cbar = self.fig.colorbar(self.im)
            
            cid = self.fig.canvas.mpl_connect('scroll_event', self.onscroll_cbar)
            cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick_cbar)
        
       
        #mpl.ticklabel_format(style='plain',useOffset=False)
                
        #ax = mpl.gca()
        self.ax.showCoord = self.showCoord(ax=self.ax) #neil

        xtickBase = 10**(np.round(np.log10((self.gridRA[-1]-self.gridRA[0])*180./np.pi/10.)))
        ytickBase = 10**(np.round(np.log10((self.gridDec[-1]-self.gridDec[0])*180./np.pi/10.)))
        self.ax.xaxis.set_major_locator(mpl.MultipleLocator(base=xtickBase))
        self.ax.yaxis.set_major_locator(mpl.MultipleLocator(base=ytickBase))

        
        
        #### Neil adding stuff from Popup.
       
        #ax.xaxis.set_major_locator(mpl.MultipleLocator(base=0.001))
        #ax.yaxis.set_major_locator(mpl.MultipleLocator(base=0.001))
        #mpl.show()
       
    def draw(self):
        self.fig.canvas.draw()    

    def onscroll_cbar(self, event):
        '''
        added by Neil for adjustable color bar on display
        '''
        if event.inaxes is self.fig.cbar.ax:
            increment=0.05
            currentClim = self.fig.cbar.mappable.get_clim()
            currentRange = currentClim[1]-currentClim[0]
            if event.button == 'up':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]+increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]+increment*currentRange)
            if event.button == 'down':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]-increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]-increment*currentRange)
            self.fig.cbar.mappable.set_clim(newClim)
            self.fig.canvas.draw()

    def onclick_cbar(self,event):
        '''
        added by Neil for adjustable color bar on display
        '''
        if event.inaxes is self.fig.cbar.ax:
            self.fig.currentClim = self.fig.cbar.mappable.get_clim()
            lower = self.fig.currentClim[0]
            upper = self.fig.currentClim[1]
            fraction = event.ydata
            currentRange = upper-lower
            clickedValue = lower+fraction*currentRange
            extrapolatedValue = lower+event.ydata*currentRange
            if event.button == 1:
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (clickedValue,upper)
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (lower,clickedValue)
            if event.button == 3:
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = ((lower-fraction*upper)/(1.-fraction),upper)
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (lower,lower+currentRange/fraction)
            self.fig.cbar.mappable.set_clim(newClim)
            self.fig.canvas.draw()

    def showCoord(self, ax):
        '''
        added by neil
        ax - this is the ax instance used in display as fig, ax = mpl.subplots()
             ax.imshow....
        '''
        def format_coord(x,y):
            try:
                im = mpl.gca().get_images().get_array().data
                nrow,ncol = numpy.shape(im)
                row,col = int(y+0.5),int(x+0.5)
                if col>=0 and col<=ncol and row>=0 and row<nrow:
                    z = im[row,col]
                    return 'x=%1.7f, y=%1.7f, z=%1.7f'%(x, y, z)
                else:
                    return 'x=%1.7f, y=%1.7f, --'%(x, y)     
            except:
                return 'x=%1.7f, y=%1.7f, --'%(x, y)
        
        ax.format_coord = format_coord
               
    
    
        
    def writeFits(self, fileName='RADecImage.fits', expWeight=True):
        '''
        Write the current image (if present) to a fits file.
        
        INPUTS:
            fileName - name of output file to write
            expWeight - if True, scale by the per-pixel exposure time weights.
        '''
        
        if self.image is None:
            raise RuntimeError, 'No current image to write to FITS file'
        
        if expWeight is True:
            hdu = pyfits.PrimaryHDU(self.image*self.expTimeWeights)
        else:
            hdu = pyfits.PrimaryHDU(self.image)
        
        hdu.writeto(fileName)
        
    def getApertureSpectrum(self, cenRA, cenDec, nPixRA, nPixDec, radius1, radius2=None, degrees=False, error = False, offSet = False, radRA = None, radDec = None,  radius3 = None, phase = False, crab = False)                                                                         
        
        '''
        Added by Neil - 9/1/14

        cenRa/cenDec - the location of the center of the aperture
        nPixRa/nPixDec - the number of pixels in the virtual image 
        radius1/2/3 - depending on if the object of interest contains multiple objects these can be used to create an aperture on one part of the object 
                      and still have the sky subtraction be done using the radius of the whole thing. (must use in order)
        
        degrees - True if cenRA/cenDec are in degrees. This is fed to util.aperture
        error -  If true, include the count error on spectrum. Note this assumes the only error is in the photon counts.
        offSet - If True, use radRA/radDec and radius 2 and radius 3 to create sky suptracted spectrum (so far only really need this for EinsteinCross where multiple objects are within a small radius)
        radRa/radDec - position for the offset rings
        phase - This should only be used if the virtual image object contains phase informatio
        
        crab - this is somewhat preliminary. If this is True then that sky subtracted array containing the counts per wavelength bin per phase bin will be sorted by the phase range. This is hard coded to be the 
               phase ranges of the pulse and interpulse of Crab. This option can only be use if phase is True.
        

        note that as of right now the error analysis is a bit spotty.
        '''
        limits = [self.gridRA[0], self.gridRA[-1], self.gridDec[0], self.gridDec[-1]] #these are the limits of the grid in RA/Dec. They are in radians as is.
                                                                                      #RADecImage._init_ returns these in radians inherintly.
        
        print '----inner aperture----'
        aperture1 = utils.aperture(cenRA=cenRA,cenDec=cenDec,nPixRA=nPixRA,nPixDec=nPixDec,limits=limits,radius=radius1,degrees=degrees) #see utils.aperture
        innerAp = aperture1
        
        print '----annulus----'
        #aperture2 = utils.aperture(cenRA=cenRA,cenDec=cenDec,nPixRA=nPixRA,nPixDec=nPixDec,limits=limits,radius=radius2,degrees=degrees)
        
         
        #skyMask = aperture2 - aperture1
        
        if offSet is not False:   #for things like the einstein cross that have multiple objects that apppear as 1 object in the image
            if radRA and radDec is not None: #these are the offset RA and Dec. must be provided to get an aperture for offset object
                if radius2 and radius3 is not None: #must be provided inner and outer radius for sky subtraction
                    aperture2 = utils.aperture(cenRA=radRA,cenDec=radDec,nPixRA=nPixRA,nPixDec=nPixDec,limits=limits,radius=radius2,degrees=degrees)
                    aperture3 = utils.aperture(cenRA=radRA,cenDec=radDec,nPixRA=nPixRA,nPixDec=nPixDec,limits=limits,radius=radius3,degrees=degrees)
                    skyMask = aperture3 - aperture2 #sky subtraction. Note that this would be used if there was a small object in a cluster or something where you
                                                    #would want to use the objects center for the actual aperture and then set the center of the cluster as the offset
                                                    #RA and Dec and set the radius large enough so the sky spectrum can be obtained while still only having the aperture
                                                    #around the object of interest.
                else: 
                    ValueError('cannot do sky subtraction') #If offset is passed as a parameter and RA/Dec position is given but not enough radii are given for sky subtraction.
            else:
                raise ValueError('Need to provide center position for offset aperture') #If offset is passed as a parameter but no ofset position is provided.
        
        else:
            if radius2 and radius3 is not None:  #If offset is not indicated, and 3 radii given, there will be two anuli around the aperture.
                aperture2 = utils.aperture(cenRA=cenRA,cenDec=cenDec,nPixRA=nPixRA,nPixDec=nPixDec,limits=limits,radius=radius2,degrees=degrees) 
                aperture3 = utils.aperture(cenRA=cenRA,cenDec=cenDec,nPixRA=nPixRA,nPixDec=nPixDec,limits=limits,radius=radius3,degrees=degrees)
                skyMask = aperture3 - aperture2 #outer most ring subtracted by second ring to provide sky spectrum.
            elif radius2 is not None: #this is classic aperture and anulus. This will be used most of the time.
                aperture2 = utils.aperture(cenRA=cenRA,cenDec=cenDec,nPixRA=nPixRA,nPixDec=nPixDec,limits=limits,radius=radius2,degrees=degrees)
                skyMask = aperture2 - aperture1 #sky subtraction.
            else:
                raise ValueError('additional radius must be given to get sky subtracted spectrum') #if no aditional radii given there can be no sky subtraction.

      
        
        print '----locating pixel positions within aperture----'
        y1_values, x1_values = np.where(innerAp==0) #This gathers the pixel locations contained in the aperture
        self.display(image = innerAp) #quick check to see how the aperture position looks
        
        print '----locating pixel positions within annulus----'
        y_sky, x_sky = np.where(skyMask==0) #This gathers the pixel locations contained in the annulus. 
        self.display(image = skyMask) #quick look to see how the sky anulus looks
       
        #print y1_values, np.shape(y1_values)
        #print x1_values, np.shape(x1_values)
        #print y_sky, np.shape(y_sky)
        #print x_sky, np.shape(x_sky)
                


        print '----gathering virtual image information----'
        print '---finding wavelngth bins---'        
        wvlBinEdges = self.wvlEdges      #gets the wavelength binning information generated when making the image.
        if phase is not False:
            print '---finding phase bins---'
            phaseBinEdges = self.phaseEdges #gets phase binning information generated when making the image
        
        print '---making image---'
        self.scaledSpectrum = self.image*self.expTimeWeights #this is the virtual image Cube(counts/pixel location) scaled by the exposure time weightings!!.
                                                          #the units of self.scaledSpectrum are counts/pixel
                                                          #There may be a significant amount of NANS here so be careful.
        
        print 'self.scaledSpectrum', np.shape(self.scaledSpectrum), self.scaledSpectrum

        print '---finding dimensions of image---'
        #this is to make sure that the information we are trying to pull out of the image is consistant with the # of dimensions of the virtual image cube. 
        if phase is False: #image could have a phase dimension or not, but if it does and we have chosen not to pull it out, the image is summed along this dimension.
            if np.shape(self.scaledSpectrum) == 4:
                self.scaledSpectrum = np.sum(self.scaledSpectrum, axis = 3)
            else:
                assert np.shape(self.scaledSpectrum) == 3   #This is because if phase is not true all we will pull out is spectra which is the 3rd dimension. (#pixelsX,#pixelsY,#wavelength bins)       
        
        
        if error is not False:  #initialize error analysis
            self.err = np.sqrt(self.image) #the error on counting events is sqrt(N).
            self.spectrumErr = self.err/self.effIntTimes #scale error by integration times. There may be NANS here as well. 
            
        
            print 'self.spectrumErr', self.spectrumErr, np.shape(self.spectrumErr)
        
        

        if phase is not False:
            skyspectrum = np.zeros((len(x_sky), len(self.scaledSpectrum[x_sky[0]][y_sky[0]][:]), len(self.scaledSpectrum[x_sky[0]][y_sky[0]][0][:]))) #initializes array with dimensions (#pixels,#wavelengthbins,#phaseBins)
        else:
            skyspectrum = np.zeros((len(x_sky), len(self.scaledSpectrum[x_sky[0]][y_sky[0]][:]))) #initializes array that has dimensions of (#pixels,#wavelength bins)
                                                                                                  #This is a (#of pixels in annulus, #wavelenth bin per pixel) array. 
                                                                                                  #Note that len(x_sky)=len(y_sky)= #pixels
        print '---finding skyspectrum---'
        for i in range(len(x_sky)):     
            for j in range(len(self.scaledSpectrum[x_sky[i]][y_sky[i]][:])):
                if phase is not False: 
                    for k in range(len(self.scaledSpectrum[x_sky[i]][y_sky[i]][j][:])):
                        skyspectrum[i][j][k] = self.scaledSpectrum[x_sky[i]][y_sky[i]][j][k] #iterates through x_sky and builds new array that has spectra and phase corresponding to 
                                                                                             #number of pixels in annulus. dimensions (#pixels in ann, #wvlBins, #phasebins)
                else:
                    skyspectrum[i][j] = self.scaledSpectrum[x_sky[i]][y_sky[i]][j]   #iterates through x_sky and builds new array that has spectra corresponding to number of pixels in
                #print skyspectrum[i][j]                                             #in annulus on virtual grid. Array has dimensions (#of pixels in annulus, #wavelength bin per pixel)

        #print 'skyspectrum', skyspectrum, np.shape(skyspectrum)       
        if phase is not False:
            sky_array = np.zeros((len(self.scaledSpectrum[x_sky[0]][y_sky[0]][:]),len(self.scaledSpectrum[x_sky[0]][y_sky[0]][0][:]))) #initializes array for skysubtraction (#wvlbins, #phasebins)
        else: 
            sky_array = np.zeros(len(self.scaledSpectrum[x_sky[0]][y_sky[0]][:]))    #initialize array to be used for sky subtraction with dimensions of (#wavelength bins)

        if error is not False:
            sky_arrayErr = np.zeros(len(self.scaledSpectrum[x_sky[0]][y_sky[0]][:])) #initialize array of errors for count values per wavelength bin.
                                                                                   #has dimensions of (#wavelength bins per pixel)
        
        skyspectrum = skyspectrum.transpose()  #transpose skyspectrum so that dimensions are now (#wavelenth bin per pixel, #of pixels in annulus) -IF SPECTRUM ONLY IS USED
                                               #this makes it easier to take the median with NANs values. Now medians will be taken across a single wavelength bin.
                                               #IF INCLUDING PHASE INFORMATION!!!! ----this array will now be (#phasebins, #wvlbins, #pixels) where again the medians are taken for the pixels
                                               # across a wavelength bin.
        
        print '---finding median sky count for sky subtrction---'
        for i in range(len(skyspectrum)):
            if phase is not False:
                for j in range(len(skyspectrum[j])):
                    nanNot = ~np.isnan(skyspectrum[i][j]) #finds where the NANS are NOT.
                    sky_array[j][i] = np.median(skyspectrum[i][j][nanNot]) #takes median of counts across all positions in sky annulus per wavelength bin, excluding NANS. 
                                                                           #has dimensions (#wvlnbins, #phasebins) where each element contains the median count of all pixels. 
               
            else:
                nanNot = ~np.isnan(skyspectrum[i])  #finds where the NANS are NOT.
                sky_array[i] = np.median(skyspectrum[i][nanNot]) #takes median of counts across all positions in sky annulus per wavelength bin, excluding NANS. 
                                                                 #has dimensions (#wavelength bins per pixel) where each element contains the median count of all pixels.     
            if error is not False:
                sky_arrayErr[i] = np.std(skyspectrum[i][nanNot])  #takes the standard deviation of the values that went into finding the medium sky value for a wavelength bin
                                                                #and sets that as the error for that median count for that bin. 
                        

        print 'skyspectrum', skyspectrum, np.shape(skyspectrum)
        print 'sky_array', sky_array, np.shape(sky_array)

        #print 'creating sky subtracted spectrum'
        '''
        if error is not False: 
            skyspectrumErr = np.zeros((len(x_sky), len(self.spectrumErr[x_sky[0]][y_sky[0]][:]))) #initializes array with dimensions (#pixels in annulus, #wavelength bin per pixel)
            for i in range(len(x_sky)):
                for j in range(len(self.spectrumErr[x_sky[i]][y_sky[i]][:])):               #iterates through x_sky and builds array with dimension (#pixels in anulus, #wavelengh bin)
                    skyspectrumErr[i][j] = self.spectrumErr[x_sky[i]][y_sky[i]][j]
            
            sky_arrayErr = np.zeros(len(self.spectrumErr[x_sky[0]][y_sky[0]][:]))       #initializes array for sky ubtraction errors with dimension (#wavelength bin)
            
            skyspectrumErr = skyspectrumErr.transpose()     #transpose skyspectrumErr so dimensions are (#wavelength bins per pixel, #pixels in annulus)
            
            for i in range(len(skyspectrumErr)):        `    
                nanNot = ~np.isnan(skyspectrumErr[i])      #find where the NANS are not
                sky_arrayErr[i] = np.median(skyspectrumErr[i][nanNot])    #takes median
        '''
        #print 'sky_arrayErr', sky_arrayErr, np.shape(sky_arrayErr)

        if phase is not False:
            spectrumIn = np.zeros((len(x1_values),len(self.scaledSpectrum[x1_values[0]][y1_values[0]][:]),len(self.scaledSpectrum[x1_values[0]][y1_values[0]][0][:])))  # initializes arr (#pixels,#wvlbins,#phasebins)
        else:
            spectrumIn = np.zeros((len(x1_values),len(self.scaledSpectrum[x1_values[0]][y1_values[0]][:]))) #Initializes array of dimension (#pixels in aperture,#wavelength bin per pixel)
        
        print '---gathering aperture spectrum---'
        if phase is not False:
            for i in range(len(x1_values)):
                specIn = np.zeros((len(self.scaledSpectrum[x1_values[0]][y1_values[0]][:]), len(self.scaledSpectrum[x1_values[0]][y1_values[0]][0][:]))) #array is (#wvlbins, #phasebins)
                for j in range(len(self.scaledSpectrum[x1_values[i]][y1_values[i]][:])):
                    
                    for k in range(len(self.scaledSpectrum[x1_values[i]][y1_values[i]][j][:])):
                        specIn[j][k] = self.scaledSpectrum[x1_values[i]][y1_values[i]][j][k]  # fills in the array (#wvlbins, #phase bins)
                spectrumIn[i] = specIn - sky_array  #subtracts the median sky counts from the counts per wavelength bin for a single pixel
                                                #Note that spectrum in will have dimensions (#pixels in aperture, #wavelength bins per pixel, #phasebins)
                   
                
        else: 
            for i in range(len(x1_values)):  #iterates through each pixel contained in the aperture
                specIn = np.zeros(len(self.scaledSpectrum[x1_values[0]][y1_values[0]][:])) #initializes array of dimensions (#wavelength bins per pixel)
                for j in range(len(self.scaledSpectrum[x1_values[i]][y1_values[i]][:])):   #iterates through all wavelength bins for a single pixel
                    specIn[j] = self.scaledSpectrum[x1_values[i]][y1_values[i]][j]        #contains the counts per wavelength bin per pixel with dimensions (#wavelength bins per pixel)
                #print specIn
                #print 'specIn', specIn, np.shape(specIn)
                spectrumIn[i] = specIn - sky_array  #subtracts the median sky counts from the counts per wavelength bin for a single pixel
                                                #Note that spctrum in will have dimensions (#pixels in aperture, #wavelength bins per pixel)
        
        #print 'spectrumIn', spectrumIn, np.shape(spectrumIn)           
        print '---summing spectra for all pixels in aperture---'

        summed_arrayIn = np.sum(spectrumIn, axis = 0) #sums across all the pixels the number of counts per wavelength bin
                                                      #this array has dimension (#wavelength bins) where each element contains the total count of all positions for that bin
                                                      # If phase is True then this will have dimensions (#wavelengthbins,#phasebins)
        

        print 'summed_arrayIn', summed_arrayIn, np.shape(summed_arrayIn)
        
        if error is not False:
            summed_arrayInErr = np.std(spectrumIn, axis = 0)
        
            print 'summed_arrayInErr', summed_arrayInErr, np.shape(summed_arrayInErr)
        
        

        if phase is not False:
            spectrum = np.sum(summed_arrayIn, axis = 1)   #this takes the (#wvlbin,#phasebin) array and sums it along the phase dimension to get the spectrum
            phaseprofile = np.sum(summed_arrayIn,axis = 0) # this sums across the wavlength bins to get the phase profile

            if crab is not False:   #specific for crab, may need changing         
                #this finds the center of the phase bins. This is done to specify the phase ranges for the pulse and interpulse.
                nphaseBins = len(phaseBinEdges) - 1
   
                phases = np.empty((nphaseBins), dtype = float)
                for n in xrange(nphaseBins):                              
                    pbinsize = phaseBinEdges[n+1]-phaseBinEdges[n]
                    phases[n] = (phaseBinEdges[n] + (pbinsize/2.0))
                
                print '---determining spectrum of mainPulse and interpulse----'
                spectrumP = np.zeros((len(summed_arrayIn), len(np.where(np.logical_and(.60 <= phases,.72>=phases))[0])))   #initializes array of counts for the specified phase range             
                spectrumIP = np.zeros((len(summed_arrayIn), len(np.where(np.logical_and(.02 <= phases,.14>=phases))[0])))  #these arrays should have dimension (#wvlbins, #phasebins in phase range)            
                print 'whereP: ', np.where(np.logical_and(.60 <= phases,.72>=phases))[0]
                print 'whereIP: ', np.where(np.logical_and(.02 <= phases,.14>=phases))[0]                
                
                for i in range(len(summed_arrayIn)):

                    spectrumP[i][:] = summed_arrayIn[i][np.where(np.logical_and(.6 <= phases,.72>=phases))[0]]   #these srrays contain only the spectral information for the specified phase ranges.
                    spectrumIP[i][:] = summed_arrayIn[i][np.where(np.logical_and(.02 <= phases,.14>=phases))[0]] #arrays should have dimension (#wvlbins, #phasebins in phase range)
                    #print 'spectrumP[i] = ', spectrumP[i]
                    #print 'spectrumIP[i] = ', spectrumIP[i]

                print 'spectrumP before summing: ', np.shape(spectrumP), spectrumP
                print 'spectrumIP before summing: ', np.shape(spectrumIP), spectrumIP
                
                spectrumP = np.sum(spectrumP, axis = 1)  #sums across the phase bins to get the spectrum of the pulse in counts. dimensions (#wavlength bins)
                spectrumIP = np.sum(spectrumIP, axis = 1)      
                print 'spectrumP', np.shape(spectrumP), spectrumP
                print 'spectrumIP', np.shape(spectrumIP), spectrumIP
                
                for i in range(len(spectrumP)):
                    spectrumP[i]/=(wvlBinEdges[i+1]-wvlBinEdges[i])  #divides the counts at each wavelength bin by the binwidth => units are now counts/A
                    spectrumIP[i]/=(wvlBinEdges[i+1]-wvlBinEdges[i])
                    spectrum[i]/=(wvlBinEdges[i+1]-wvlBinEdges[i])
            else: #if Crab is false
    
                for i in range(len(spectrum)):
                    spectrum[i]/=(wvlBinEdges[i+1]-wvlBinEdges[i])  #divides the counts at each wavelength bin by the binwidth => units are now counts/A
                                    
                
        else: #if phase is False
            spectrum = summed_arrayIn
            for i in range(len(spectrum)):
                spectrum[i]/=(wvlBinEdges[i+1]-wvlBinEdges[i])  #divides the counts at each wavelength bin by the binwidth => units are now counts/A
                
                #summed_arrayInErr[i]/=(wvlBinEdges[i+1]-wvlBinEdges[i])
        #area = len(x1_values)*.02**2  #radius of aperture in cm
        #r = (len(x1_values)/2)*.02
        #area = np.pi*r**2

        #summed_arrayIn = summed_arrayIn/area  #now in counts/second/angstrom/area


        print '----completing aperture spectroscopy----'
        #print np.shape(wvlBinEdges)
        if phase is not False:
            if crab is not False:

                if error is not False:
                
                    return spectrum, spectrumP, spectrumIP, phaseprofile, wvlBinEdges, phaseBinEdges, summed_arrayInErr, innerAp 
                else:
                    return spectrum, spectrumP, spectrumIP,  phaseprofile, wvlBinEdges, phaseBinEdges, innerAp
            
            elif error is not False:

                return spectrum, phaseprofile, wvlBinEdges, phaseBinEdges, summed_arrayInErr, innerAp 
            
            else:
                return spectrum, phaseprofile, wvlBinEdges, phaseBinEdges, innerAp 
            
        elif error is not False:
    
            return spectrum, wvlBinEdges, summed_arrayInErr, innerAp 
        else:
            return spectrum, wvlBinEdges, innerAp




    @staticmethod
    def makeWvlBins(energyBinWidth=.1, wvlStart=3000, wvlStop=11000):
        """
        returns an array of wavlength bin edges, with a fixed energy bin width
        withing the limits given in wvlStart and wvlStop
        Args:
            energyBinWidth: bin width in eV
            wvlStart: Lower wavelength edge in Angstrom
            wvlStop: Upper wavelength edge in Angstrom
        Returns:
            an array of wavelength bin edges that can be used with numpy.histogram(bins=wvlBinEdges)
        
        Neil 8/1/14:
        Adding this in so thta wavelength bins can be found for wavelengths in virtual iomage photonlists
        should work fine.....
        """

        #Calculate upper and lower energy limits from wavelengths
        #Note that start and stop switch when going to energy
        energyStop = RADecImage.h * RADecImage.c * RADecImage.angstromPerMeter / wvlStart
        energyStart = RADecImage.h * RADecImage.c * RADecImage.angstromPerMeter / wvlStop
        nWvlBins = int((energyStop - energyStart) / energyBinWidth)
        #Construct energy bin edges
        energyBins = np.linspace(energyStart, energyStop, nWvlBins + 1)
        #Convert back to wavelength and reverse the order to get increasing wavelengths
        wvlBinEdges = np.array(RADecImage.h * RADecImage.c * RADecImage.angstromPerMeter / energyBins)
        wvlBinEdges = wvlBinEdges[::-1]
        return wvlBinEdges
        
    @staticmethod
    def makePhaseBins(binwidth = .02, pStart = 0, pStop = 1):
        nphaseBins = int((pStop - pStart)/binwidth)
        phaseBinEdges = np.linspace(pStart, pStop, nphaseBins + 1)
        
        return phaseBinEdges

def test(photListFileName='/Users/vaneyken/Data/UCSB/ARCONS/Palomar2012/corot18/testPhotonList-blosc.h5',
         vPlateScale=0.1, integrationTime=-1,firstSec=0):
    photList = tables.openFile(photListFileName,mode='r')
    try:
        im = RADecImage(photList,vPlateScale=vPlateScale,firstSec=firstSec,
                        integrationTime=integrationTime)
    finally:
        print 'Closing phot. list file.'
        photList.close()
    im.display()
    return im
        
