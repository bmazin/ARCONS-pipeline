'''
Author: Julian van Eyken                    Date: May 15 2013

Package/class for handling of images created from photon lists that are derotated and
mapped to sky coordinates.

Under construction....
'''
import numpy as np
import tables
import matplotlib.pyplot as mpl
import hotpix.hotPixels as hp
from util import utils
import photonlist.photlist as pl
from photonlist import boxer
from astrometry import CalculateRaDec as crd
from headers import pipelineFlags

class RADecImage(object): 
    '''
    Class to hold derotated, integrated and possibly stacked images, in the sky coordinate
    frame.
    '''
    
    def __init__(self,photList=None,nPixRA=None,nPixDec=None,cenRA=None,cenDec=None,
                 vPlateScale=0.1, detPlateScale=None, firstSec=0, integrationTime=-1,
                 expWeightTimeStep=1.0):
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
            expWeightTimeStep: float, time step to use when calculating exposure
                        time weights for the virtual pixels (seconds).
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
            self.expWeightTimeStep = expWeightTimeStep
        else:
            self.image = None
            self.effIntTimes = None
            self.expTimeWeights = None
            self.gridRA = None
            self.gridDec = None
            self.totExpTime = None
            self.expWeightTimeStep = expWeightTimeStep
        if (cenRA is not None and cenDec is not None and vPlateScale is not None
            and nPixRA is not None and nPixDec is not None):
            self.setCoordGrid()            
        if photList is not None:
            self.loadImage(photList,firstSec=firstSec,integrationTime=integrationTime)     
    
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
                  doStack=False,expWeightTimeStep=None, 
                  savePreStackImage=None, doWeighted=True):  #savePreStackImage is temporary for test purposes
        '''
        
        Build a de-rotated stacked image from a photon list (PhotList) object.
        If the RADecImage instance already contains an image, the new image is added to it.
        
        INPUTS:
            photList - a PhotList object from which to construct the image.
            firstSec - time from start of exposure to start the 'integration' for the image (seconds)
            integrationTime - duration of integration time to include in the image (in seconds; -1 or NaN => to end of exposure)
            wvlMin, wvlMax - min and max wavelengths of photons to include in the image (Angstroms).
            doStack - boolean; if True, then stack the image to be loaded on top of any image data already present.
            expWeightTimeStep - see __init__. If set here, overrides any value already set in the RADecImage object.
                                If the new image is being stacked on top of a current image, a new value can be
                                supplied that is different from the current image's value; but only the last value used
                                (i.e. the one supplied) will be stored in the class attribute.
            wvlMin, wvlMax - set min and max wavelength cutoffs for photons to be loaded in.
            savePreStackImage - temporary fudge, set to a file-name to save the image out to a file prior to stacking.
            doWeighted - if True, includes flatfield weighting factors from photons, and rejects photons from pixels where the
                                flatfield is bad at any wavelength within the requested wavelength range (all if wvlMin/wvl Max
                                not specified). 
        '''
        
        #posErr = 0.8    #Approx. position error in arcsec (just a fixed estimate for now, will improve later)
        #posErr *= 2*np.pi/(60.*60.*360.)  #Convert to radians
        
        photTable = photList.file.root.photons.photons   #Shortcut to table
        if expWeightTimeStep is not None:
            self.expWeightTimeStep=expWeightTimeStep
        
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
        deadPixMask = np.where(wvlCalFlagImage == pipelineFlags.waveCal['good'], 1.0, 0.0)   #1.0 where flag is good; 0.0 otherwise. (Straight boolean mask would work, but not guaranteed for Python 4....)

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
            photons = photTable.readWhere('(arrivalTime>=firstSec) & (arrivalTime<=lastSec)')
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
        toReject = np.in1d(allPhotXY,badXY)
        #Chuck out the bad photons
        print 'Rejecting photons from bad pixels...'
        photons = np.delete(photons,toReject)
        #########################################################################
        
        photRAs = photons['ra']       #Read all photon coords into an RA and a dec array.
        photDecs = photons['dec']
        photHAs = photons['ha']       #Along with hour angles...
        photWeights = photons['flatWeight']     #*photons['fluxWeight']
        photWavelengths = photons['wavelength']
        if wvlMin is not None or wvlMax is not None:
            assert all(photWavelengths>=wvlMin) and all(photWavelengths<=wvlMax)
        print 'Min, max photon wavelengths found: ', np.min(photWavelengths), np.max(photWavelengths)
        nPhot = len(photRAs)
        
        
        #Add uniform random dither to each photon, distributed over a square 
        #area of the same size and orientation as the originating pixel at 
        #the time of observation.
        xRand = np.random.rand(nPhot)*self.detPlateScale-self.detPlateScale/2.0
        yRand = np.random.rand(nPhot)*self.detPlateScale-self.detPlateScale/2.0       #Not the same array!
        ditherRAs = xRand*np.cos(photHAs) - yRand*np.sin(photHAs)
        ditherDecs = yRand*np.cos(photHAs) + xRand*np.sin(photHAs)
        
        photRAs=photRAs+ditherRAs
        photDecs=photDecs+ditherDecs
        
        #Make the image for this integration
        if doWeighted:
            print 'Making weighted image'
            thisImage,thisGridDec,thisGridRA = np.histogram2d(photDecs,photRAs,[self.gridDec,self.gridRA],
                                                                      weights=photWeights)
        else:
            print 'Making unweighted image'        
            thisImage,thisGridDec,thisGridRA = np.histogram2d(photDecs,photRAs,[self.gridDec,self.gridRA])
                
        if 1==1:
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
            tStartFramesAll = np.array(photList.file.root.centroidList.times.read()) #Convert to array, since it's saved as a list.
            tEndFramesAll = np.append(tStartFramesAll[1:], np.inf)                   #Last frame goes on forever as far as we know at the moment
            withinIntegration = ((tStartFramesAll < lastSec) & (tEndFramesAll > firstSec))
            tStartFrames = tStartFramesAll[withinIntegration].clip(min=firstSec)   #Now clip so that everything is within the requested integration time.
            tEndFrames = tEndFramesAll[withinIntegration].clip(max=lastSec)
            nFrames = len(tStartFrames)
            assert nFrames > 0      #Otherwise we have a problem....
            assert np.all(tStartFrames <= lastSec) and np.all(tEndFrames >= firstSec)
            
            #Get x,y locations of detector pixel corners (2D array of each x,y value, in detector space)
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
                detExpTimes = (hp.getEffIntTimeImage(photList.hotPixTimeMask, integrationTime=tEndFrames[iFrame]-tStartFrames[iFrame],
                                                     firstSec=tStartFrames[iFrame]) * deadPixMask * flatCalMask).flatten()
                
                #Loop over the virtual pixels and accumulate the exposure time that falls in each
                for iVDec in np.arange(nVPixDec):
                    for iVRA in np.arange(nVPixRA):
                        maybeOverlapping = np.where((dPixRANormMax > vPixRANormMin[iVRA]) & 
                                                 (dPixRANormMin < vPixRANormMax[iVRA]) &
                                                 (dPixDecNormMax > vPixDecNormMin[iVDec]) &
                                                 (dPixDecNormMin < vPixDecNormMax[iVDec]))[0]       #[0] since 'where' returns an array in a 1-element tuple
            
                        #Loop over the detector pixels which may be overlapping the current virtual pixel
                        for overlapLoc in maybeOverlapping:
                            #Calculate overlap fraction for given virtual pixel with given detector pixel
                            overlapFrac = boxer.boxer(iVDec,iVRA,dPixCornersDec[:,overlapLoc],dPixCornersRA[:,overlapLoc])      #Feed one set of pixel corners at a time to boxer.
                            #And add the contributing exposure time to vexptimes:
                            expTimeToAdd = overlapFrac*detExpTimes[overlapLoc]
                            vExpTimesStack[iVDec,iVRA,iFrame] += expTimeToAdd        #vExpTimes is 2D and it should all magically work out.
        
            #------------ End loop through time steps ----------
        
        #Sum up the exposure times from each frame:
        vExpTimes = np.sum(vExpTimesStack,axis=2)
        
        #------ Done with exposure time weighting, in principle ------
        
        
        #Temporary for testing-------------
        if savePreStackImage is not None:
            print 'Saving pre-stacked image to '+savePreStackImage
            mpl.imsave(fname=savePreStackImage,arr=thisImage,origin='lower',cmap=mpl.cm.gnuplot2,vmin=np.percentile(thisImage, 0.5), vmax=np.percentile(thisImage,99.5))
        #---------------------------------
        
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

        
        print 'Done.'



    def display(self,normMin=None,normMax=None,expWeight=True,pclip=None,colormap=mpl.cm.gnuplot2,
                image=None, logScale=False):
        '''
        Display the current image. Currently just a short-cut to utils.plotArray,
        but needs updating to mark RA and Dec on the axes.
        '''
        if expWeight:
            toDisplay = np.copy(self.image*self.expTimeWeights)
        else:
            toDisplay = np.copy(self.image)
        
        if logScale is True: toDisplay = np.log10(toDisplay)
        
        if image is not None: toDisplay = image
        
        if pclip:
            normMin = np.percentile(toDisplay[np.isfinite(toDisplay)],q=pclip)
            normMax = np.percentile(toDisplay[np.isfinite(toDisplay)],q=100.0-pclip)

        #Display NaNs as zeros so it looks better
        toDisplay[np.isnan(toDisplay)] = 0
        
        #Find the coordinates of the centers of the virtual pixels in degrees
        ra = (self.gridRA[0:-1] + self.gridRA[1:])/2.0 / np.pi * 180.
        dec = (self.gridDec[0:-1] + self.gridDec[1:])/2.0 / np.pi * 180.
        
        utils.plotArray(toDisplay,cbar=True,normMin=normMin,normMax=normMax,colormap=colormap)
        #mpl.imshow(toDisplay,vmin=normMin,vmax=normMax,extent=(ra[0],ra[-1],dec[0],dec[-1])) 

        


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
        
