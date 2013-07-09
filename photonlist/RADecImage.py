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
from astrometry.CalculateRaDec import CalculateRaDec
import photonlist.boxer.boxer

class RADecImage(object): 
    '''
    Class to hold derotated, integrated and possibly stacked images, in the sky coordinate
    frame.
    '''
    
    def __init__(self,photList=None,nPixRA=None,nPixDec=None,cenRA=None,cenDec=None,
                 vPlateScale=0.1, detPlatescale=None, firstSec=0, integrationTime=-1,
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
        if detPlateScale is None:
            self.detPlatescale = CalculateRaDec.platescale*2*np.pi/1296000      #Radians per detector pixel. ******For now - but this really needs reading in from the photon list file.
        else:
            self.detPlatescale = detPlatescale
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
        else:
            self.image = None
            self.effIntTimes = None
            self.expTimeWeights = None
            self.gridRA = None
            self.gridDec = None
            self.totExpTime = None
        if cenRA is not None and cenDec is not None and vPixPerArcsec is not None:
            self.setCoordGrid(cenRA,cenDec,vPlateScale)            
        if photList is not None:
            self.loadImage(photList,firstSec=firstSec,integrationTime=integrationTime)     
    
    def setCoordGrid(self):
        '''
        #Establish RA and dec coordinates grid for virtual pixel grid 
        '''
        #self.gridRA = np.empty((self.nPixDec,self.nPixRA),dtype=float)
        #self.gridRA.fill(np.nan)
        #self.gridDec = np.empty((self.nPixDec,self.nPixRA),dtype=float)
        #self.gridDec.fill(np.nan)
        self.gridRA = self.cenRA + (self.vPlateScale*(np.arange(self.nPixRA) - (self.nPixRA//2)))
        self.gridDec = self.cenDec + (self.vPlateScale*(np.arange(self.nPixDec) - (self.nPixDec//2)))
    
    def loadImage(self,photList,firstSec=0,integrationTime=-1,wvlMin=-np.inf,wvlMax=np.inf,
                  stack=False,expWeightTimeStep=None,
                  savePreStackImage=None):  #savePreStackImage is temporary for test purposes
        '''
        Build a de-rotated stacked image from a photon list (PhotList) object.
        If the RADecImage instance already contains an image, the new image is added to it.
        
        INPUTS:
            photList - a PhotList object from which to construct the image.
            firstSec - time from start of exposure to start the 'integration' for the image (seconds)
            integrationTime - duration of integration time to include in the image (in seconds; -1 or NaN => to end of exposure)
            wvlMin, wvlMax - min and max wavelengths of photons to include in the image (Angstroms).
            stack - boolean; if True, then stack the image to be loaded on top of any image data already present.
            expWeightTimeStep - see __init__. If set here, overrides any value already set in the RADecImage object.
            savePreStackImage - temporary fudge, set to a file-name to save the image out to a file prior to stacking.
        '''
        
        #posErr = 0.8    #Approx. position error in arcsec (just a fixed estimate for now, will improve later)
        #posErr *= 2*np.pi/(60.*60.*360.)  #Convert to radians
        
        photTable = photList.file.root.photons.photons   #Shortcut to table
        if expWeightTimeStep is not None:
            self.expWeightTimeStep=expWeightTimeStep
        
        #Figure out last second of integration
        obsFileExpTime = photList.header.header.cols.exptime[0]
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
        
        print 'Getting photon coords'
        photons = photTable.readWhere('(arrivalTime>=firstSec) & (arrivalTime<=lastSec)')
        photRAs = photons['ra']       #Read all photon coords into an RA and a dec array.
        photDecs = photons['dec']
        photHAs = photons['ha']       #Along with hour angles...
        nPhot = len(photRAs)
        
        #Add uniform random dither to each photon, distributed over a square 
        #area of the same size and orientation as the originating pixel at 
        #the time of observation.
        xRand = np.random.rand(nPhot)*self.detPlatescale-self.detPlatescale/2.0
        yRand = np.random.rand(nPhot)*self.detPlatescale-self.detPlatescale/2.0       #Not the same array!
        ditherRAs = xRand*np.cos(photHAs) - yRand*np.sin(photHAs)
        ditherDecs = yRand*np.cos(photHAs) + xRand*np.sin(photHAs)
        
        photRAs=photRAs+ditherRAs
        photDecs=photDecs+ditherDecs
        
        #Make the image for this integration
        print 'Making image'
        thisImage,thisGridRA,thisGridDec = np.histogram2d(photRAs,photDecs,[self.gridRA,self.gridDec])
                
        if 1==0:
            #And now figure out the exposure time weights....
            
            print 'Calculating effective exposure times'
            tStartFrames = np.arange(start=firstSec,stop=lastSec,
                                     step=self.expWeightTimeStep)
            tEndFrames = (tStartFrames+self.expWeightTimeStep).clip(max=lastSec)    #Clip so that the last value doesn't go beyond the end of the exposure.
            nFrames = len(tStartFrames)
            
            #Get x,y locations of detector pixel corners (2D array of each x,y value, in detector space)
            dPixXmin = np.indices((nDPixRow,nDPixCol))[1] - 0.5
            dPixXmax = np.indices((nDPixRow,nDPixCol))[1] + 0.5
            dPixYmin = np.indices((nDPixRow,nDPixCol))[0] - 0.5
            dPixYmax = np.indices((nDPixRow,nDPixCol))[0] + 0.5
            
            #Create (1D) arrays for normalised center locations of virtual pixel grid (=index numbers, representing location of unit squares)
            vPixRANormCen = np.arange(nPixRA)   #np.indices(nVPixDec,nVPixRA)[1]
            vPixDecNormCen = np.arange(nPixDec) #np.indices(nVPixDec,nVPixRA)[0]
            
            #Create 1D arrays marking edges of virtual pixels (in 'normalised' space...)
            vPixRANormMin = np.arange(nVPixRA)-0.5
            vPixRANormMax = np.arange(nVPixRA)+0.5
            vPixDecNormMin = np.arange(nVPixDec)-0.5
            vPixDecNormMax = np.arange(nVPixDec)+0.5
            
            #Find origin of virtual array (center of virtual pixel 0,0) in RA/dec space.
            vPixOriginRA = np.mean(self.gridRA[0:1])     
            vPixOriginDec = np.mean(self.gridDec[0:1])
            vPixSize = self.vPlateScale       #Short hand, Length of side of virtual pixel in radians (assume square pixels)
            
            #Make array to take the total exposure times for each virtual pixel at each time step
            vExpTimes = np.zeros((nVPixDec,nVPixRA,nFrames))
            
            #Array to hold list of (equal) timestamps for each pixel at each timestep (just for calculating the RA/dec coordinates of the pixel corners)
            frameTime = np.zeros((nDPixRow,nDPixCol))
            frameTime.fill(np.nan)
            
            #------------ Loop through the time steps ----------
            for iFrame in range(nFrames):
                
                #Calculate detector pixel corner locations in RA/dec space    - NEEDS TO BE CLOCKWISE IN RA/DEC SPACE -- CHECK ON THIS!!
                frameTime.fill(tStartFrames[iFrame])
                dPixRA1,dPixDec1,dummy = raDecCalcObject.getRaDec(frameTime,dPixXmin,dPixYmin)   #May work with 2D arrays, need to find out....
                dPixRA2,dPixDec2,dummy = raDecCalcObject.getRaDec(frameTime,dPixXmin,dPixYmax)   #dPix* should all be 2D.
                dPixRA3,dPixDec3,dummy = raDecCalcObject.getRaDec(frameTime,dPixXmax,dPixYmax)
                dPixRA4,dPixDec4,dummy = raDecCalcObject.getRaDec(frameTime,dPixXmax,dPixYmin)
                
                #Normalise to scale where virtual pixel size=1 and origin is the origin of the virtual pixel grid
                dPixNormRA1 = (dPixRA1 - vPixOriginRA)/vPixSize     #dPixNorm* should all be 2D.
                dPixNormRA2 = (dPixRA2 - vPixOriginRA)/vPixSize
                dPixNormRA3 = (dPixRA3 - vPixOriginRA)/vPixSize
                dPixNormRA4 = (dPixRA4 - vPixOriginRA)/vPixSize
                dPixNormDec1 = (dPixDec1 - vPixOriginDec)/vPixSize
                dPixNormDec2 = (dPixDec2 - vPixOriginDec)/vPixSize
                dPixNormDec3 = (dPixDec3 - vPixOriginDec)/vPixSize
                dPixNormDec4 = (dPixDec4 - vPixOriginDec)/vPixSize
                    
                #Get min and max RA/decs for each of the detector pixels    
                dPixCornersRA = np.array([dPixNormRA1,dPixNormRA2,dPixNormRA3,dPixNormRA4])      #3D array, nRow x nCol x 4
                dPixCornersDec = np.array([dPixNormDec1,dPixNormDec2,dPixNormDec3,dPixNormDec4])
                dPixRANormMin = dPixCornersRA.min(axis=0)     #2D array, nRow x nCol
                dPixRANormMax = dPixCornersRA.max(axis=0)
                dPixDecNormMin = dPixCornersDec.min(axis=0)
                dPixDecNormMax = dPixCornersDec.max(axis=0)
                
                #Get array of effective exposure times for each detector pixel (2D array, nRow x nCol).
                detExpTimes = hp.getEffIntTimeImage(hotPixDict, integrationTime=tEndFrames[iFrame]-tStartFrames[iFrame],
                                                     firstSec-tStartFrames[iFrame])
        
                #Loop over the virtual pixels and accumulate the exposure time that falls in each
                for iVDec in arange(nVPixDec):
                    for iVRA in arange(nVPixRA):
                        maybeOverlapping = where(dPixRANormMax > vPixRANormMin[iVRA] & 
                                                 dPixRANormMin < vPixRANormMax[iVRA] &
                                                 dPixDecNormMax > vPixDecNormMin[iVDec] &
                                                 dPixDecNormMin < vPixDecNormMax[iVDec])
            
                        #Loop over the detector pixels which may be overlapping the current virtual pixel
                        for overlapLoc in maybeOverlapping:
                            #Calculate overlap fraction for given virtual pixel with given detector pixel
                            overlapFrac = photonlist.boxer.boxer(iVDec,iVRA,dPixCornersDec[:,:,overlapLoc],dPixCornersRA[:,:,overlapLoc])
                            #And add the contributing exposure time to vexptimes:
                            expTimeToAdd = overlapFrac*vdPixAreaRatio*detExpTimes*detExpTimes[overlapLoc]   #[overlapLoc] indexing *should* work...
                            vExpTimes[iVDec,iVRA,iFrame] += expTimeToAdd
        
            #------------ End loop through time steps ----------
        
        #------ Done with exposure time weighting, in principle ------
        
        
        #Temporary for testing-------------
        if savePreStackImage is not None:
            print 'Saving pre-stacked image to '+savePreStackImage
            mpl.imsave(fname=savePreStackImage,arr=thisImage,origin='lower',colormap=mpl.cm.gnuplot2)  #,vmin=np.percentile(thisImage, 0.5), vmax=np.percentile(thisImage,99.5))
        #---------------------------------
        
        if self.image is None or stack is False:
            self.image = thisImage
            self.effIntTimes = vExpTimes
            self.totExpTime = lastSec-firstSec
            self.expTimeWeights = self.totExpTime/self.EffIntTimes
        else:
            print 'Stacking'
            self.image += thisImage
            self.effIntTimes += vExpTimes
            self.totExpTime += lastSec-firstSec
            self.expTimeWeights = self.totExpTime/self.effIntTimes
        
        print 'Done.'



    def display(self,normMin=None,normMax=None,expWeight=True):
        '''
        Display the current image. Currently just a short-cut to utils.plotArray,
        but needs updating to mark RA and Dec on the axes.
        '''
        utils.plotArray(self.image*self.expTimeWeights,cbar=True,normMin=normMin,normMax=normMax)


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
        