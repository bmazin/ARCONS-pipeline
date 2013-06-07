'''
Author: Julian van Eyken                    Date May 15 2013

Package/class for handling of images created from photon lists that are derotated and
mapped to sky coordinates.

Under construction....
'''
from util import utils
from astrometry.CalculateRaDec import CalculateRaDec
import numpy as np
import tables
import matplotlib.pyplot as mpl

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
                        virtual image pixel)
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
        if detPlatescale is None:
            self.detPlatescale = CalculateRaDec.platescale*2*np.pi/1296000      #******For now - but this really needs reading in from the photon list file.
        else:
            self.detPlatescale = detPlatescale
        if nPixRA is not None and nPixDec is not None:
            self.image = np.empty((self.nPixDec,self.nPixRA),dtype=float)
            self.image.fill(np.nan)
            self.expTimeWeights = np.empty((self.nPixDec,self.nPixRA),dtype=float)
            self.expTimeWeights.fill(np.nan)
            self.gridRA = np.empty((self.nPixRA),dtype=float)
            self.gridRA.fill(np.nan)
            self.gridDec = np.empty((self.nPixDec),dtype=float)
            self.gridDec.fill(np.nan)
        else:
            self.image = None
            self.expTimeWeights = None
            self.gridRA = None
            self.gridDec = None
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
                  vPlateScale=None,stack=False,expWeightTimeStep=None,
                  savePreStackImage=None):  #savePreStackImage is temporary for test purposes
        '''
        Build a de-rotated stacked image from a photon list (PhotList) object.
        If the RADecImage instance already contains an image, the new image is added to it.
        
        INPUTS:
            photList - a PhotList object from which to construct the image.
            firstSec - time from start of exposure to start the 'integration' for the image (seconds)
            integrationTime - duration of integration time to include in the image (in seconds; -1 => to end of exposure)
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
        
        #Get RA/dec range:
        if integrationTime==-1:
            lastSec = np.inf
        else:
            lastSec = firstSec+integrationTime
       
        print 'Finding RA/dec ranges' 
        #Take advantage of the fact that the ra/dec columns are (or should be) indexed....
        self.raMin = photTable.cols.ra[photTable.colindexes['ra'][0]]
        self.raMax = photTable.cols.ra[photTable.colindexes['ra'][-1]]
        self.decMin = photTable.cols.dec[photTable.colindexes['dec'][0]]
        self.decMax = photTable.cols.dec[photTable.colindexes['dec'][-1]]
        self.cenRA = (self.raMin+self.raMax)/2.0
        self.cenDec = (self.decMin+self.decMax)/2.0
        
        #Set size of virtual grid to accommodate.
        if self.nPixRA is None:
            self.nPixRA = int((self.raMax-self.raMin)//self.vPlateScale + 2)     #+1 for round up; +1 because coordinates are the boundaries of the virtual pixels, not the centers.
        if self.nPixDec is None:
            self.nPixDec = int((self.decMax-self.decMin)//self.vPlateScale + 2)
        self.setCoordGrid()
                
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
        thisImage,thisGridRA,thisGridDec = np.histogram2d(photRAs-self.cenRA,photDecs-self.cenDec,
                                                          [self.gridRA-self.cenRA,self.gridDec-self.cenDec])
        
        if 1==0:
            #And now figure out the exposure time weights....
            tStartFrames = np.arange(start=firstSec,stop=lastSec,
                                     step=self.expWeightTimeStep)
            tEndFrames = (tStartFrames+self.expWeightTimeStep).clip(max=lastSec)    #Clip so that the last value doesn't go beyond the end of the exposure.
        
        #Temporary for testing-------------
        if savePreStackImage is not None:
            print 'Saving pre-stacked image to '+savePreStackImage
            mpl.imsave(fname=savePreStackImage,arr=thisImage)  #,vmin=np.percentile(thisImage, 0.5), vmax=np.percentile(thisImage,99.5))
        #---------------------------------
        
        if self.image is None or stack is False:
            self.image = thisImage
        else:
            print 'Stacking'
            self.image+=thisImage
        
        assert all(thisGridRA==self.gridRA-self.cenRA) and all(thisGridDec==self.gridDec-self.cenDec)
        
        print 'Done.'



    def display(self,normMin=None,normMax=None):
        utils.plotArray(self.image,cbar=True,normMin=normMin,normMax=normMax)


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
        