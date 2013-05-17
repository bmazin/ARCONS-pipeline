'''
Author: Julian van Eyken                    Date May 15 2013

Package/class for handling of images created from photon lists that are derotated and
mapped to sky coordinates.

Under construction....
'''
from util import utils
import numpy as np
import tables

class RADecImage(object): 
    '''
    Class to hold derotated, integrated and possibly stacked images, in the sky coordinate
    frame.
    '''
    
    def __init__(self,photList=None,nPixRA=None,nPixDec=None,cenRA=None,cenDec=None,
                 vPlateScale=0.1):
        '''
        Initialise an empty RA-dec coordinate frame image.
        vPlateScale in arcsec per pix
        '''
        self.nPixRA = nPixRA                    #No. of virtual pixels in RA direction
        self.nPixDec = nPixDec                  #No. of virtual pixels in dec. direction
        self.cenRA = cenRA                      #RA location of center of field (radians)
        self.cenDec = cenDec                    #Dec location of center of field (rad)
        self.vPlateScale = vPlateScale*2*np.pi/1296000          #No. of radians on sky per virtual pixel.
        if nPixRA is not None and nPixDec is not None:
            self.image = np.empty((self.nPixDec,self.nPixRA),dtype=float)
            self.image.fill(np.nan)
            self.expTimeWeights = np.empty((self.nPixDec,self.nPixRA),dtype=float)
            self.expTimeWeights.fill(np.nan)
            #self.gridRA = np.empty((self.nPixDec,self.nPixRA),dtype=float)
            #self.gridRA.fill(np.nan)
            #self.gridDec = np.empty((self.nPixRA,self.nPixDec),dtype=float)
            #self.gridDec.fill(np.nan)
        else:
            self.image = None
            self.expTimeWeights = None
            self.gridRA = None
            self.gridDec = None
        #if cenRA is not None and cenDec is not None and vPixPerArcsec is not None:
        #    self.setCoordGrid(cenRA,cenDec,vPlateScale)            
        if photList is not None:
            self.loadImage(photList)     
    
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
                  vPlateScale=None):
        '''
        Build a de-rotated stacked image from a photon list (PhotList) object
        
        INPUTS:
            photList - a PhotList object from which to construct the image.
            firstSec - time from start of exposure to start the 'integration' for the image (seconds)
            integrationTime - duration of integration time to include in the image (in seconds; -1 => to end of exposure)
            wvlMin, wvlMax - min and max wavelengths of photons to include in the image (Angstroms).
        
        '''
        
        posErr = 0.8    #Approx. position error in arcsec (just a fixed estimate for now, will improve later)
        posErr *= 2*np.pi/(60.*60.*360.)  #Convert to radians
        
        photTable = photList.root.photons.photons   #Shortcut to table
        
        #Get RA/dec range:
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
        photRAs = photTable.col('ra')       #Read all photon coords into an RA and a dec array.
        photDecs = photTable.col('dec')
        nPhot = len(photRAs)
        
        #For now, calculate and add uniform random dither to each photon (uniform circular distribution)
        ditherDists = np.random.rand(nPhot)*posErr
        ditherAngles = np.random.rand(nPhot)*2*np.pi
        ditherRAs=ditherDists*np.cos(ditherAngles)
        ditherDecs=ditherDists*np.sin(ditherAngles)
        
        photRAs=photRAs+ditherRAs
        photDecs=photDecs+ditherDecs
        
        print 'Making image'
        self.image, self.gridRA, self.gridDec = np.histogram2d(photRAs-self.cenRA,photDecs-self.cenDec,[self.gridRA-self.cenRA,self.gridDec-self.cenDec])
        print 'Done.'

    def display(self,normMin=None,normMax=None):
        utils.plotArray(self.image,cbar=True,normMin=normMin,normMax=normMax)


def test(photListFileName='/Users/vaneyken/Data/UCSB/ARCONS/Palomar2012/corot18/testPhotonList-blosc.h5',
         vPlateScale=0.1):
    photList = tables.openFile(photListFileName,mode='r')
    try:
        im = RADecImage(photList,vPlateScale=vPlateScale)
    finally:
        print 'Closing phot. list file.'
        photList.close()
    im.display()
    return im
        