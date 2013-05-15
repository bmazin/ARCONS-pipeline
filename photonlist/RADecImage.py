class RADecImage(object): 
    '''
    Class to hold derotated, integrated and possibly stacked images, in the sky coordinate
    frame.
    '''
    
    def __init__(self,nPixRA=200,nPixDec=200,cenRA=None,cenDec=None,
                 vPlateScale=None):
        '''
        Initialise an empty RA-dec coordinate frame image.
        '''
        self.nPixRA = nPixRA                    #No. of virtual pixels in RA direction
        self.nPixDec = nPixDec                  #No. of virtual pixels in dec. direction
        self.cenRA = cenRA                      #RA location of center of field (radians)
        self.cenDec = cenDec                    #Dec location of center of field (rad)
        self.vPlateScale = vPlateScale          #No. of radians on sky per virtual pixel.
        self.image = np.empty((self.nPixDec,self.nPixRA),dtype=float)
        self.image.fill(np.nan)
        self.expTimeWeights = np.empty((self.nPixDec,self.nPixRA),dtype=float)
        self.expTimeWeights.fill(np.nan)
        self.gridRA = np.empty((self.nPixDec,self.nPixRA),dtype=float)
        self.gridRA.fill(np.nan)
        self.gridDec = np.empty((self.nPixRA,self.nPixDec),dtype=float)
        self.gridDec.fill(np.nan)
        if cenRA is not None and cenDec is not None and vPixPerArcsec is not None:
            self.setCoordGrid(cenRA,cenDec,vPlateScale)            
             
    
    def setCoordGrid(self):
        '''
        #Establish RA and dec coordinates grid for virtual pixel grid 
        '''
        for iDec in range(nPixDec):
            self.gridRA[nPixDec,:] = self.cenRA + (vPlateScale*np.arange(nPixRA)-(nPixRA//2))
            self.gridDec[nPixDec,:] = self.cenDec + (vPlateScale*iDec*np.ones((nPixRA)) - (nPixDec//2))
    
    
    def loadImage(self,photList,firstSec=0,integrationTime=-1,wvlMin=-np.inf,wvlMax=np.inf,
                  vPlateScale=0.1):
        '''
        Build a de-rotated stacked image from a photon list (PhotList) object
        
        INPUTS:
            photList - a PhotList object from which to construct the image.
            firstSec - time from start of exposure to start the 'integration' for the image (seconds)
            integrationTime - duration of integration time to include in the image (in seconds; -1 => to end of exposure)
            wvlMin, wvlMax - min and max wavelengths of photons to include in the image (Angstroms).
        
        '''
        
        #Get RA/dec range:
        self.raMin = photTable.ra[photTable.colindexes['ra'][0]]
        self.raMax = photTable.ra[photTable.colindexes['ra'][-1]]
        self.decMin = photTable.dec[photTable.colindexes['dec'][0]]
        self.decMax = photTable.dec[photTable.colindexes['dec'][-1]]
        self.cenRA = (self.raMin+self.raMax)/2.0
        self.cenDec = (self.decMin+self.decMax)/2.0
        
        #Set size of virtual grid to accommodate.
        self.nPixRA = (raMax-raMin)/vPlateScale + 1     #+1 because coordinates are the boundaries of the virtual pixels, not the centers.
        self.nPixDec = (decMax-decMin)/vPlateScale + 1
        self.setCoordGrid()
        
        photRAs = photTable.col('ra')       #Read all photon coords into an RA and a dec array.
        photDecs = photTable.col('dec')
        
        self.image, self.gridRA, self.gridDec = histogram2d(photRAs,photDecs,[self.gridRA,self.gridDec])


    def displayImage(self):
        plotArray(self.image)
         
        