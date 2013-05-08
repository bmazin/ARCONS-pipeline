'''
Author: Julian van Eyken                Date: May 7 2013
For handling calibrated output photon lists.
'''

import numpy as np
import os.path
import tables

class PhotList:
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
        self.file = tables.openFile(self.fullFileName, mode='r')

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
        
    
    
    def getPhotonsForPixel(self):
        '''
        Space holder for now
        '''
        pass
    
    
    def getImage(self,firstSec=-np.Inf,integrationTime=np.Inf,wvlMin=-np.Inf,
                 wvlMax=np.Inf):
        '''
        Return a 2D image consisting of photon counts in each pixel.
        '''
        
        lastSec = firstSec+integrationTime
        
        #Initialise an empty image
        image = np.empty((self.nRow,self.nCol),dtype=float)
        image.fill(np.NaN)
        
        #Step through pixels and fill em up.
        for iCol in range(self.nCol):
            print iCol
            for iRow in range(self.nRow):
                #image[iRow,iCol] = len(self.photTable.getWhereList('(Ypix == iRow) & (Xpix == iCol) & (ArrivalTime >= firstSec) & (ArrivalTime < lastSec) & (Wavelength >= wvlMin) & (Wavelength < wvlMax)' ))
                image[iRow,iCol] = len(self.photTable.getWhereList('(Ypix == iRow)' ))
                
        #That should be it....
        return image
    
    