#!/bin/python

'''
Author: Paul Szypryt		Date: May 8, 2013

Based on RaDecMap.py.  Takes a list of centroids and hour angles and outputs the ra and dec for
an individual photon.
'''

import numpy as np
import tables
import os
from util.FileName import FileName
import ephem
from time import time


class CalculateRaDec:

    degreesToRadians = np.pi/180.0
    radiansToDegrees = 180.0/np.pi
    platescale=0.44 # arcsec/pixel

    def __init__(self,centroidListFile):
        '''
        centroidListFile may either be a pathname (string) to a centroid list file, a PyTables instance
        of a centroid list file, or a PyTables instance of the actual centroid list node (table) within
        a centroid list file (as stored in PhotList instances).
        '''
        
        if type(centroidListFile) is str:
            # Assume input is a pathname string - check centroidListFileName
            if os.path.isabs(centroidListFile) == True:
                fullCentroidListFileName = centroidListFile
            else:
                scratchDir = os.getenv('INTERM_PATH')
                centroidDir = os.path.join(scratchDir,'centroidListFiles')
                fullCentroidListFileName = os.path.join(centroidDir,centroidListFileName)
            if (not os.path.exists(fullCentroidListFileName)):
                raise IOError, 'Centroid list file does not exist: '+fullCentroidListFileName
                #print 'Centroid list file does not exist: ', fullCentroidListFileName
                #return
            # Load centroid positions, center of rotation, hour angle, and RA/DEC data from centroidListFile.
            clFile = tables.openFile(fullCentroidListFileName, mode='r')
            centroidListNode = clFile.root.centroidlist
        elif type(centroidListFile) is tables.file.File:
            clFile = None
            centroidListNode = centroidListFile.root.centroidlist
        elif type(centroidListFile) is tables.group.Group:
            clFile = None
            centroidListNode = centroidListFile
        else:
            raise ValueError('Input parameter centroidListFile must be either a pathname string or a PyTables file instance.')
        
        
        self.params = np.array(centroidListNode.params.read())
        '''
        self.xCenterOfRotation = float(self.params[0])
        self.yCenterOfRotation = float(self.params[1])
        '''
        self.centroidRightAscension = self.params[2]
        self.centroidDeclination = self.params[3]
        self.times = np.array(centroidListNode.times.read())
        self.hourAngles = np.array(centroidListNode.hourAngles.read())
        self.xCentroids = np.array(centroidListNode.xPositions.read())
        self.yCentroids = np.array(centroidListNode.yPositions.read())
        self.centroidFlags = np.array(centroidListNode.flags.read())
        if clFile is not None: clFile.close()
        self.frames=len(self.times)
        
        '''
        # Set origin to center of rotation. Calculate x and y centroid position in this new coordinate system.
        self.xCentroidOffset = self.xCentroids - self.xCenterOfRotation
        self.yCentroidOffset = self.yCentroids - self.yCenterOfRotation
        '''

        # Account for 0.5 offset in centroiding algorithm
        self.xCentroids -= 0.5
        self.yCentroids -= 0.5
        
        #Note - changed rotation direction here - JvE 6/5/2013
        self.rotationMatrix = np.array([[np.cos(self.hourAngles),-np.sin(self.hourAngles)],[np.sin(self.hourAngles),np.cos(self.hourAngles)]]).T  
        self.centroidRotated = np.dot(self.rotationMatrix,np.array([self.xCentroids,self.yCentroids])).diagonal(axis1=0,axis2=2)
      
        self.pixelCount = 44*46
        self.values = np.zeros((2,self.pixelCount))
        for i in range(self.pixelCount):
            self.values[0,i] = i/46
            self.values[1,i] = i%46
        self.rotatedValues = np.dot(self.rotationMatrix,self.values)
      
    def getRaDec(self,timestamp,xPhotonPixel,yPhotonPixel):
        self.timestamp = np.array(timestamp)
        self.xPhotonPixel = np.array(xPhotonPixel).astype('int')
        self.yPhotonPixel = np.array(yPhotonPixel).astype('int')

        self.inputLength = len(self.timestamp)
                
        self.deltaTime = self.times[1]-self.times[0]
        self.binNumber = np.array(self.timestamp/self.deltaTime).astype('int')

        self.photonHourAngle = self.hourAngles[self.binNumber]

        self.indexArray = np.array(self.xPhotonPixel*46+self.yPhotonPixel)

        self.xPhotonRotated=np.zeros(len(self.timestamp))
        self.yPhotonRotated=np.zeros(len(self.timestamp))
        
        # Find better way to do this, taking majority of the time currently
        for i in range(len(self.timestamp)):
            self.xPhotonRotated[i] = self.rotatedValues[self.binNumber[i]][0][self.indexArray[i]]
            self.yPhotonRotated[i] = self.rotatedValues[self.binNumber[i]][1][self.indexArray[i]]           

        # Use the centroid as the zero point for ra and dec offsets
        self.rightAscensionOffset = -CalculateRaDec.platescale*(self.xPhotonRotated - self.centroidRotated[0][self.binNumber])
        self.declinationOffset = CalculateRaDec.platescale*(self.yPhotonRotated - self.centroidRotated[1][self.binNumber])
        
        # Convert centroid positions in DD:MM:SS.S and HH:MM:SS.S format to radians.
        self.centroidRightAscensionRadians = ephem.hours(self.centroidRightAscension).real 
        self.centroidDeclinationRadians = ephem.degrees(self.centroidDeclination).real
                       
        # Convert centroid position radians to arcseconds.
        self.centroidDeclinationArcseconds = self.centroidDeclinationRadians * CalculateRaDec.radiansToDegrees * 3600.0
        self.centroidRightAscensionArcseconds = self.centroidRightAscensionRadians * CalculateRaDec.radiansToDegrees * 3600.0
        
        # Add the photon arcsecond offset to the centroid offset.
        self.photonDeclinationArcseconds = self.centroidDeclinationArcseconds + self.declinationOffset
        self.photonRightAscensionArcseconds = self.centroidRightAscensionArcseconds + self.rightAscensionOffset
        
        # Convert the photon positions from arcseconds to radians
        self.photonDeclinationRadians = (self.photonDeclinationArcseconds / 3600.0) * CalculateRaDec.degreesToRadians
        self.photonRightAscensionRadians = (self.photonRightAscensionArcseconds / 3600.0) * CalculateRaDec.degreesToRadians

        return self.photonRightAscensionRadians, self.photonDeclinationRadians, self.photonHourAngle


# Test script
if __name__ == "__main__":
    run = 'PAL2012'
    sunsetDate='20121208'
    utcDate='20121209'
    centroidTimestamp = '20121209-120530'
    calTimestamp = '20121209-131132'
    centroidListFileName=FileName(run=run,date=sunsetDate,tstamp=centroidTimestamp).centroidList()

    # Test photon
    xPhotonPixel=np.linspace(0,43, num =1000000).astype('int')
    yPhotonPixel=np.linspace(0,45,num = 1000000).astype('int')
    timestamp = np.linspace(0,299,num=1000000)
    
    tic = time()
    raDecObject = CalculateRaDec(centroidListFileName)
    outs = raDecObject.getRaDec(timestamp=timestamp,xPhotonPixel=xPhotonPixel,yPhotonPixel=yPhotonPixel)
    print time()-tic

    #print dec, ra

