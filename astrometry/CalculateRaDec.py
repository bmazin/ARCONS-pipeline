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

    def __init__(self,centroidListFileName):
        # Check centroidListFileName
        if os.path.isabs(centroidListFileName) == True:
            fullCentroidListFileName = centroidListFileName
        else:
            scratchDir = os.getenv('INTERM_PATH')
            centroidDir = os.path.join(scratchDir,'centroidListFiles')
            fullCentroidListFileName = os.path.join(centroidDir,centroidListFileName)
        if (not os.path.exists(fullCentroidListFileName)):
            print 'centroid list file does not exist: ', fullCentroidListFileName
            return

        # Load centroid positions, center of rotation, hour angle, and RA/DEC data from centroidListFile.
        self.centroidListFile = tables.openFile(fullCentroidListFileName, mode='r')
        self.params = np.array(self.centroidListFile.root.centroidlist.params.read())
        self.xCenterOfRotation = float(self.params[0])
        self.yCenterOfRotation = float(self.params[1])
        self.centroidRightAscension = self.params[2]
        self.centroidDeclination = self.params[3]
        self.times = np.array(self.centroidListFile.root.centroidlist.times.read())
        self.hourAngles = np.array(self.centroidListFile.root.centroidlist.hourAngles.read())
        self.xCentroids = np.array(self.centroidListFile.root.centroidlist.xPositions.read())
        self.yCentroids = np.array(self.centroidListFile.root.centroidlist.yPositions.read())
        self.centroidFlags = np.array(self.centroidListFile.root.centroidlist.flags.read())
        self.centroidListFile.close()
        self.frames=len(self.times)

        # Set origin to center of rotation. Calculate x and y centroid position in this new coordinate system.
        self.xCentroidOffset = self.xCentroids - self.xCenterOfRotation
        self.yCentroidOffset = self.yCentroids - self.yCenterOfRotation 

        # Calculate rotation matrix for each frame
        self.rotationMatrix = np.zeros(((self.frames,2,2)))
        self.xCentroidRotated = np.zeros(self.frames)
        self.yCentroidRotated = np.zeros(self.frames)
        for iFrame in range(self.frames):
            self.rotationMatrix[iFrame][0][0] = self.rotationMatrix[iFrame][1][1] = np.cos(self.hourAngles[iFrame])
            self.rotationMatrix[iFrame][0][1] = -np.sin(self.hourAngles[iFrame])
            self.rotationMatrix[iFrame][1][0] = -self.rotationMatrix[iFrame][0][1]
            self.xCentroidRotated[iFrame] = self.xCentroidOffset[iFrame]*self.rotationMatrix[iFrame][0][0] + self.yCentroidOffset[iFrame]*self.rotationMatrix[iFrame][0][1]
            self.yCentroidRotated[iFrame] = self.xCentroidOffset[iFrame]*self.rotationMatrix[iFrame][1][0] + self.yCentroidOffset[iFrame]*self.rotationMatrix[iFrame][1][1]
      
    def getRaDec(self,timestamp,xPhotonPixel,yPhotonPixel):
        self.timestamp = np.array(timestamp)
        self.xPhotonPixel = np.array(xPhotonPixel)
        self.yPhotonPixel = np.array(yPhotonPixel)

        self.inputLength = len(self.timestamp)
                
        self.deltaTime = self.times[1]-self.times[0]
        self.binNumber = (self.timestamp/self.deltaTime).astype('int')

        self.photonHourAngle = np.zeros(self.inputLength)
        for i in range(self.inputLength):
            self.photonHourAngle[i] = self.hourAngles[self.binNumber[i]]

        # Set origin to center of rotation.  Calculate x and y photon position in this new coordinate system.
        self.xPhotonOffset = self.xPhotonPixel - self.xCenterOfRotation
        self.yPhotonOffset = self.yPhotonPixel - self.yCenterOfRotation

        # Rotate by the hour angle at the origin.  This will lead to DEC = y and RA = -x?
        self.xPhotonRotated = np.zeros(self.inputLength)
        self.yPhotonRotated = np.zeros(self.inputLength)
        for i in range(self.inputLength):
            self.xPhotonRotated[i] = self.xPhotonOffset[i]*self.rotationMatrix[self.binNumber[i]][0][0] + self.yPhotonOffset[i]*self.rotationMatrix[self.binNumber[i]][0][1]
            self.yPhotonRotated[i] = self.xPhotonOffset[i]*self.rotationMatrix[self.binNumber[i]][1][0] + self.yPhotonOffset[i]*self.rotationMatrix[self.binNumber[i]][1][1]

        # Use the centroid as the zero point for ra and dec offsets
        self.declinationOffset = CalculateRaDec.platescale*(self.yPhotonRotated - self.yCentroidRotated[self.binNumber])
        self.rightAscensionOffset = -CalculateRaDec.platescale*(self.xPhotonRotated - self.xCentroidRotated[self.binNumber])

        # Convert centroid positions in DD:MM:SS.S and HH:MM:SS.S format to radians.
        self.centroidDeclinationRadians = ephem.degrees(self.centroidDeclination).real
        self.centroidRightAscensionRadians = ephem.hours(self.centroidRightAscension).real        
        
        # Convert centroid position radians to arcseconds.
        self.centroidDeclinationArcseconds = self.centroidDeclinationRadians * CalculateRaDec.radiansToDegrees * 3600.0
        self.centroidRightAscensionArcseconds = self.centroidRightAscensionRadians * CalculateRaDec.radiansToDegrees * 3600.0
        
        # Add the photon arcsecond offset to the centroid offset.
        self.photonDeclinationArcseconds = self.centroidDeclinationArcseconds + self.declinationOffset
        self.photonRightAscensionArcseconds = self.centroidRightAscensionArcseconds + self.rightAscensionOffset
        
        # Convert the photon positions from arcseconds to radians
        self.photonDeclinationRadians = (self.photonDeclinationArcseconds / 3600.0) * CalculateRaDec.degreesToRadians
        self.photonRightAscensionRadians = (self.photonRightAscensionArcseconds / 3600.0) * CalculateRaDec.degreesToRadians

        # Return the right ascension and declination, in radians               
        return self.photonDeclinationRadians, self.photonRightAscensionRadians, self.photonHourAngle





        

    


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
    print raDecObject.getRaDec(timestamp=timestamp,xPhotonPixel=xPhotonPixel,yPhotonPixel=yPhotonPixel)
    print (time()-tic)

