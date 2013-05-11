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


class calculateRaDec:

    degreesToRadians = np.pi/180.0
    radiansToDegrees = 180.0/np.pi
    platescale=0.44 # arcsec/pixel

    def __init__(self,centroidListFileName,timestamp,xPhotonPixel,yPhotonPixel):
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

        # Calculate total number of frames, time step durations, and bin a given timestamp.
        self.frames=len(self.times)
        self.deltaTime = self.times[1]-self.times[0]
        self.binNumber = int(timestamp/self.deltaTime)

        # Set origin to center of rotation.  Calculate x and y photon position in this new coordinate system.
        self.xPhotonOffset = xPhotonPixel - self.xCenterOfRotation
        self.yPhotonOffset = yPhotonPixel - self.yCenterOfRotation

        # Calculate x and y centroid position in same new coordinate system.
        self.xCentroidOffset = self.xCentroids[self.binNumber] - self.xCenterOfRotation
        self.yCentroidOffset = self.yCentroids[self.binNumber] - self.yCenterOfRotation

        # Calculate rotation matrix
        self.rotationMatrix = np.zeros((2,2))
        self.rotationMatrix[0][0] = self.rotationMatrix[1][1] = np.cos(self.hourAngles[self.binNumber])
        self.rotationMatrix[0][1] = -np.sin(self.hourAngles[self.binNumber])
        self.rotationMatrix[1][0] = -self.rotationMatrix[0][1]

        # Rotate by the hour angle at the origin.  This will lead to DEC = y and RA = -x?
        self.xPhotonRotated = self.xPhotonOffset*self.rotationMatrix[0][0] + self.yPhotonOffset*self.rotationMatrix[0][1]
        self.yPhotonRotated = self.xPhotonOffset*self.rotationMatrix[1][0] + self.yPhotonOffset*self.rotationMatrix[1][1]
        self.xCentroidRotated = self.xCentroidOffset*self.rotationMatrix[0][0] + self.yCentroidOffset*self.rotationMatrix[0][1]
        self.yCentroidRotated = self.xCentroidOffset*self.rotationMatrix[1][0] + self.yCentroidOffset*self.rotationMatrix[1][1]

        # Use the centroid as the zero point for ra and dec offsets
        self.declinationOffset = calculateRaDec.platescale*(self.yPhotonRotated - self.yCentroidRotated)
        self.rightAscensionOffset = -calculateRaDec.platescale*(self.xPhotonRotated - self.xCentroidRotated)
        
        # Convert centroid positions in DD:MM:SS.S and HH:MM:SS.S format to radians.
        self.centroidDeclinationRadians = ephem.degrees(self.centroidDeclination).real
        self.centroidRightAscensionRadians = ephem.hours(self.centroidRightAscension).real        
        
        # Convert centroid position radians to arcseconds.
        self.centroidDeclinationArcseconds = self.centroidDeclinationRadians * calculateRaDec.radiansToDegrees * 3600.0
        self.centroidRightAscensionArcseconds = self.centroidRightAscensionRadians * calculateRaDec.radiansToDegrees * 3600.0
        
        # Add the photon arcsecond offset to the centroid offset.
        self.photonDeclinationArcseconds = self.centroidDeclinationArcseconds + self.declinationOffset
        self.photonRightAscensionArcseconds = self.centroidRightAscensionArcseconds + self.rightAscensionOffset
        
        # Convert the photon positions from arcseconds to radians
        self.photonDeclinationRadians = (self.photonDeclinationArcseconds / 3600.0) * calculateRaDec.degreesToRadians
        self.photonRightAscensionRadians = (self.photonRightAscensionArcseconds / 3600.0) * calculateRaDec.degreesToRadians

        # Print RA in more readable HH:MM:SS.S and DEC in DD:MM:SS.S format. Values not returned.
        self.photonDeclination = ephem.degrees(self.photonDeclinationRadians)
        self.photonRightAscension = ephem.hours(self.photonRightAscensionRadians)
        print 'Centroid RA: ' + self.centroidRightAscension + ', Centroid DEC: ' + self.centroidDeclination
        print 'Photon RA: ' + str(self.photonRightAscension) + ', Photon DEC: ' + str(self.photonDeclination)

        # Return the right ascension and declination, in radians
    def getRaDec(self):
        return self.photonDeclinationRadians, self.photonRightAscensionRadians

# Test Script, will eventually want to load this data with a params file
run = 'PAL2012'
sunsetDate='20121208'
utcDate='20121209'
centroidTimestamp = '20121209-120530'
calTimestamp = '20121209-131132'
centroidListFileName=FileName(run=run,date=sunsetDate,tstamp=centroidTimestamp).centroidList()

# Test photon
xPhotonPixel=25.0
yPhotonPixel=25.0
timestamp = 12.35223

print calculateRaDec(centroidListFileName,timestamp,xPhotonPixel,yPhotonPixel).getRaDec()

