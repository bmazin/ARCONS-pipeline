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


class calculateRaDec:

    degreesToRadians = np.pi/180.0
    radiansToDegrees = 180.0/np.pi
    platescale=0.44 # arcsec/pixel

    def __init__(self,centroidListFileName,xPhotonPixel,yPhotonPixel,xCenterOfRotation,yCenterOfRotation):
        if os.path.isabs(centroidListFileName) == True:
            fullCentroidListFileName = centroidListFileName
        else:
            scratchDir = os.getenv('INTERM_PATH')
            centroidDir = os.path.join(scratchDir,'centroidListFiles')
            fullCentroidListFileName = os.path.join(centroidDir,centroidListFileName)
        if (not os.path.exists(fullCentroidListFileName)):
            print 'centroid list file does not exist: ', fullCentroidListFileName
            return
        self.centroidListFile = tables.openFile(fullCentroidListFileName, mode='r')
        self.times = np.array(self.centroidListFile.root.centroidlist.times.read())
        self.hourAngles = calculateRaDec.degreesToRadians*np.array(self.centroidListFile.root.centroidlist.hourAngles.read())
        self.xCentroids = np.array(self.centroidListFile.root.centroidlist.xPositions.read())
        self.yCentroids = np.array(self.centroidListFile.root.centroidlist.yPositions.read())
        self.centroidFlags = np.array(self.centroidListFile.root.centroidlist.flags.read())
        self.frames=len(self.times)

        # Set origin to center of rotation.  Calculate x and y photon position in this new coordinate system.
        self.xPhotonOffset = xPhotonPixel - xCenterOfRotation
        self.yPhotonOffset = yPhotonPixel - yCenterOfRotation

        # Currently using data for frame 0, will eventually bin timestamp into a corresponding time range

        # Calculate x and y centroid position in same new coordinate system.
        self.xCentroidOffset = self.xCentroids[0] - xCenterOfRotation
        self.yCentroidOffset = self.yCentroids[0] - yCenterOfRotation

        # Calculate rotation matrix
        self.rotationMatrix = np.zeros((2,2))
        self.rotationMatrix[0][0] = self.rotationMatrix[1][1] = np.cos(self.hourAngles[0])
        self.rotationMatrix[0][1] = -np.sin(self.hourAngles[0])
        self.rotationMatrix[1][0] = -self.rotationMatrix[0][1]

        # Rotate by the hour angle at the origin.  This will lead to DEC = y and RA = -x?
        self.xPhotonRotated = self.xPhotonOffset*self.rotationMatrix[0][0] + self.yPhotonOffset*self.rotationMatrix[0][1]
        self.yPhotonRotated = self.xPhotonOffset*self.rotationMatrix[1][0] + self.yPhotonOffset*self.rotationMatrix[1][1]
        self.xCentroidRotated = self.xCentroidOffset*self.rotationMatrix[0][0] + self.yCentroidOffset*self.rotationMatrix[0][1]
        self.yCentroidRotated = self.xCentroidOffset*self.rotationMatrix[1][0] + self.yCentroidOffset*self.rotationMatrix[1][1]

        # Use the centroid as the zero point for ra and dec offsets
        self.declinationOffset = calculateRaDec.platescale*(self.yPhotonRotated - self.yCentroidRotated)
        self.rightAscensionOffset = -calculateRaDec.platescale*(self.xPhotonRotated - self.xCentroidRotated)

        print self.declinationOffset
        print self.rightAscensionOffset
        

# Test Script, will eventually want to load this data with a params file
run = 'PAL2012'
sunsetDate='20121208'
utcDate='20121209'
centroidTimestamp = '20121209-120530'
calTimestamp = '20121209-131132'
centroidListFileName=FileName(run=run,date=sunsetDate,tstamp=centroidTimestamp).centroidList()
xPhotonPixel=10.0
yPhotonPixel=10.0

calculateRaDec(centroidListFileName,xPhotonPixel,yPhotonPixel,xCenterOfRotation=30.0,yCenterOfRotation=30.0)


'''
def arcsec_to_radians(total_arcsec):
    total_degrees = total_arcsec/3600.0
    total_radians = total_degrees*d2r
    return total_radians

def radians_to_arcsec(total_radians):
    total_degrees = total_radians*r2d
    total_arcsec = total_degrees*3600.0
    return total_arcsec

# Create the h5 output file
out_file = 'coords_' + obs_file
h5out = openFile(path + out_file, mode = 'w')
ragroup = h5out.createGroup('/','ra', 'RA Map of Array')
decgroup = h5out.createGroup('/','dec', 'DEC Map of Array')
filt1 = Filters(complevel=0, complib='blosc', fletcher32=False)   

# Extract relevant header information from the h5 file
original_lst = h5file.root.header.header.col('lst')[0]
exptime = h5file.root.header.header.col('exptime')[0]
try:
    ts = h5file.root.header.header.col('unixtime')[0]
except KeyError:
    print 'Using "ut" instead of "unixtime" from header'
    ts = h5file.root.header.header.col('ut')[0]
print 'Original LST from telescope:', original_lst
print 'Exptime:', exptime

# Initial RA and LST
centroid_RA_radians = ephem.hours(centroid_RA).real
centroid_RA_arcsec = radians_to_arcsec(centroid_RA_radians)

centroid_DEC_radians = ephem.degrees(centroid_DEC).real
centroid_DEC_arcsec = radians_to_arcsec(centroid_DEC_radians)

original_lst_radians = ephem.hours(original_lst).real
original_lst_seconds = radians_to_arcsec(original_lst_radians)/15.0

rotation_matrix =np.zeros((2,2))
offsets_hypercube = np.zeros(((((grid_height,grid_width,2,exptime)))), dtype = '|S10')
# Calculations done for each second interval
for elapsed_time in range(exptime):
    # Load the image to find the star centroid        
    ra_array = h5out.createCArray(ragroup, 't%i' %elapsed_time, StringAtom(itemsize=10), (grid_height,grid_width), filters = filt1)
    dec_array = h5out.createCArray(decgroup, 't%i' %elapsed_time, StringAtom(itemsize=10), (grid_height,grid_width), filters = filt1)
    # Find an hour angle for each frame, assume center does not move
    current_lst_seconds = original_lst_seconds + elapsed_time
    current_lst_radians = arcsec_to_radians(current_lst_seconds*15.0)
    HA_variable = current_lst_radians - centroid_RA_radians
    HA_static = HA_offset*d2r
    HA_current = HA_variable + HA_static
    # Calculate rotation matrix elements
    rotation_matrix[0][0] = rotation_matrix[1][1] = np.cos(HA_current)
    rotation_matrix[0][1] = -np.sin(HA_current)
    rotation_matrix[1][0] = np.sin(HA_current)

    # Calculate the offsets from center
    x_offsets = np.zeros((grid_height,grid_width))
    y_offsets = np.zeros((grid_height,grid_width))
    rotated_x_offsets = np.zeros((grid_height,grid_width))
    rotated_y_offsets = np.zeros((grid_height,grid_width))
    for y in range(grid_height):
        for x in range(grid_width):
            # Unrotated matrix elements, multiplied by plate scale
            x_offsets[y][x] = plate_scale*(crpix1-x)
            y_offsets[y][x] = plate_scale*(crpix2-y)
            # Apply rotation by hour angle
            rotated_x_offsets[y][x] = centroid_RA_arcsec - plate_scale*(centroid_x[elapsed_time]-crpix1) + rotation_matrix[0][0]*x_offsets[y][x] + rotation_matrix[0][1]*y_offsets[y][x]
            rotated_y_offsets[y][x] = centroid_DEC_arcsec - plate_scale*(centroid_y[elapsed_time]-crpix2) + rotation_matrix[1][0]*x_offsets[y][x] + rotation_matrix[1][1]*y_offsets[y][x]   
            ra_array[y,x] = str(ephem.hours(arcsec_to_radians(rotated_x_offsets[y][x])))
            dec_array[y,x] = str(ephem.degrees(arcsec_to_radians(rotated_y_offsets[y][x])))
            h5out.flush()
'''
