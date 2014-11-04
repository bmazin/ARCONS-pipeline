#
# Look at a .h5 file from the Palomar 2011 run
# Make a FITS image of the photons
# Use the /beammap/beamimage information to gt the list of roach board/pixel/time
#
# Set the environment variable MKID_RAW_PATH to point to the data location
#
# Example use:
# 
# $ export MKID_RAW_PATH=/Volumes/data/Palomar2011/Pal20110728
# python palomar-2011.py obs_20110729-151443.h5

import sys, os
import tables
import pyfits
import numpy as np

if (len(sys.argv) < 2):
    print "Usage:  ",sys.argv[0]," hdf5FileName"
    sys.exit(1)

# make the full file name by joining the input name to the MKID_RAW_PATH (or .)
hdf5FileName = sys.argv[1]
dataDir = os.getenv('MKID_RAW_PATH','.')
hdf5FullFileName = os.path.join(dataDir,hdf5FileName)
print "full file name is ",hdf5FullFileName
if (not os.path.exists(hdf5FullFileName)):
    print "file does not exist: ",hdf5FullFileName
    sys.exit(1)

# open the actual file.  This might take a while
fid = tables.openFile(hdf5FullFileName, mode='r')

# get the beam image.  This is a 2d array of roach board/pixel/time locations
beamImage = fid.getNode("/beammap/beamimage")

# Make a 2d array to hold the image
shape = beamImage.shape
pixels = np.zeros(shape,dtype=np.uint32)

# count the total number of photons in the file
nPhoton = 0
iRow = -1
for rows in beamImage:
    iRow += 1
    print "Begin iRow = ",iRow
    iCol = -1
    for pixel in rows:
        iCol += 1
        print "       iCol = ",iCol
        print "       pixel=",pixel
        # so now we have a roach board/pixel/time. 
        sum = 0
        for sec in fid.getNode(pixel):
            for packet in sec:
                # here is the 64-bit number.  
                packet = int(packet)
                nPhoton += 1
                sum += 1
        pixels[iRow][iCol] = sum
        print "iRow=",iRow,"  iCol=",iCol,"  sum=",sum
print "nPhoton=",nPhoton
hdu = pyfits.PrimaryHDU(pixels)
hdu.writeto('new.fits')
