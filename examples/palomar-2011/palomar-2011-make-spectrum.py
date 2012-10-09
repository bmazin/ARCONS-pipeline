#
# Look at a .h5 file from the Palomar 2011 run
# Make a FITS image of the photons
# Also, make a histogram of all the peak fits, which cna then be converted to an energy

# Use the /beammap/beamimage information to gt the list of roach board/pixel/time
#
# Set the environment variable MKID_DATA_DIR to point to the data location
#
# Example use:
# 
# $ export MKID_DATA_DIR=/Volumes/data/Palomar2011/Pal20110728
# python palomar-2011.py obs_20110729-151443.h5

import sys, os
import tables
import pyfits
import numpy as np
import matplotlib.pyplot as plt

def peakfit(y1,y2,y3):
	#returns the y value of the peak of the parabola defined by the given 3 points
	#where 3 points are separated by 1 in x
    y4=y2-0.125*((y3-y1)**2)/(y3+y1-2*y2)
    return y4

if (len(sys.argv) < 2):
    print "Usage:  ",sys.argv[0]," hdf5FileName"
    sys.exit(1)

# make the full file name by joining the input name to the MKID_DATA_DIR (or .)
hdf5FileName = sys.argv[1]
dataDir = os.getenv('MKID_DATA_DIR','.')
hdf5FullFileName = os.path.join(dataDir,hdf5FileName)
print "full file name is ",hdf5FullFileName
if (not os.path.exists(hdf5FullFileName)):
    print "file does not exist: ",hdf5FullFileName
    sys.exit(1)

# open the actual file.  This might take a while
fid = tables.openFile(hdf5FullFileName, mode='r')

# get the header information
header = fid.root.header.header
titles = header.colnames
info = header[0]

# get the beam image.  This is a 2d array of roach board/pixel/time locations
beamImage = fid.getNode("/beammap/beamimage")

# Make a 2d array to hold the image
shape = beamImage.shape
pixels = np.zeros(shape,dtype=np.uint32)

# Make a 1d array to hold the spectrum
hgFit = np.zeros(1000+2**12, dtype=np.uint32)
hgAt = np.zeros(1000+2**12, dtype=np.uint32)
hgTime = np.zeros(2**13, dtype=np.uint32)
pulseMask = int(12*'1',2) #bitmask of 12 ones
timeMask = int(20*'1',2)#bitmask of 20 ones

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
            timeVals = []
            for packet in sec:
                # here is the 64-bit number.  
                packet = int(packet)
                beforePeak = packet>>44 & pulseMask
                atPeak = packet>>32 & pulseMask
                afterPeak = packet>>20 & pulseMask
                peak = peakfit(beforePeak, atPeak, afterPeak)
                hgFit[peak] += 1
                hgAt[atPeak] += 1
                nPhoton += 1
                sum += 1
                time = packet & timeMask
                timeVals.append(time)
            timeVals.sort()
            

            previousTime = -9999
            for time in timeVals:
                if (previousTime != -9999) :
                    dt = time - previousTime
                    if (dt < hgTime.size-1):
                        hgTime[dt] += 1
                previousTime = time


        pixels[iRow][iCol] = sum
        print "iRow=",iRow,"  iCol=",iCol,"  sum=",sum
print "nPhoton=",nPhoton
hdu = pyfits.PrimaryHDU(pixels)
hdu.writeto('new.fits', clobber=True)

plt.clf()
plt.semilogy(hgFit, drawstyle='steps', label="fit peak")
plt.semilogy(hgAt, drawstyle='steps', label="at peak")
plt.xlabel("peak height")
plt.ylabel("dN")
plt.legend()
plt.title(hdf5FileName)
plt.savefig("peak-hg.png")

plt.clf()
plt.semilogy(hgTime, drawstyle='steps', label="within one pixel")
plt.xlabel("delta time (ADU)")
plt.ylabel("dN")
plt.legend()
plt.title(hdf5FileName)
plt.savefig("time-hg.png")
