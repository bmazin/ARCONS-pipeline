# Demonstrate how to use ObsFileSeq to make mosaic images
import math
import numpy as np
from util import ObsFileSeq as ObsFileSeq
from util import utils
import pyfits

# Define a set of observation files for this mosaic run
name='ring-20141020'
run = "PAL2014"
date = "20141020"
timeStampList = [
'20141021-033954',
'20141021-034532',
'20141021-035035',
'20141021-035538',
'20141021-040041',
'20141021-040544',
'20141021-041047',
]
dt = 200
ofs = ObsFileSeq.ObsFileSeq(name,run,date,timeStampList,dt)

# The ofs is now a sequence of frames.  A new frame starts
# whenever the telescope moves to a new location.
# If the telescope stays at one location for a time longer
# than two times dt, then that time is broken up to multiple frames,
# so dt is the longest exposure time for any one frame.


# Now make a plot of the ra,dec offset for each frame.
lfn = ofs.name+"-locations.png"
ofs.plotLocations(lfn)

# This is where the ObsFiles get calibrations applied.
# Look in util/ObsFileSeq.py to see all of the settings
# that are passed along to the calibration.

# This step seems to be really slow.  So, once it is done
# it writes out a pickle file with all of the spectral cubes
# in it.  The next time it is called, it simply reads back
# from that pickle file.

# WARNING -- if you change options to this call, (for example, 
# weighted) then delete or rename the pickle file!
ofs.loadSpectralCubes()

# Work out the plate scale (degrees per pixel)
# and how the row,col axes are rotated with respect to ra,dec axes

# Pick two frames where one object exists
iFrame0 = 28
iFrame1 = 29

# Define the row,col of the object in that frame.  Here, approximate
# by looking at images, to within 0.5 pixels.  
# It would be better to clean up bad pixels, run an object finding
# algoritm, and fit the position of the star.  But, for this first
# demonstration, I just looked at a couple of images and guessed
# the position to within a 1/2 pixel.
rc0 = np.array([10.,18.])
rc1 = np.array([9,36.5])
dCol = rc1[1]-rc0[1]
dRow = rc1[0]-rc0[0]
dPix = math.sqrt(((rc1-rc0)**2).sum())

# Get the ra,dec offsets (in arcsec) from the ObsFileSeq
dRa = ofs.tcsDict['raOffset'][iFrame1]-ofs.tcsDict['raOffset'][iFrame0]
dDec = ofs.tcsDict['decOffset'][iFrame1]-ofs.tcsDict['decOffset'][iFrame0]
dDeg = math.sqrt(dRa**2+dDec**2)/3600

# plate scale
degPerPix = dDeg/dPix
# rotation
thetaPix = math.atan2(dCol,dRow) # angle from vertical
thetaSky = math.atan2(dRa,dDec) # angle from North
theta = thetaPix-thetaSky

# The telecope is pretty accurate in dec, but
# it drifts a little bit in ra.  (The equatorial drive motor
# does not go at exactly the right speed.)

# Pick two frames where the ra,dec offset is zero,
# usually the beginning and ending frames
iFrameA = 0
iFrameB = 65
# Define the object in that frame.  Here, approximate
# by looking at images, to within 0.5 pixels.
rcA = np.array([22.,27.])
rcB = np.array([16.,27.])
sct = math.cos(theta)*degPerPix
sst = math.sin(theta)*degPerPix
# This rotation matrix converts from row,col to ra,dec in degrees
rm = np.array([[sct,-sst],[sst,sct]])

rdA = rm.dot(rcA)
rdB = rm.dot(rcB)
deltaRa = rdB[0]-rdA[0]
deltaTime = ofs.getTimeBySeq(iFrameB) - ofs.getTimeBySeq(iFrameA)
raArcsecPerSec = 3600*deltaRa/deltaTime

# Load the plate scale, rotation, and ra drift.
# Note that, in principle, we should not have to calculate
# these things for each mosaic.  The plate scale should not change
# during an entire run.  Theta should be the hour angle plus a
#constant offset.  The ra drift is probably constant, but it
# could depend on the hour angle or zenith distance.
# Once we calculate these at a number of different times
# during a run, we can use the general solution.  But, this
# gets us started.

ofs.setRm(degPerPix,
          math.degrees(theta),
          raArcsecPerSec,
        )

#
# Make a FITS file of each frame of the cube
for iFrame in range(66):
    frame = ofs.cubes[iFrame]['cube'].sum(axis=2)
    hdu = pyfits.PrimaryHDU(frame)
    fn = "%s-%02d.fit"%(ofs.name,iFrame)
    print "now make fn=",fn
    hdu.writeto(fn)

# Whew!  Now do the coaddition for all of the sequences.
# This uses all wavelengths.  The first improvement will be to
# have it use a subset of the wavelength bins.
mosaic = ofs.makeMosaicImage(range(66))

# Make a simple plot for now.  You can also save data as a FITS file, or combine
# the frames for three different wavelengths into a fabulous color picture!
# But right now, let's just dump out a heat map so we something to show off.
utils.plotArray(mosaic,cbar=True,plotTitle=ofs.name,showMe=False,plotFileName=ofs.name+"-all.png")

# Write it out as a FITS file, too
hdu = pyfits.PrimaryHDU(mosaic)
fn = "%s-all.fit"%ofs.name
hdu.writeto(fn)
del ofs
