# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:59:51 2015

@author: stoughto
"""
import math
import numpy as np
from util import ObsFileSeqV2 as ObsFileSeq
from util.readDict import readDict
import pyfits

d = readDict()
d.read_from_file("mosaic.dict")
name = 'PSN234416a'
dt = 200
ofs = ObsFileSeq.ObsFileSeq(name, d['run'], d['date'], d['timeStampList'], dt)

ofs.plotLocations("plotLocations.png")

# Pick two frames where one object exists
iFrames = [0, 45, 46, 55, 101]
iRows = [24, 16, 30, 12, 21]
iCols = [21, 13, 16, 28, 25]

# print out ra, dec, and time of the frames
for iFrame in range(ofs.nFrames):
    ra0 = ofs.tcsDict['raOffset'][iFrame]
    dec0 = ofs.tcsDict['decOffset'][iFrame]
    t0 = ofs.tcsDict['timeOffset'][iFrame]
    print "frame=%3d ra=%5.1f dec=%5.1f time=%9.4f"%(iFrame, ra0, dec0, t0)

# look for pairs of frames with identical ra,dec offsets
rowDot = 0
colDot = 0
for iFrame0 in iFrames:
    ra0 = ofs.tcsDict['raOffset'][iFrame0]
    dec0 = ofs.tcsDict['decOffset'][iFrame0]
    t0 = ofs.tcsDict['timeOffset'][iFrame0]
    for iFrame1 in iFrames:
        ra1 = ofs.tcsDict['raOffset'][iFrame1]
        dec1 = ofs.tcsDict['decOffset'][iFrame1]
        t1 = ofs.tcsDict['timeOffset'][iFrame1]
        if iFrame0 < iFrame1:
            if ra0 == ra1 and dec0 == dec1:
                print "====== iFrame0, iFrame1", iFrame0, iFrame1
                index0 = iFrames.index(iFrame0)
                row0, col0 = iRows[index0], iCols[index0]
                index1 = iFrames.index(iFrame1)
                row1, col1 = iRows[index1], iCols[index1]
                dr = float(row1 - row0)
                dc = float(col1 - col0)
                dt = t1 - t0
                rowDot = dr/dt
                colDot = dc/dt
                print "dt, rowDot, colDot", dt, rowDot, colDot
                #theta = math.atan2(dr, dc)
                #print 'theta=',math.degrees(theta)

thetas = dict(yes=[], no=[])
scales = dict(yes=[], no=[])

flipSign = dict(yes=-1, no=1)
for iFrame0 in iFrames:
    for iFrame1 in iFrames:
        if iFrame0 < iFrame1:
            print "\n=================", iFrame0, iFrame1
            ra0 = ofs.tcsDict['raOffset'][iFrame0]
            dec0 = ofs.tcsDict['decOffset'][iFrame0]
            dt0 = ofs.tcsDict['timeOffset'][iFrame0]
            ra1 = ofs.tcsDict['raOffset'][iFrame1]
            dec1 = ofs.tcsDict['decOffset'][iFrame1]
            dt1 = ofs.tcsDict['timeOffset'][iFrame1]
            index0 = iFrames.index(iFrame0)
            row0, col0 = iRows[index0]-dt0*rowDot, iCols[index0]-dt0*colDot
            index1 = iFrames.index(iFrame1)
            row1, col1 = iRows[index1]-dt1*rowDot, iCols[index1]-dt1*colDot
            if not (ra0 == ra1 and dec0 == dec1):
                matchList = [
                    {"ra": ra0, "dec": dec0, "row": row0, "col": col0},
                    {"ra": ra1, "dec": dec1, "row": row1, "col": col1}
                    ]
                for flip in ("yes", "no"):
                    scaleTheta = ObsFileSeq.ObsFileSeq.getScaleTheta(matchList, flip=flipSign[flip])
                    print "flip=%2d theta=%8.3f scale=%5.3f" % \
                    (flipSign[flip], math.degrees(scaleTheta['theta']), \
                        scaleTheta['scale'])
                    thetas[flip].append(scaleTheta['theta'])
                    scales[flip].append(scaleTheta['scale'])
            else:
                print "check how the drift correction worked"
                print "row0, row1, delta", row0, row1, row1-row0
                print "col0, col1, delta", col0, col1, col1-col0

bestThetaStd = 9999.99
for flip in ("yes", "no"):
    ta = np.asarray(thetas[flip])
    print "theta %3s %6.2f +/- %5.2f"%(flip, math.degrees(ta.mean()), math.degrees(ta.std()))
    sa = np.asarray(scales[flip])
    print "scale %3s %6.3f +/- %5.3f"%(flip, sa.mean(), sa.std())
    if ta.std() < bestThetaStd:
        bestThetaStd = ta.std()
        bestTheta = ta.mean()
        bestScale = sa.mean()
        bestFlip = flip

print "best flip, theta, scale=",bestFlip, bestTheta, bestScale

ofs.setTransform(bestScale, math.degrees(bestTheta), flipSign[bestFlip], rowDot, colDot)

iFrame0 = 0
iFrame1 = 101

for iFrame in iFrames:
    ofs.getR0C0(iFrame)

cps = ofs.makeMosaicImage([0], verbose=True)
pyfits.PrimaryHDU(cps).writeto('frame-000.fits', clobber=True)
cps = ofs.makeMosaicImage([101], verbose=True)
pyfits.PrimaryHDU(cps).writeto('frame-101.fits', clobber=True)
cps = ofs.makeMosaicImage([0, 101], verbose=True)
pyfits.PrimaryHDU(cps).writeto('frame-000-101.fits', clobber=True)
cps = ofs.makeMosaicImage(iFrames, verbose=True)
pyfits.PrimaryHDU(cps).writeto('frame-iFrames.fits', clobber=True)
cps = ofs.makeMosaicImage()
pyfits.PrimaryHDU(cps).writeto('frame-all.fits', clobber=True)



del ofs
