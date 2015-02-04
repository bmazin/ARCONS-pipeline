import math
import numpy as np
from util import ObsFileSeq as ObsFileSeq
from util import utils
import pyfits
from util.popup import PopUp, plotArray
import matplotlib.pyplot as plt
import pickle
import astrometry.CentroidCalc as cc
#import obsfileViewerTest as ovt

# Define a set of observation files for this mosaic run
name='ObjectFinderDemoMosaic'
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

RA = 283.3961625 #degrees
Dec = 33.029175

fd = ofs.getFrameDict()

ofs.loadSpectralCubes()

rmi = open('ObjectFinderDemoMosaic-cubes.pkl', 'rb')
data = pickle.load(rmi)

iFrameList = range(len(ofs.frameObsInfos))

image_list = []

for iFrame in iFrameList:
    c = data[iFrame]['cube']
    c = np.sum(c, axis = 2)
    t = data[iFrame]['effIntTime']
    image = c/t
    nanspot = np.isnan(image)
    image[nanspot] = 0
    image_list.append(image)

#these are all the matching stars i found in the 66 frames
#frame_list = np.array([image_list[0], image_list[10], image_list[11], image_list[18], image_list[19], image_list[28], image_list[29], image_list[34], image_list[35], image_list[36], image_list[37], image_list[42], image_list[43], image_list[65]])

#these are the frames previously used in chris's example
frame_list = np.array([image_list[0], image_list[28], image_list[29], image_list[65]])

#degPerPix, theta, raArcsecPerSec = ObsFileSeq.ObjectFinder(image_list, fd, RA, Dec)
degPerPix, theta, raArcsecPerSec = ObsFileSeq.ObjectFinder(frame_list, fd, RA, Dec)

ofs.setRm(degPerPix,
          math.degrees((theta)),
          raArcsecPerSec,
        )

mosaic = ofs.makeMosaicImage(range(66))

plotArray(mosaic)
